# Segmentation
# Update date : Feb. 7, 2025

#' Segment Audio into Syllables
#'
#' @description
#' Segments audio recordings into syllables using dynamic thresholding of spectrograms.
#'
#' @param x An object to segment, either a file path or SAP object
#' @param start_time Start time in seconds
#' @param end_time End time in seconds
#' @param wl Window length for spectrogram (default: 256)
#' @param ovlp Overlap percentage (0-100) (default: 80)
#' @param flim Frequency limits in kHz (default: c(1, 10))
#' @param silence_threshold Threshold for silence detection (0-1) (default: 0.05)
#' @param min_syllable_ms Minimum syllable length in milliseconds (default: 50)
#' @param max_syllable_ms Maximum syllable length in milliseconds (default: 200)
#' @param min_level_db Minimum threshold level in dB (Default is 0)
#' @param max_level_db Maximum threshold level in dB (Default is 40)
#' @param db_delta Step size for threshold search in dB (default: 5)
#' @param search_direction Direction for threshold search: "up" or "down".
#'   \itemize{
#'     \item "up": Starts from min_level_db and increases (recommended for quiet or variable recordings)
#'     \item "down": Starts from max_level_db and decreases (recommended for loud, clear recordings)
#'   }
#' @param verbose Print progress messages (default: TRUE)
#' @param plot Display detection plot (default: TRUE)
#' @param smooth Smooth spectrogram visualization (default: FALSE)
#' @param save_plot Save detection plot (default: FALSE)
#' @param plot_dir Directory to save plots
#' @param segment_type For SAP objects: Type of segments ('bouts', 'motifs')
#' @param day For SAP objects: Days to process
#' @param indices For SAP objects: Specific indices to process
#' @param cores For SAP objects: Number of processing cores
#' @param plot_percent For SAP objects: Percentage of files to plot (default: 10)
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For WAV files:
#' \itemize{
#'   \item Reads and validates the audio file
#'   \item Computes spectrogram using Short-Time Fourier Transform
#'   \item Performs adaptive thresholding to detect syllables
#'   \item Validates detected segments against duration constraints
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Supports batch processing with parallel execution
#'   \item Processes specific days or indices
#'   \item Organizes results by source type
#'   \item Maintains metadata relationships
#' }
#'
#' dB Scale Conversion:
#' \itemize{
#'   \item User input: 0 to 40 dB (intuitive positive scale)
#'   \item Internal conversion: Subtracts reference level (20 dB)
#'   \item Actual dBFS: -60 to -20 dB relative to full scale
#' }
#'
#' Search Direction Guidelines:
#' \itemize{
#'   \item Use "up" when:
#'     \itemize{
#'       \item Recording has variable amplitude
#'       \item Background noise is significant
#'       \item Want to detect quieter syllables
#'     }
#'   \item Use "down" when:
#'     \itemize{
#'       \item Recording is clean with good SNR
#'       \item Want to avoid false positives
#'       \item Syllables are consistently loud
#'     }
#' }
#'
#' @return
#' Returns a data frame containing syllable information:
#' \itemize{
#'   \item filename: Name of the audio file
#'   \item selec: Sequential number for each syllable
#'   \item threshold: Final threshold used for detection
#'   \item .start: Start time relative to the analyzed segment
#'   \item .end: End time relative to the analyzed segment
#'   \item start_time: Absolute start time in the original audio file
#'   \item end_time: Absolute end time in the original audio file
#'   \item duration: Duration of syllable in seconds
#'   \item silence_gap: Gap to next syllable in seconds (NA for last syllable)
#' }
#'
#' If no syllables are detected, returns NULL.
#'
#' For SAP objects: Updated object with syllable information in segments slot
#'
#' @examples
#' \dontrun{
#' # Basic segmentation of WAV file
#' syllables <- segment("song.wav")
#'
#' # Custom parameters for clean recording
#' syllables <- segment("clean_song.wav",
#'                      search_direction = "down",
#'                      min_syllable_ms = 30,
#'                      max_syllable_ms = 150)
#'
#' # Process specific days in SAP object
#' sap_obj <- segment(sap_object,
#'                    segment_type = "bouts",
#'                    day = c(30, 40),
#'                    cores = 4)
#'
#' # Process with custom detection parameters
#' sap_obj <- segment(sap_object,
#'                    segment_type = "motifs",
#'                    min_level_db = 10,
#'                    max_level_db = 30,
#'                    save_plot = TRUE)
#' }
#'
#' @rdname segment
#' @export
segment <- function(x, ...) {
  UseMethod("segment")
}

#' @rdname segment
#' @export
segment.default <- function(x,  # x is wav file path
                            start_time = NULL,
                            end_time = NULL,
                            wl = 256,
                            ovlp = 80,
                            flim = c(1, 10),
                            silence_threshold = 0.05,
                            min_syllable_ms = 50,
                            max_syllable_ms = 200,
                            min_level_db = 0,
                            max_level_db = 40,
                            db_delta = 5,
                            search_direction = c("up", "down"),
                            verbose = TRUE,
                            plot = TRUE,
                            smooth = FALSE,
                            save_plot = FALSE,
                            plot_dir = NULL,
                            ...) {
  # Validate file path
  if (!file.exists(x)) {
    stop("File does not exist: ", x)
  }

  # Always read header first to get duration
  header <- tuneR::readWave(x, header = TRUE)
  duration <- header$samples / header$sample.rate

  # Handle start_time
  if (is.null(start_time)) {
    if (verbose) warning("'start_time' not provided - using 0s", call. = FALSE)
    start_time <- 0
  }

  # Handle end_time
  if (is.null(end_time)) {
    end_time <- duration
  } else {
    # Validate user-provided end_time against actual duration
    if (end_time > duration) {
      if (verbose) warning("end_time (", end_time, "s) exceeds file duration (",
                           round(duration, 2), "s) - using duration", call. = FALSE)
      end_time <- duration
    }
  }

  # Final boundary checks
  if (start_time < 0) {
    stop("start_time cannot be negative (got ", start_time, "s)")
  }

  if (start_time >= end_time) {
    stop("Invalid time window: start_time (", start_time,
         "s) >= end_time (", end_time, "s)")
  }

  # Validate window length and overlap
  if (!is.numeric(wl) || wl <= 0) {
    stop("Window length (wl) must be a positive number")
  }
  if (!is.numeric(ovlp) || ovlp < 0 || ovlp >= 100) {
    stop("Overlap percentage (ovlp) must be between 0 and 100")
  }

  # Validate flim parameter
  if (length(flim) != 2) {
    stop("flim must be a vector of length 2 (minimum and maximum frequency)")
  }

  # Set default plot directory if not provided
  if (save_plot && is.null(plot_dir)) {
    plot_dir <- file.path(dirname(x), "plots", "segment_detection")
  }

  # Validate search direction
  search_direction <- match.arg(search_direction)

  # Internal reference level (not exposed to users)
  ref_level_db <- 20

  # Pre-calculate constants
  min_syllable_length <- min_syllable_ms/1000
  max_syllable_length <- max_syllable_ms/1000
  filename <- basename(x)

  # Read wave file
  wave <- tuneR::readWave(x,
                          from = start_time,
                          to = end_time,
                          units = "seconds")

  # sampling_rate <- wave@samp.rate

  # Compute Short-Time Fourier Transform
  stft_result <- seewave::spectro(
    wave,
    wl = wl,
    ovlp = ovlp,
    fftw = TRUE,
    flimd = flim,
    plot = FALSE,
    osc = TRUE,
    cont = TRUE
  )

  times <- stft_result$time
  time_step <- mean(diff(times))
  frequencies <- stft_result$freq
  sp <- matrix(stft_result[["amp"]],
               nrow = length(frequencies),
               ncol = length(times))

  if(verbose) {
    cat(sprintf("Spectrogram summary:\nDimensions: %s\nRange: %.2f to %.2f\nTime range: %.2f to %.2f seconds\nFrequency range: %.2f to %.2f Hz\nTime step: %.2f ms\n",
                paste(dim(sp), collapse=" x "),
                range(sp)[1], range(sp)[2],
                range(times)[1], range(times)[2],
                range(frequencies)[1], range(frequencies)[2],
                time_step*1000))
  }

  # Initialize variables
  envelope_is_good <- FALSE
  final_envelope <- NULL
  final_threshold <- NULL
  final_regions <- NULL
  segments <- NULL

  # Calculate thresholds sequence based on direction
  thresholds <- if(search_direction == "up") {
    seq(min_level_db, max_level_db, by = db_delta)
  } else {
    seq(max_level_db, min_level_db, by = -db_delta)
  }

  if(verbose) {
    cat(sprintf("\nSearching thresholds %s: %.1f to %.1f dB\n",
                search_direction,
                thresholds[1],
                tail(thresholds, 1)))
  }

  # Dynamic threshold loop to find good envelope
  for(mldb in thresholds) {
    if(verbose) cat("\nTrying threshold:", mldb, "dB\n")

    # Normalize spectrogram
    spec_norm <- normalize_spec(sp, max_level_db = mldb, ref_level_db = ref_level_db)
    spec_median <- apply(spec_norm, 1, median)
    spec_norm <- sweep(spec_norm, 1, spec_median)  # More efficient than matrix replication
    spec_norm[spec_norm < 0] <- 0

    # Calculate vocal envelope more efficiently
    max_vals <- apply(spec_norm, 2, max)
    mean_vals <- colMeans(spec_norm)
    vocal_envelope <- max_vals * sqrt(mean_vals)

    if(all(vocal_envelope == 0)) {
      if(verbose) cat("Warning: All zero envelope\n")
      next
    }

    vocal_envelope <- vocal_envelope / max(vocal_envelope)

    # Find regions for vocalization
    vocal_regions <- find_continuous_regions(vocal_envelope > silence_threshold)

    if(nrow(vocal_regions) > 0) {
      # Calculate durations vectorized
      vocal_durations <- (vocal_regions[,2] - vocal_regions[,1]) * time_step
      max_vocal_len <- max(vocal_durations)

      if(verbose) {
        cat(sprintf("Longest vocalization: %.2f ms\n", max_vocal_len*1000))
      }

      if(max_vocal_len <= max_syllable_length) {
        envelope_is_good <- TRUE
        final_envelope <- vocal_envelope
        final_threshold <- mldb
        final_regions <- vocal_regions
        break
      }
    }
  }

  if(!envelope_is_good) {
    if(verbose) cat("Could not find suitable threshold\n")
    return(NULL)
  }

  if(verbose) cat("Found suitable threshold:", final_threshold, "dB\n")

  # Calculate durations vectorized
  region_durations <- times[final_regions[,2]] - times[final_regions[,1]]
  valid_idx <- region_durations >= min_syllable_length

  if(sum(valid_idx) > 0) {
    segments <- data.frame(
      filename = filename,
      selec = seq_len(sum(valid_idx)),
      threshold = final_threshold,
      .start = times[final_regions[valid_idx,1]],
      .end = times[final_regions[valid_idx,2]],
      start_time = start_time + times[final_regions[valid_idx,1]],
      end_time = start_time + times[final_regions[valid_idx,2]],
      duration = region_durations[valid_idx]
    )

    # Calculate silence gaps
    if(nrow(segments) > 1) {
      segments$silence_gap <- c(NA, segments$start[-1] - segments$end[-nrow(segments)])
    } else {
      segments$silence_gap <- NA
    }

    if(verbose) {
      cat(sprintf("\nFinal results:\nTotal segments found: %d\nDuration range: %.2f to %.2f ms\n",
                  nrow(segments),
                  min(segments$duration)*1000,
                  max(segments$duration)*1000))

      if(nrow(segments) > 1) {
        cat(sprintf("Silence gaps range: %.2f to %.2f ms\n",
                    min(segments$silence_gap[-1], na.rm=TRUE)*1000,
                    max(segments$silence_gap[-1], na.rm=TRUE)*1000))
      }
    }
  }
  # Plotting section
  if (!is.null(segments)) {
    plot_fn <- plot_syllable_detection(times = times,
                                       frequencies = frequencies,
                                       sp = sp,
                                       syllables = segments,
                                       max_level_db = max_level_db,
                                       ref_level_db = ref_level_db,
                                       final_threshold = final_threshold,
                                       final_envelope = final_envelope,
                                       silence_threshold = silence_threshold,
                                       smooth = smooth)

    if(plot) plot_fn()

    if(save_plot) {
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

      plot_filename <- file.path(plot_dir,
                                 paste0(tools::file_path_sans_ext(basename(x)),
                                        "_", round(start_time, 2),
                                        ".png"))

      png(filename = plot_filename, width = 1600, height = 1200, res = 300)
      plot_fn()
      dev.off()
    }
  } else {
    if (verbose) message("No segments detected - skipping plot")
  }

  return(segments)
}

#' @rdname segment
#' @export
segment.Sap <- function(x,  # x is SAP object
                        day = NULL,
                        indices = NULL,
                        segment_type = c("bouts", "motifs"),
                        cores = NULL,
                        save_plot = FALSE,
                        plot_percent = 10,
                        wl = 256,
                        ovlp = 80,
                        flim = c(1, 10),
                        silence_threshold = 0.05,
                        min_syllable_ms = 50,
                        max_syllable_ms = 200,
                        min_level_db = 0,
                        max_level_db = 40,
                        db_delta = 5,
                        search_direction = c("up", "down"),
                        verbose = TRUE,
                        ...) {

  if(verbose) message(sprintf("\n=== Starting Syllable Segmentation ==="))

  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  # Match arguments
  segment_type <- match.arg(segment_type)

  # Determine which data to process
  if (!is.null(segment_type)) {
    if (segment_type == "bouts" && !is.null(x$bouts)) {
      process_data <- x$bouts
    } else if (segment_type == "motifs" && !is.null(x$motifs)) {
      process_data <- x$motifs
    }
  } else {
    process_data <- x$metadata
  }

  # Filter data based on day
  if (!is.null(day)) {
    process_data <- process_data[process_data$day_post_hatch %in% day, ]
    days_to_process <- day

    if (nrow(process_data) == 0) {
      stop("No files found for specified day(s)")
    }
  } else {
    days_to_process <- unique(process_data$day_post_hatch)
  }

  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }

  # Create plot directories if save_plot is TRUE
  if (save_plot) {
    plots_dir <- file.path(x$base_path, "plots", "syllable_detection")
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Process each day
  all_results <- list()

  for (current_day in days_to_process) {
    # Get day-specific data
    day_data <- process_data[process_data$day_post_hatch == current_day, ]

    # Apply indices if specified
    if (!is.null(indices)) {
      valid_indices <- indices[indices <= nrow(day_data)]
      if (length(valid_indices) > 0) {
        day_data <- day_data[valid_indices, ]
      } else {
        cat(sprintf("\nNo valid indices for day %s.\n", current_day))
        next
      }
    }

    # Process ALL files
    file_indices <- 1:nrow(day_data)

    cat(sprintf("\nProcessing %d audio segments for day %s using %d cores.\n",
                length(file_indices), current_day, cores))

    # Determine which files to plot
    files_to_plot <- numeric(0)

    if (!is.null(indices) || save_plot) {
      if (!is.null(indices)) {
        files_to_plot <- file_indices
      } else if (save_plot) {
        n_plots <- ceiling(length(file_indices) * plot_percent / 100)
        files_to_plot <- sort(sample(file_indices, n_plots))
      }
    }

    # Define processing function
    process_segment <- function(i) {
      # Determine if this segment should be plotted
      should_plot_file <- i %in% files_to_plot

      # Set plot directory if needed
      plot_dir <- NULL
      if (should_plot_file && save_plot) {
        plot_dir <- file.path(plots_dir, paste0("day_", current_day))
        dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
      }

      # Construct file path
      wavfile <- construct_wav_path(day_data[i,], wav_dir = x$base_path)

      # Use segment.default with specific time window
      segments <- segment.default(
        x = wavfile,
        start_time = day_data$start_time[i],
        end_time = day_data$end_time[i],
        wl = wl,
        ovlp = ovlp,
        flim = flim,
        silence_threshold = silence_threshold,
        min_syllable_ms = min_syllable_ms,
        max_syllable_ms = max_syllable_ms,
        min_level_db = min_level_db,
        max_level_db = max_level_db,
        db_delta = db_delta,
        search_direction = search_direction,
        verbose = FALSE,
        plot = should_plot_file,
        save_plot = save_plot,
        plot_dir = plot_dir,
        ...
      )

      if (!is.null(segments) && nrow(segments) > 0) {
        current_source_type <- segment_type %||% "metadata"

        # Add base columns from day_data
        segments <- segments |>
          dplyr::mutate(
            day_post_hatch = day_data$day_post_hatch[i],
            label = day_data$label[i],
            .after = filename
          )

        # Handle selec merging only if source is bouts
        segments <- segments |>
          dplyr::mutate(
            source_type = current_source_type,
            source_selec = if ("selec" %in% names(day_data)) day_data$selec[i] else NA,
            selec = dplyr::case_when(
              current_source_type == "bouts" & !is.na(source_selec) ~
                paste(source_selec, selec, sep = "-"),
              TRUE ~ as.character(selec)
            )
          )
        return(segments)
      }
      return(NULL)
    }

    # Parallel processing
    if (cores > 1) {
      if (Sys.info()["sysname"] == "Linux") {
        day_results <- pbmcapply::pbmclapply(
          file_indices,
          process_segment,
          mc.cores = cores,
          mc.preschedule = FALSE
        )
      } else {
        day_results <- pbapply::pblapply(
          file_indices,
          process_segment,
          cl = cores
        )
      }
    } else {
      day_results <- pbapply::pblapply(
        file_indices,
        process_segment
      )
    }

    # Combine results for this day
    valid_detections <- day_results[!sapply(day_results, is.null)]
    if (length(valid_detections) > 0) {
      day_segments <- do.call(rbind, valid_detections)
      all_results[[as.character(current_day)]] <- day_segments

      cat(sprintf("\nProcessed segments in day %s. Total segments: %d\n",
                  current_day, nrow(day_segments)))
    } else {
      cat(sprintf("\nNo segments found in day %s.\n", current_day))
    }
  }

  # Combine all results and store in SAP object
  if (length(all_results) > 0) {
    final_results <- do.call(rbind, all_results)
    row.names(final_results) <- NULL

    # Store results in segments
    x$segments <- as_segment(final_results)

    # Update last modified time
    x$misc$last_modified <- Sys.time()

    message(sprintf("\nTotal segments detected across all days: %d", nrow(final_results)))
    message("Access results via: sap$segments")
    if (save_plot && !is.null(plots_dir)) {
      message("Review plots via: ", plots_dir)
    }
  } else {
    warning("No segments detected")
  }

  # Return the modified SAP object invisibly
  invisible(x)
}

#' Normalize Spectrogram
#'
#' @description
#' Internal function to normalize spectrogram using different methods.
#'
#' @param spec Input spectrogram matrix
#' @param max_level_db Maximum threshold level
#' @param ref_level_db Reference level
#' @param method Normalization method
#'
#' @details
#' Supports multiple normalization methods:
#' \itemize{
#'   \item db: dB-scale normalization
#'   \item minmax: Min-max scaling
#'   \item zscore: Z-score normalization
#' }
#'
#' @return
#' Normalized spectrogram matrix
#'
#' @keywords internal
normalize_spec <- function(spec,
                          max_level_db = NULL,
                          ref_level_db = NULL,
                          method = "db") {

  if (!is.matrix(spec)) {
    stop("Input must be a matrix")
  }

  switch(method,
         "db" = {
           # Convert to dB scale
           spec_db <- 20 * log10(abs(spec) + 1e-6)

           # Normalize with reference level
           spec_norm <- (spec_db - ref_level_db - max_level_db) / -max_level_db

           # Clip values
           spec_norm[spec_norm < 0] <- 0
           spec_norm[spec_norm > 1] <- 1

           return(spec_norm)
         },

         "minmax" = {
           # Simple min-max normalization
           spec_min <- min(spec)
           spec_max <- max(spec)
           spec_norm <- (spec - spec_min) / (spec_max - spec_min)
           return(spec_norm)
         },

         "zscore" = {
           # Z-score normalization
           spec_mean <- mean(spec)
           spec_sd <- sd(spec)
           spec_norm <- (spec - spec_mean) / spec_sd

           # Optional: rescale to 0-1
           spec_norm <- (spec_norm - min(spec_norm)) /
             (max(spec_norm) - min(spec_norm))
           return(spec_norm)
         },

         stop("Unknown normalization method")
  )
}

#' Find Continuous Regions
#'
#' @description
#' Internal function to identify continuous regions in a binary signal.
#'
#' @param signal Binary signal vector
#'
#' @details
#' Processes binary signal to:
#' \itemize{
#'   \item Find transitions between states
#'   \item Identify continuous regions
#'   \item Handle edge cases
#' }
#'
#' @return
#' Matrix of start and end indices
#'
#' @keywords internal
find_continuous_regions <- function(signal) {
  # Find all transitions
  transitions <- diff(c(0, signal, 0))
  # Get rising edges (0->1)
  starts <- which(transitions == 1)
  # Get falling edges (1->0)
  ends <- which(transitions == -1) - 1

  # Ensure valid regions (ends > starts)
  valid_regions <- ends > starts

  if(sum(valid_regions) == 0) {
    return(matrix(nrow=0, ncol=2))
  }

  regions <- cbind(starts[valid_regions], ends[valid_regions])

  # Remove first row if signal starts with 1
  if(signal[1] == 1 && nrow(regions) > 0) {
    regions <- regions[-1,, drop=FALSE]
  }

  # Remove last row if signal ends with 1
  if(signal[length(signal)] == 1 && nrow(regions) > 0) {
    regions <- regions[-nrow(regions),, drop=FALSE]
  }

  return(regions)
}

#' Plot Syllable Detection Results
#'
#' @description
#' Internal function to create visualization of syllable detection results.
#'
#' @param times Time points vector
#' @param frequencies Frequency points vector
#' @param sp Spectrogram matrix
#' @param syllables Detected syllables data frame
#' @param max_level_db Maximum threshold level
#' @param ref_level_db Reference level
#' @param final_threshold Final detection threshold
#' @param final_envelope Final detection envelope
#' @param silence_threshold Silence threshold
#' @param smooth Whether to smooth visualization
#'
#' @details
#' Creates visualization with:
#' \itemize{
#'   \item Spectrogram with syllable boundaries
#'   \item Detection envelope trace
#'   \item Threshold indicators
#'   \item Optional smoothing
#' }
#'
#' @return
#' Function that creates the plot
#'
#' @keywords internal
plot_syllable_detection <- function(times, frequencies, sp, syllables, max_level_db,
                                    ref_level_db, final_threshold, final_envelope, silence_threshold,
                                    smooth = FALSE) {
  # Create color palette
  mou_palette <- colorRampPalette(c("black", "red", "yellow", "white"))

  # Normalize spectrogram with max_level_db
  spec_norm <- normalize_spec(sp, max_level_db = max_level_db, ref_level_db = ref_level_db)
  spec_norm <- spec_norm - matrix(
    rep(apply(spec_norm, 1, median), ncol(spec_norm)),
    nrow = nrow(spec_norm)
  )
  spec_norm[spec_norm < 0] <- 0

  # If smooth is TRUE, interpolate for smoother visualization using base R
  if(smooth) {
    # Create finer grid
    times_plot <- seq(min(times), max(times), length.out = length(times) * 3)
    frequencies_plot <- seq(min(frequencies), max(frequencies), length.out = length(frequencies) * 3)

    # Initialize the smoothed matrix with correct dimensions
    spec_smooth <- matrix(0, nrow = length(frequencies_plot), ncol = length(times_plot))

    # First interpolate in time dimension
    temp_spec <- matrix(0, nrow = nrow(spec_norm), ncol = length(times_plot))
    for(i in 1:nrow(spec_norm)) {
      temp_spec[i,] <- approx(x = times, y = spec_norm[i,], xout = times_plot)$y
    }

    # Then interpolate in frequency dimension
    for(j in 1:ncol(temp_spec)) {
      spec_smooth[,j] <- approx(x = frequencies, y = temp_spec[,j], xout = frequencies_plot)$y
    }

    spec_plot <- t(spec_smooth)
  } else {
    times_plot <- times
    frequencies_plot <- frequencies
    spec_plot <- t(spec_norm)
  }

  # Function to create the actual plot
  create_plot <- function() {
    # Set up the plotting layout with 2 panels, detection trace on top
    layout(matrix(1:2, nrow=2), heights=c(1.5, 2))

    # Set consistent margins for both panels
    top_mar <- c(2,4,2,1)
    bottom_mar <- c(4,4,0,1)

    # Calculate x-axis limits
    xlim <- range(times)

    # Plot detection trace with black background
    par(mar=top_mar, bg='black')
    # Set plot region to match data coordinates exactly
    plot.new()
    plot.window(xlim=xlim, ylim=c(0, 1))
    lines(times, final_envelope, col='white')
    axis(1, col='white', col.axis='white')  # Add x-axis
    axis(2, col='white', col.axis='white')  # Add y-axis
    title(xlab="Time (s)", ylab="Amplitude", col.lab='white')
    abline(h=silence_threshold, col='red', lty=2)
    box(col='white')

    # Add legend for threshold
    legend("topright",
           legend=c(paste("Threshold =", round(silence_threshold, 3),
                   "(",final_threshold, "dB)")),
           col=c("red", "white"),
           lty=c(2, 0, 0, 0),  # only threshold line has a line style
           text.col="white",
           bg="transparent",
           box.col="transparent",
           cex=0.8,
           inset=0.02)

    # Plot spectrogram with black background
    par(mar=bottom_mar, bg='black')
    # Set plot region to match data coordinates exactly
    plot.new()
    plot.window(xlim=xlim, ylim=range(frequencies))

    # Use more colors and raster interpolation if smooth is TRUE
    n_colors <- if(smooth) 200 else 100
    image(times_plot, frequencies_plot, spec_plot,
          col=mou_palette(n_colors),
          add=TRUE,
          useRaster=smooth,
          interpolate=smooth)

    axis(1, col='white', col.axis='white')  # Add x-axis
    axis(2, col='white', col.axis='white')  # Add y-axis
    # Add titles and labels separately
    title(ylab="Frequency (kHz)", col.lab='white')
    title(xlab="Time (s)", col.lab='white')
    box(col='white')

    # Add syllable boxes and labels overlapping the spectrogram
    if(!is.null(syllables)) {
      # Get the plot dimensions in user coordinates
      usr <- par("usr")
      freq_range <- usr[4] - usr[3]

      # Define box parameters relative to frequency range
      box_height <- freq_range * 0.05
      box_top <- usr[4] - freq_range * 0.05
      box_bottom <- box_top - box_height
      label_y <- (box_top + box_bottom)/2

      for(i in 1:nrow(syllables)) {
        rect(syllables$.start[i], box_bottom,
             syllables$.end[i], box_top,
             col=rgb(1,0,0,0.6), border='white', lwd=2)

        text(x = (syllables$.start[i] + syllables$.end[i])/2,
             y = label_y,
             labels = i,
             col = 'white',
             cex = 0.8)
      }
    }

    # Reset layout
    layout(1)
  }

  # Return a function that creates the plot
  return(create_plot)
}

