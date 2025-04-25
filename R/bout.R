# Find Bout
# Update date : Feb. 7, 2025

#' Detect Song Bouts in Audio Recordings
#'
#' @description
#' Detects and analyzes song bouts using RMS amplitude thresholding with bandpass filtering.
#'
#' @param x An object to analyze, either a file path or SAP object
#' @param wl Window length for RMS calculation (default: 1024)
#' @param ovlp Overlap percentage between windows (default: 50)
#' @param norm_method Method for normalizing RMS values ("quantile" or "max")
#' @param rms_threshold Threshold for bout detection (default: 0.1)
#' @param min_duration Minimum bout duration in seconds (default: 0.5)
#' @param gap_duration Minimum gap between bouts (default: 0.3)
#' @param edge_window Time window for edge effects (default: 0.05)
#' @param freq_range Frequency range for bandpass filter (default: c(3, 5))
#' @param plot Whether to display visualization (default: TRUE)
#' @param save_plot Whether to save plots to file (default: FALSE)
#' @param plot_dir Directory for saving plots
#' @param day For SAP objects: Days to process
#' @param indices For SAP objects: Specific indices to process
#' @param segment_type For SAP objects: Type of segments (default: "motifs")
#' @param cores For SAP objects: Number of processing cores
#' @param plot_percent For SAP objects: Percentage of files to plot (default: 10)
#' @param summary For SAP objects: Include additional statistics (default: FALSE)
#' @param verbose For SAP objects: Print progress messages (default: TRUE)
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For WAV files:
#' \itemize{
#'   \item Applies bandpass filtering to focus on vocalization frequencies
#'   \item Calculates RMS envelope with specified window parameters
#'   \item Detects bouts using adaptive thresholding
#'   \item Handles edge cases and minimum duration constraints
#'   \item Creates optional visualizations
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Processes multiple recordings in parallel
#'   \item Validates bouts against existing motif detections
#'   \item Provides optional summary statistics
#'   \item Maintains metadata relationships
#'   \item Supports selective plotting
#' }
#'
#' When summary = TRUE for SAP objects with motifs:
#' \itemize{
#'   \item n_motifs: Count of motifs per bout
#'   \item align_time: First motif time for alignment
#'   \item bout_number_day: Sequential numbering
#'   \item bout_gap: Time from previous bout
#' }
#'
#' @return
#' For default method: Data frame containing:
#' \itemize{
#'   \item filename: Name of WAV file
#'   \item selec: Bout number
#'   \item start_time: Onset time
#'   \item end_time: Offset time
#' }
#'
#' For SAP objects: Updated object with bout information in bouts slot
#'
#' @examples
#' \dontrun{
#' # Basic bout detection from file
#' bouts <- find_bout("song.wav",
#'                    rms_threshold = 0.1,
#'                    min_duration = 0.7)
#'
#' # Custom parameters with visualization
#' bouts <- find_bout("song.wav",
#'                    freq_range = c(2, 8),
#'                    plot = TRUE,
#'                    save_plot = TRUE)
#'
#' # Process SAP object with summary
#' sap_obj <- find_bout(sap_object,
#'                      segment_type = "motifs",
#'                      day = c(30, 40),
#'                      summary = TRUE)
#'
#' # Process specific files with plots
#' sap_obj <- find_bout(sap_object,
#'                      indices = 1:5,
#'                      save_plot = TRUE,
#'                      cores = 4)
#' }
#'
#' @seealso \code{\link{segment}} for syllable-level segmentation
#'
#' @rdname find_bout
#' @export
find_bout <- function(x, ...) {
  UseMethod("find_bout")
}

#' @rdname find_bout
#' @export
find_bout.default <- function(x,  # x is wav file path
                              wl = 1024,
                              ovlp = 50,
                              norm_method = c("quantile", "max"),
                              rms_threshold = 0.1,
                              min_duration = 0.5,
                              gap_duration = 0.3,
                              edge_window = 0.05,
                              freq_range = c(3, 5),
                              plot = TRUE,
                              save_plot = FALSE,
                              plot_dir = NULL,
                              ...) {

  # Validate file path
  if (!file.exists(x)) {
    stop("File does not exist: ", x)
  }

  # Set default plot directory if not provided
  if (save_plot && is.null(plot_dir)) {
    plot_dir <- file.path(dirname(x), "plots", "bout_detection")
  }

  # Read entire wave file
  wv <- readWave(x)

  # Apply bandpass filter using seewave
  filtered_wave <- seewave::bwfilter(wave = wv,
                                     f = wv@samp.rate,
                                     n = 2,  # Filter order
                                     from = freq_range[1] * 1000,  # Convert kHz to Hz
                                     to = freq_range[2] * 1000,
                                     bandpass = TRUE)

  # Compute RMS envelope using filtered data
  rms_env <- compute_rms(filtered_wave, wl = wl, ovlp = ovlp)

  # Normalize RMS envelope
  norm_method <- match.arg(norm_method)

  if (norm_method == "quantile") {
    rms_env <-rms_env / quantile(rms_env, 0.95, na.rm = TRUE)
  } else { # "max"
    rms_env <-rms_env / max(rms_env, na.rm = TRUE)
  }

  # Calculate stride (hop size) in samples and time
  stride <- round(wl * (1 - ovlp / 100))
  stride_time <- stride / wv@samp.rate  # Convert stride to time

  # Calculate time points for RMS values
  num_windows <- length(rms_env)
  time_points <- seq(from = (wl / 2) / wv@samp.rate,
                     by = stride_time,
                     length.out = num_windows)

  # Threshold the RMS envelope
  above_thresh <- rms_env > rms_threshold

  # Handle edge effects at the start
  edge_window_samples <- ceiling(edge_window / stride_time)
  if (above_thresh[1]) {
    # Check if signal drops below threshold within the edge window
    edge_region <- above_thresh[1:min(edge_window_samples, length(above_thresh))]
    first_below <- which(!edge_region)[1]

    if (!is.na(first_below)) {
      # Set everything before the first drop to below threshold
      above_thresh[1:(first_below - 1)] <- FALSE
    }
  }

  # Recalculate crossings after edge adjustment
  crossings <- diff(above_thresh)

  # Find all onsets and offsets
  all_onsets <- which(crossings == 1) # Correct index after diff
  all_offsets <- which(crossings == -1) + 1L

  # # If signal starts above threshold, add onset at beginning
  # if (above_thresh[1]) {
  #   all_onsets <- c(1, all_onsets)
  # }

  # If signal ends above threshold, add offset at end
  if (above_thresh[length(above_thresh)]) {
    all_offsets <- c(all_offsets, length(rms_env))
  }

  # Convert durations from seconds to number of RMS samples
  gap_samples <- ceiling(gap_duration / stride_time)
  min_duration_samples <- ceiling(min_duration / stride_time)

  # Initialize vectors for bout detection
  bout_onsets <- numeric()
  bout_offsets <- numeric()

  # Detect bouts
  if (length(all_onsets) > 0 && length(all_offsets) > 0) {
    current_onset <- all_onsets[1]
    current_offset <- all_offsets[1]

    # Check if there are multiple onsets to iterate through
    if (length(all_onsets) >= 2) {
      for (idx in 2:length(all_onsets)) {
        # Calculate gap between current offset and next onset
        gap <- all_onsets[idx] - current_offset

        if (gap >= gap_samples) {
          # Found a gap, check if current bout meets minimum duration
          if ((current_offset - current_onset) >= min_duration_samples) {
            bout_onsets <- c(bout_onsets, current_onset)
            bout_offsets <- c(bout_offsets, current_offset)
          }
          # Start new bout
          current_onset <- all_onsets[idx]
          next_offset_idx <- which(all_offsets > current_onset)[1]
          if (!is.na(next_offset_idx)) {
            current_offset <- all_offsets[next_offset_idx]
          } else {
            current_offset <- length(rms_env)
          }
        } else {
          # Continue current bout
          next_offset_idx <- which(all_offsets > all_onsets[idx])[1]
          if (!is.na(next_offset_idx)) {
            current_offset <- all_offsets[next_offset_idx]
          } else {
            current_offset <- length(rms_env)  # Handle case with no subsequent offsets
          }
        }
      }
    }

    # Handle the last (or only) bout
    if ((current_offset - current_onset) >= min_duration_samples) {
      bout_onsets <- c(bout_onsets, current_onset)
      bout_offsets <- c(bout_offsets, current_offset)
    }
  }

  # Early return if no valid bouts
  if (length(bout_onsets) == 0 || length(bout_offsets) == 0) {
    return(NULL)
  }

  # Create bout data frame
  bout_df <- data.frame(
    filename = basename(x),
    selec = seq_along(bout_onsets),
    start_time = time_points[bout_onsets],
    end_time = time_points[bout_offsets],
    stringsAsFactors = FALSE
  )

  # Plotting section
  plot_fn <- plot_bout_detection(
    wav_file=x,
    time_points=time_points,
    rms_env=rms_env,
    rms_threshold=rms_threshold,
    bout_df=bout_df,
    bout_onsets=bout_onsets,
    bout_offsets=bout_offsets,
    wl=wl,
    ovlp=ovlp
  )

  if (plot) {
    plot_fn()
  }


  if(save_plot) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    plot_filename <- file.path(plot_dir,
                               paste0(tools::file_path_sans_ext(basename(x)),
                                      ".png"))

    png(filename = plot_filename, width = 1600, height = 1200, res = 300)
    plot_fn()
    dev.off()
  }

  return(bout_df)
}

#' @rdname find_bout
#' @export
find_bout.Sap <- function(x,  # x is SAP object
                          day = NULL,
                          indices = NULL,
                          segment_type = "motifs",
                          cores = NULL,
                          save_plot = FALSE,
                          plot_percent = 10,
                          wl = 1024,
                          ovlp = 50,
                          norm_method = c("quantile", "max"),
                          rms_threshold = 0.1,
                          min_duration = 0.5,
                          gap_duration = 0.3,
                          edge_window = 0.05,
                          freq_range = c(3, 5),
                          summary = FALSE,
                          verbose = TRUE,
                          ...) {

  if(verbose) message(sprintf("\n=== Starting Bout Detection ==="))

  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  # Determine which data to process
  if (!is.null(segment_type) && segment_type == "motifs") {
    if (is.null(x$motifs) || nrow(x$motifs) == 0) {
      stop("No motif data available in SAP object")
    }
    process_data <- x$motifs
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
    plots_dir <- file.path(x$base_path, "plots", "bout_detection")
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

    # Get unique files for the day
    unique_files <- which(!duplicated(day_data$filename))

    cat(sprintf("\nProcessing %d files for day %s using %d cores.\n",
                length(unique_files), current_day, cores))

    # Determine which files to plot
    files_to_plot <- numeric(0)

    if (!is.null(indices) || save_plot) {
      if (!is.null(indices)) {
        # Plot all specified files
        files_to_plot <- unique_files
      } else if (save_plot) {
        # Plot percentage of files
        n_plots <- ceiling(length(unique_files) * plot_percent / 100)
        files_to_plot <- sort(sample(unique_files, n_plots))
      }
    }

    # Define processing function for this day
    process_file <- function(i) {
      # Determine if this file should be plotted
      should_plot_file <- i %in% files_to_plot

      # Set plot directory if needed
      plot_dir <- NULL
      if (should_plot_file && save_plot) {
        plot_dir <- file.path(x$base_path, "plots", "bout_detection",
                              paste0("day_",current_day))
      }

      # Construct file path
      wavfile <- construct_wav_path(day_data[i,], wav_dir = x$base_path)
      current_filename <- basename(wavfile)

      # Use find_bout.default
      bout_df <- find_bout.default(
        x = wavfile,
        wl = wl,
        ovlp = ovlp,
        norm_method = norm_method,
        rms_threshold = rms_threshold,
        min_duration = min_duration,
        gap_duration = gap_duration,
        edge_window = edge_window,
        freq_range = freq_range,
        plot = should_plot_file,
        save_plot = save_plot,
        plot_dir = plot_dir,
        ...
      )

      # Process bout_df if it exists
      if (!is.null(bout_df) && nrow(bout_df) > 0) {
        # Validate against motifs ONLY if segment_type is "motifs"
        if (!is.null(segment_type) && segment_type == "motifs") {
          # Get corresponding motifs
          file_motifs <- x$motifs[x$motifs$filename == current_filename, ]

          if (nrow(file_motifs) > 0) {
            # Vectorized validation check
            valid_bouts <- sapply(1:nrow(bout_df), function(j) {
              any(file_motifs$detection_time >= bout_df$start_time[j] &
                    file_motifs$detection_time <= bout_df$end_time[j])
            })

            # Remove bouts without matching motifs
            bout_df <- bout_df[valid_bouts, ]
          }
        }

        if (nrow(bout_df) > 0) {
          bout_df <- bout_df |>
            dplyr::mutate(
              day_post_hatch = day_data$day_post_hatch[i],
              label = day_data$label[i],
              .after = filename
            )
          return(bout_df)
        }
      }
      return(NULL)

    }

    # Choose parallel processing method based on system and cores
    if (cores > 1) {
      if (Sys.info()["sysname"] == "Linux") {
        day_results <- pbmcapply::pbmclapply(
          unique_files,
          process_file,
          mc.cores = cores,
          mc.preschedule = FALSE
        )
      } else {
        day_results <- pbapply::pblapply(
          unique_files,
          process_file,
          cl = cores
        )
      }
    } else {
      day_results <- pbapply::pblapply(
        unique_files,
        process_file
      )
    }

    # Combine results for this day
    valid_detections <- day_results[!sapply(day_results, is.null)]
    if (length(valid_detections) > 0) {
      day_detections <- do.call(rbind, valid_detections)
      all_results[[as.character(current_day)]] <- day_detections

      cat(sprintf("\nProcessed files in day %s. Total bouts: %d\n",
                  current_day, nrow(day_detections)))
    } else {
      cat(sprintf("\nNo bouts found in day %s.\n", current_day))
    }
  }

  # Combine all results and store in SAP object
  if (length(all_results) > 0) {
    final_results <- do.call(rbind, all_results)
    row.names(final_results) <- NULL

    # Add summary information if requested
    if (summary && segment_type == "motifs") {
      # Get motif data once
      dt <- x$motifs$detection_time
      fn <- x$motifs$filename

      final_results <- final_results |>
        # Calculate motif information for each bout
        dplyr::group_by(filename, day_post_hatch, selec) |>
        dplyr::mutate(
          # Count motifs within bout
          n_motifs = {
            bout_motifs <- dt[dt >= start_time &
                                dt <= end_time &
                                fn == filename]
            length(bout_motifs)
          },
          # Get first motif time for alignment
          align_time = {
            bout_motifs <- dt[dt >= start_time &
                                dt <= end_time &
                                fn == filename]
            dplyr::first(bout_motifs)
          }
        ) |>
        # Add bout numbers by day
        dplyr::arrange(day_post_hatch, filename) |>
        dplyr::group_by(day_post_hatch) |>
        dplyr::mutate(bout_number_day = dplyr::row_number()) |>
        # Add bout gaps by file
        dplyr::group_by(filename) |>
        dplyr::mutate(
          bout_gap = dplyr::case_when(
            dplyr::n() == 1 ~ NA_real_,  # Single bout in file
            dplyr::row_number() == 1 ~ NA_real_,  # First bout in file
            TRUE ~ start_time - dplyr::lag(end_time)  # Gap from previous bout
          )
        ) |>
        dplyr::ungroup()
    }

    # Store results in bouts
    x$bouts <- as_segment(final_results)

    # Update last modified time
    x$misc$last_modified <- Sys.time()

    message(sprintf("\nTotal bouts detected across all days: %d", nrow(final_results)))
    message("Access results via: sap$bouts")
    if (save_plot && !is.null(plots_dir)) {
      message("Review plots via: ", plots_dir)
    }
  } else {
    warning("No bouts detected")
  }

  # Return the modified SAP object invisibly
  invisible(x)
}

#' Plot Bout Detection Results
#'
#' @description
#' Creates a two-panel plot showing RMS envelope and spectrogram with bout boundaries.
#'
#' @inheritParams find_bout.default
#' @keywords internal
plot_bout_detection <- function(wav_file, time_points, rms_env, rms_threshold,
                                bout_df, bout_onsets, bout_offsets, wl, ovlp) {

  create_plot <- function() {
    # # Reset the graphics device first
    # graphics.off()
    # # Open a new device if none is active
    # if (!dev.cur()) dev.new()

    # Set up the plotting layout with 2 panels
    old_par <- par(no.readonly = TRUE)  # Save old parameters
    on.exit(par(old_par))  # Restore parameters when function exits

    # Set up the plotting layout with 2 panels
    layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(1.5,2))

    # Set consistent margins for both panels
    top_mar <- c(0,4,2,2)
    bottom_mar <- c(4,4,0,2)

    # Calculate x-axis limits
    xlim <- range(time_points)
    ylim_rms <- c(0, max(rms_env) * 1.1)

    # Plot RMS envelope
    par(mar=top_mar, bg='black')
    # Set plot region to match data coordinates exactly
    plot.new()
    plot.window(xlim=xlim, ylim=ylim_rms)

    # Add RMS envelope line
    lines(time_points, rms_env, col="red")

    # Add threshold line
    abline(h=rms_threshold, col="green", lty=2)

    # Add bout markers
    if(nrow(bout_df) > 0) {
      for(b in 1:nrow(bout_df)) {
        onset_idx <- bout_onsets[b]
        offset_idx <- bout_offsets[b]

        points(bout_df$start_time[b], rms_env[onset_idx],
               col="purple", pch=19, cex=1)
        points(bout_df$end_time[b], rms_env[offset_idx],
               col="orange", pch=19, cex=1)
      }
    }

    # Add axes and labels
    axis(2, col="white", col.axis="white")  # y-axis only for top panel
    title(ylab="RMS", col.lab="white")
    box(col="white")

    # Add legend with white text
    legend("topright",
           legend=c("RMS", "Threshold", "Onset", "Offset"),
           col=c("red", "green", "purple", "orange"),
           lty=c(1, 2, NA, NA),
           pch=c(NA, NA, 19, 19),
           bty="n",
           text.col="white")


    # Plot spectrogram
    par(mar=bottom_mar)
    visualize_song.default(wav_file,
                           fft_window_size=wl,
                           overlap=ovlp/100,
                           dark_mode=TRUE,
                           legend=FALSE,
                           keep.par=TRUE,
                           verbose = FALSE)

    # Add bout boundary lines to spectrogram
    if(nrow(bout_df) > 0) {
      for(b in 1:nrow(bout_df)) {
        abline(v=bout_df$start_time[b], col="white", lty=2, lwd=1.5)
        abline(v=bout_df$end_time[b], col="white", lty=2, lwd=1.5)
      }
    }

    # Reset layout
    layout(1)
  }

  # Return the plotting function
  return(create_plot)
}

#' Compute RMS Envelope
#'
#' @description
#' Internal function to calculate RMS envelope of audio signal.
#'
#' @inheritParams find_bout.default
#' @keywords internal
compute_rms <- function(data, wl, ovlp) {
  # Ensure data is a numeric vector
  data <- as.numeric(data)

  # Calculate stride in samples
  stride <- round(wl * (1 - ovlp / 100))

  # Square the data
  data_squared <- data^2

  # Create the moving average filter (window)
  window <- rep(1 / wl, wl)

  # Perform convolution using the filter
  rms_squared <- stats::filter(data_squared, filter = window, sides = 1)

  # Remove NA values introduced by the filter function
  rms_squared <- rms_squared[wl:length(rms_squared)]

  # Calculate indices to sample based on the stride
  indices <- seq(1, length(rms_squared), by = stride)

  # Get the RMS values at the specified strides
  rms_values <- sqrt(rms_squared[indices])

  return(rms_values)
}
