
# Spectrograming from audio segments -------------------------------------------------------------------
# Update date : Apr. 21, 2025

##' Extract and pad spectrograms from audio segments
#'
#' @param x Data frame with audio segment information or a SAP object
#' @param wav_dir Directory containing WAV files (default: NULL)
#' @param cores Number of CPU cores to use for parallel processing (default: NULL, uses available cores - 1)
#' @param wl Window length for spectrogram analysis (default: 256 )
#'          Time resolution: wl/sampling rate per window
#'          Frequency resolution: sampling rate/wl Hz per bin
#' @param ovlp Overlap percentage between windows (default: 20 )
#' @param wn Window function name (default: "hanning")
#' @param freq_range Frequency range to analyze in kHz (default: c(1,10) )
#' @param fftw Logical, use FFTW or not (default: TRUE)
#' @param segment_type Type of segments to process: "segments" or "syllables" (default: "segments")
#' @param sample_percent Percentage of segments to sample (default: NULL)
#' @param balanced Whether to balance samples across labels (default: FALSE)
#' @param labels Specific labels to include (default: NULL)
#' @param seed Random seed for reproducible sampling (default: 222)
#' @param verbose Whether to display progress messages (default: TRUE)
#' @param ... Additional parameters passed to spectro
#'
#' @return For default method: data frame with original metadata and flattened, padded spectrograms
#'         For SAP method: updated SAP object with spectrograms added to features
#' @export
extract_spec <- function(x, ...) {
  UseMethod("extract_spec")
}

#' Default method for extract_spec
#' @rdname extract_spec
#' @export
extract_spec.default <- function(x,
                                 wav_dir = NULL,
                                 cores = NULL,
                                 wl = 256,  # Time resolution: 256/44100 ≈ 5.8 ms per window /Frequency resolution: 44100/256 ≈ 172.2 Hz per bin
                                 ovlp = 20,
                                 wn = "hanning",
                                 freq_range = c(1,10),
                                 fftw = TRUE,
                                 ...) {

  # Check required columns
  required_cols <- c("filename", "start_time", "end_time")
  missing_cols <- required_cols[!required_cols %in% names(x)]
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Check if data frame is empty
  if (nrow(x) == 0) {
    stop("Input data frame is empty")
  }

  # Handle wav_dir
  if (is.null(wav_dir)) {
    wav_dir <- attr(x, "wav_dir")
    if (is.null(wav_dir)) {
      stop("wav_dir must be provided either as argument or attribute")
    }
  }

  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
    cores <- max(1, cores)  # Ensure at least 1 core
  }

  cat(sprintf("\nExtracting spectrograms from %d audio segments using %d cores.\n",
              nrow(x), cores))

  # Function to process a single row
  process_row <- function(i) {
    extract_single_spec(x[i, ],
                        wav_dir = wav_dir,
                        wl = wl,
                        ovlp = ovlp,
                        wn = wn,
                        freq_range = freq_range,
                        fftw = fftw,
                        ...)
  }

  # Extract spectrograms in parallel
  if (fftw) ensure_pkgs("fftw")

  spectrograms <- parallel_apply(
    seq_len(nrow(x)),
    process_row,
    cores = cores
  )

  # Filter out NULL results
  valid_indices <- which(!sapply(spectrograms, is.null))
  valid_spectrograms <- spectrograms[valid_indices]

  if (length(valid_spectrograms) == 0) {
    warning("No valid spectrograms extracted")
    return(NULL)
  }

  # Find the maximum length for padding
  max_length <- max(sapply(valid_spectrograms, length))

  cat(sprintf("\nPadding spectrograms to match length: %d\n", max_length))

  # Pad all spectrograms to the same length
  padded_spectrograms <- lapply(valid_spectrograms, function(spec) {
    pad_array(spec, max_length)
  })

  # Create a data frame with the original metadata
  df <- x[valid_indices, ]

  # Check which metadata columns actually exist in the data
  metadata_cols <- c("filename", "day_post_hatch", "label", "start_time", "end_time")
  available_cols <- intersect(metadata_cols, names(df))

  # Subset only the available columns
  metadata <- df[, available_cols, drop = FALSE]

  # Convert list of padded spectrograms to a matrix
  spec_matrix <- do.call(rbind, padded_spectrograms)
  colnames(spec_matrix) <- paste0("V_", seq_len(ncol(spec_matrix)))

  # Combine metadata with spectrogram features
  res <- cbind(metadata, as.data.frame(spec_matrix))

  cat(sprintf("\nSuccessfully processed %d of %d segments.\n",
              length(valid_spectrograms), nrow(x)))

  # Add metadata as attributes
  attr(res, "max_length") <- max_length

  return(res)
}

#' @rdname extract_spec
#' @export
extract_spec.Sap <- function(x,
                             segment_type = c("segments", "syllables"),
                             sample_percent = NULL,
                             balanced = FALSE,
                             labels = NULL,
                             seed = 222,
                             cores = NULL,
                             wl = 256,
                             ovlp = 20,
                             wn = "hanning",
                             freq_range = c(1,10),
                             fftw = TRUE,
                             verbose = TRUE,
                             ...) {
  if(verbose) message(sprintf("\n=== Starting Spectrogram Extraction ===\n"))

  # Input validation
  segment_type <- match.arg(segment_type)
  segments_df <- x[[segment_type]]

  if (!inherits(segments_df, "segment") || nrow(segments_df) == 0) {
    stop("No segments found in the specified segment type")
  }

  # Select and balance segments
  segments_df <- select_segments(segments_df,
                                 labels = labels,
                                 balanced = balanced,
                                 sample_percent = sample_percent,
                                 seed = seed)

  # Process segments using default method
  if (fftw) ensure_pkgs("fftw")

  result <- extract_spec.default(segments_df,
                                 wav_dir = x$base_path,
                                 cores = cores,
                                 wl = wl,
                                 ovlp = ovlp,
                                 wn = wn,
                                 freq_range = freq_range,
                                 fftw = fftw,
                                 ...)

  if (nrow(result) < nrow(segments_df)) {
    filtered_count <- nrow(segments_df) - nrow(result)
    filtered_percent <- round(100 * filtered_count / nrow(segments_df), 1)

    warning_msg <- sprintf(
      "\n %d segments (%.1f%%) were filtered out during processing.\nRecommendations to improve processing success rate:\n1. Try reducing window length (wl)\n2. Try increasing freq_range or set it to NULL",
      filtered_count, filtered_percent
    )

    warning(warning_msg)
  }

  # Add attributes to result
  attr(result, "segment_indices") <- which(x[[segment_type]]$filename %in% segments_df$filename &
                                             x[[segment_type]]$start_time %in% segments_df$start_time &
                                             x[[segment_type]]$end_time %in% segments_df$end_time)
  attr(result, "segment_type") <- segment_type

  # Update SAP object's features
  feature_type <- sub("s$", "", segment_type)  # Remove 's' from end
  x$features[[feature_type]][["spectrogram"]] <- result

  # Add message about data access
  if(verbose) {
    cat(sprintf("\nAccess spectrogram data via: x$features$%s$spectrogram\n", feature_type))
  }

  # Return updated SAP object invisibly
  invisible(x)
}


#' Extract spectrogram for a single audio segment
#'
#' @keywords internal
extract_single_spec <- function(x, wav_dir, wl, ovlp, wn, fftw, freq_range, ...) {
  tryCatch({
    # Get file path
    file_path <- construct_wav_path(x, wav_dir = wav_dir)
    if (!file.exists(file_path)) {
      warning(sprintf("File not found: %s", file_path))
      return(NULL)
    }

    # Read wave file for specific time segment
    wave_data <- tuneR::readWave(
      file_path,
      from = x$start_time,
      to = x$end_time,
      units = "seconds"
    )

    # Apply frequency range if provided
    sr <- wave_data@samp.rate
    if (is.null(freq_range)) {
      freq_range <- c(0.5, sr/2000)
    }

    # Compute spectrogram
    stft_result <- seewave::spectro(
      wave_data,
      f = wave_data@samp.rate,
      wl = wl,
      ovlp = ovlp,
      wn = wn,
      fftw = fftw,
      flimd = freq_range,
      plot = FALSE,
      osc = TRUE,
      cont = TRUE
    )

    # Extract amplitude array and flatten it
    amp_array <- as.vector(stft_result$amp)

    return(amp_array)
  }, error = function(e) {
    warning(sprintf("Error processing row: %s", e$message))
    return(NULL)
  })
}

#' Pad an array with zeros to match the target length
#'
#' @param arr Array to pad
#' @param target_length Target length after padding
#' @return Padded array
#'
#' @keywords internal
pad_array <- function(arr, target_length) {
  orig_length <- length(arr)

  # If already the right length, return as is
  if (orig_length == target_length) {
    return(arr)
  }

  # Calculate padding needed on each side
  pad_before <- floor((target_length - orig_length) / 2)
  pad_after <- target_length - orig_length - pad_before

  # Create padded array
  padded <- numeric(target_length)

  # Insert original array in the center
  padded[(pad_before + 1):(pad_before + orig_length)] <- arr

  return(padded)
}


# Visualize Song ----------------------------------------------------------
# Update date : Feb. 7, 2025

#' Visualize Song Data
#'
#' @description
#' Creates spectrograms and visualizes acoustic data from WAV files or SAP objects.
#'
#' @param x An object to visualize, either a file path or SAP object
#' @param start_time_in_second Numeric start time in seconds
#' @param end_time_in_second Numeric end time in seconds
#' @param fft_window_size Size of FFT window (default: 512 for default, 1024 for SAP)
#' @param overlap Overlap between windows (default: 0.5 for default, 0.75 for SAP)
#' @param dark_mode For default method: Use dark theme (default: TRUE)
#' @param legend For default method: Show spectrogram legend (default: FALSE)
#' @param indices For SAP objects: Numeric vector of specific indices to visualize
#' @param template_clips Logical. For SAP objects: whether to visualize original songs (FALSE)
#'        or template clips (TRUE) (default: FALSE)
#' @param n_samples For SAP objects: Number of samples to visualize if indices is NULL.
#'          Default is 6 or max available
#' @param random For SAP objects: Randomly sample songs if TRUE
#' @param keep.par Preserve plotting parameters
#' @param verbose Print processing messages
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For WAV files:
#' \itemize{
#'   \item Creates single spectrogram using FFmpeg's FFT
#'   \item Customizable time range and FFT settings
#'   \item Optional dark mode and legend
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Creates multi-panel spectrograms
#'   \item Supports random or sequential sampling
#'   \item Maintains plotting state for sequential viewing
#'   \item Adds day and label information to plots
#' }
#'
#' @return
#' Generates spectrogram plot(s) and returns the input object invisibly.
#'
#' @examples
#' \dontrun{
#' # Visualize a single WAV file
#' visualize_song("path/to/song.wav",
#'                start_time_in_second = 10,
#'                end_time_in_second = 20)
#'
#' # Basic visualization from SAP object
#' visualize_song(sap_object, n_sample = 4)
#'
#' # Visualize specific indices with custom FFT settings
#' visualize_song(sap_object,
#'                indices = c(1, 3, 5),
#'                fft_window_size = 2048,
#'                overlap = 0.8)
#'
#' # Sequential visualization with time ranges
#' visualize_song(sap_object,
#'                n_sample = 6,
#'                random = FALSE,
#'                start_time_in_second = rep(0, 6),
#'                end_time_in_second = rep(5, 6))
#' }
#'
#' @rdname visualize_song
#' @export
visualize_song <- function(x, ...) {
  UseMethod("visualize_song")
}

#' @rdname visualize_song
#' @export
visualize_song.default <- function(x,  # wav file path
                                   start_time_in_second = NULL,
                                   end_time_in_second = NULL,
                                   fft_window_size = 1024,
                                   overlap = 0.5,
                                   dark_mode = TRUE,
                                   legend = FALSE,
                                   keep.par = FALSE,
                                   verbose = TRUE,
                                   ...) {
  # Normalize base path
  wav_path <- normalizePath(x, mustWork = TRUE)

  # Validate file path
  if (!file.exists(wav_path)) {
    stop("File does not exist: ", x)
  }

  # Determine start and end times
  start_time <- if (!is.null(start_time_in_second)) start_time_in_second else 0
  end_time <- if (!is.null(end_time_in_second)) end_time_in_second else NULL

  # Read FFT data and create plot
  fft_data <- av::read_audio_fft(
    wav_path,  # individual wav file path
    window = av::hanning(fft_window_size),
    overlap = overlap,
    start_time = start_time,
    end_time = end_time
  )

  # Create plot
  plot(fft_data, dark = dark_mode, legend = legend, keep.par = keep.par)
  if (verbose) {
    cat("Song visualization completed for:", basename(x), "\n")
  }
}

#' @rdname visualize_song
#' @export
visualize_song.Sap <- function(x,  # sap object
                              template_clips = FALSE,
                              indices = NULL,
                              n_samples = NULL,
                              random = TRUE,
                              start_time_in_second = NULL,
                              end_time_in_second = NULL,
                              fft_window_size = 1024,
                              overlap = 0.75,
                              keep.par = TRUE,
                              verbose = FALSE,
                              ...) {
  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  # Check if clips are requested but not available
  if ( template_clips) {
    if (is.null(x$templates$template_info)) {
      stop("No templates found in this SAP object. Run create_audio_clip() first.")
    }

    # Use template data
    total_samples <- nrow(x$templates$template_info)
    source_label <- "templates"
  } else {
    # Use regular metadata
    total_samples <- nrow(x$metadata)
    source_label <- "songs"
  }

  # Set default n_samples if null
  if (is.null(n_samples)) {
    # Default is 6, but limit to total available if less than 6
    n_samples <- min(6, total_samples)
  } else {
    # User specified n_samples, ensure it doesn't exceed available samples
    n_samples <- min(n_samples, total_samples)
  }

  # Handle indices or sampling
  if (is.null(indices)) {
    # Validate n_samples
    if (n_samples < 1) {
      stop("n_samples must be at least 1")
    }

    # Select indices based on random parameter
    if (random) {
      indices <- sample(1:total_samples, size = n_samples)
    } else {
      # Use the SAP object itself to store state instead of global environment
      if (is.null(x$._state)) {
        x$._state <- list()
      }

      # Create a unique key for tracking indices
      state_key <- paste0("idx_", make.names(source_label))

      # Get or initialize the last index
      if (is.null(x$._state[[state_key]])) {
        x$._state[[state_key]] <- 0
      }
      last_idx <- x$._state[[state_key]]

      # Calculate next set of indices
      start_idx <- last_idx + 1
      end_idx <- min(start_idx + n_samples - 1, total_samples)

      indices <- start_idx:end_idx

      # Update last plotted index
      x$._state[[state_key]] <- end_idx

      # Reset if we've reached the end
      if (end_idx >= total_samples) {
        x$._state[[state_key]] <- 0
        message("Reached the end of samples. Next call will start from the beginning.")
      }
    }
  } else {
    # Validate that indices are in range
    if ( template_clips) {
      if (any(indices > total_samples) || any(indices < 1)) {
        stop(sprintf("Indices must be between 1 and %d for templates", total_samples))
      }
    } else {
      if (any(indices > nrow(x$metadata)) || any(indices < 1)) {
        stop(sprintf("Indices must be between 1 and %d for songs", nrow(x$metadata)))
      }
    }
  }

  # Validate time vectors if provided
  if (!is.null(start_time_in_second) && length(start_time_in_second) != length(indices)) {
    stop("The length of start_time_in_second must match the length of indices.")
  }
  if (!is.null(end_time_in_second) && length(end_time_in_second) != length(indices)) {
    stop("The length of end_time_in_second must match the length of indices.")
  }

  # Print index information before plotting
  cat(sprintf("\nPlotting the following %s:\n", source_label))
  for (i in seq_along(indices)) {
    index <- indices[i]

    if ( template_clips) {
      # Display clip/template information
      clip_name <- x$templates$template_info$clip_name[index]
      cat(sprintf("Index %d: Clip: %s\n", index, clip_name))
    } else {
      # Display song information
      filename <- x$metadata$filename[index]
      day_post_hatch <- x$metadata$day_post_hatch[index]
      label <- if ("label" %in% names(x$metadata)) {
        sprintf(" (Label: %s)", x$metadata$label[index])
      } else {
        ""
      }

      cat(sprintf("Index %d: Day %s, File: %s%s\n",
                  index,
                  day_post_hatch,
                  filename,
                  label))
    }
  }
  cat("\n")

  # Set up the plotting panel
  n_plots <- length(indices)
  n_rows <- ceiling(sqrt(n_plots))
  n_cols <- ceiling(n_plots/n_rows)

  # Store original par settings if keep.par is TRUE
  if (keep.par) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
  }

  # Set up the panel layout
  par(mfrow = c(n_rows, n_cols))

  # Process each index
  for (i in seq_along(indices)) {
    index <- indices[i]

    if ( template_clips) {
      # Get clip path directly from template_info
      clip_name <- x$templates$template_info$clip_name[index]
      song_path <- x$templates$template_info$clip_path[index]

      # Set title text for clip
      title_text <- sprintf("Template clips: %s", clip_name)
    } else {
      # Construct path to original song file
      song_path <- file.path(x$base_path,
                             x$metadata$day_post_hatch[index],
                             x$metadata$filename[index])

      # Set title text for original file
      title_text <- sprintf("Day %s%s",
                            x$metadata$day_post_hatch[index],
                            if ("label" %in% names(x$metadata))
                              sprintf(" (%s)", x$metadata$label[index]) else "")
    }

    # Determine start and end times
    start_time <- if (!is.null(start_time_in_second)) start_time_in_second[i] else 0
    end_time <- if (!is.null(end_time_in_second)) end_time_in_second[i] else NULL

    # Call default method for each file
    visualize_song.default(
      x = song_path,
      fft_window_size = fft_window_size,
      overlap = overlap,
      start_time_in_second = start_time,
      end_time_in_second = end_time,
      dark_mode = TRUE,
      legend = FALSE,
      keep.par = TRUE,
      verbose = verbose,
      ...
    )

    # Add title to each plot
    title(title_text, line = 0.5)
  }

  # Return the modified object with state information
  invisible(x)
}


# Visualize Song Segments -------------------------------------------------
# Update date : Feb. 7, 2025

#' Visualize Song Segments
#'
#' @description
#' Creates multi-panel spectrogram visualizations of audio segments from various sources.
#'
#' @param x An object to visualize, either a data frame or SAP object
#' @param wav_file_dir For default method: Directory containing WAV files
#' @param n_samples Number of samples to display
#' @param seed Random seed for sample selection
#' @param fft_window_size Size of FFT window (default: 1024)
#' @param overlap Overlap between windows (default: 0.75)
#' @param dark_mode Use dark theme (default: TRUE)
#' @param legend Show spectrogram legend (default: FALSE)
#' @param segment_type For SAP objects: Type of segments ('motifs', 'bouts', 'segments')
#' @param labels For SAP objects: Labels to include
#' @param clusters For SAP objects: Specific clusters to visualize
#' @param by_column For SAP objects: Arrange by columns (default: TRUE)
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For data frames:
#' \itemize{
#'   \item Requires columns: filename, start_time, end_time
#'   \item Optional column: day_post_hatch for hierarchical file structure
#'   \item Supports random sampling of segments
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Supports visualization by labels and/or clusters
#'   \item Flexible sampling within groups
#'   \item Customizable layout (by row or column)
#'   \item Automatic handling of file paths
#' }
#'
#' @return
#' Generates a multi-panel spectrogram plot and returns the input object invisibly.
#'
#' @examples
#' \dontrun{
#' # Visualize from data frame
#' song_df <- data.frame(
#'   filename = c("song1.wav", "song2.wav"),
#'   start_time = c(0, 10),
#'   end_time = c(30, 40)
#' )
#' visualize_segments(song_df,
#'                    wav_file_dir = "path/to/wav/files",
#'                    n_samples = 5)
#'
#' # Basic SAP object visualization
#' visualize_segments(sap_object,
#'                    segment_type = "motifs",
#'                    n_samples = 3)
#'
#' # Cluster-specific visualization
#' visualize_segments(sap_object,
#'                    segment_type = "segments",
#'                    clusters = c(1, 2),
#'                    labels = c("a", "b"),
#'                    n_samples = 4)
#'
#' # Custom layout
#' visualize_segments(sap_object,
#'                    segment_type = "motifs",
#'                    by_column = FALSE,
#'                    fft_window_size = 2048)
#' }
#'
#' @seealso \code{\link{visualize_song}} for single file visualization
#'
#' @rdname visualize_segments
#' @export
visualize_segments <- function(x, ...) {
  UseMethod("visualize_segments")
}

#' @rdname visualize_segments
#' @export
visualize_segments.default  <- function(x,
                                        wav_file_dir,
                                        n_samples = NULL,
                                        seed = NULL,
                                        fft_window_size = 1024,
                                        overlap = 0.75,
                                        dark_mode = TRUE,
                                        legend = FALSE,
                                        ...) {
  # Input validation
  required_cols <- c("filename", "start_time", "end_time")
  if (!all(required_cols %in% colnames(x))) {
    stop("Data frame must contain columns: ",
         paste(required_cols, collapse = ", "))
  }

  # Validate directory
  if (!dir.exists(wav_file_dir)) {
    stop("WAV file directory does not exist: ", wav_file_dir)
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Determine number of samples to plot
  total_rows <- nrow(x)
  n_to_plot <- if (is.null(n_samples)) total_rows else min(n_samples, total_rows)

  # Sample rows if necessary
  if (n_to_plot < total_rows) {
    selected_rows <- sample(1:total_rows, n_to_plot)
    x <- x[selected_rows, ]
  }

  # Set up multi-panel plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(n_to_plot, 1))

  # Process each selected file
  for (i in 1:n_to_plot) {
    # Use the path constructor function
    if (is.null(x$day_post_hatch)) {
      full_path <- file.path(wav_file_dir, x$filename[i])
    } else {
      full_path <- file.path(wav_file_dir, x$day_post_hatch[i], x$filename[i])
    }

    visualize_song.default(
      x = full_path,
      start_time_in_second = x$start_time[i],
      end_time_in_second = x$end_time[i],
      fft_window_size = fft_window_size,
      overlap = overlap,
      dark_mode = dark_mode,
      legend = legend,
      keep.par = TRUE,
      ...
    )
  }
}

#' @rdname visualize_segments
#' @export
visualize_segments.Sap <- function(x,
                                   segment_type = c("motifs", "bouts", "segments"),
                                   labels = NULL,
                                   clusters = NULL,
                                   n_samples = NULL,
                                   by_column = TRUE,
                                   seed = NULL,
                                   fft_window_size = 1024,
                                   overlap = 0.75,
                                   dark_mode = TRUE,
                                   legend = FALSE,
                                   ...) {
  # Input validation
  segment_type <- match.arg(segment_type)

  # Get appropriate data frame
  if(segment_type == "segments") {
    if(is.null(x$features$segment$feat.embeds)) {
      stop("No feature embeddings found for segments")
    }
    segments_df <- x$features$segment$feat.embeds
    cluster_col <- "cluster"

    # Convert cluster factor to numeric
    segments_df[[cluster_col]] <- as.numeric(as.character(segments_df[[cluster_col]]))
  } else {
    segments_df <- x[[segment_type]]
    cluster_col <- NULL
  }

  # Validate segments dataframe
  if (!inherits(segments_df, "data.frame") || nrow(segments_df) == 0) {
    stop(sprintf("No %s found in SAP object",
                 ifelse(segment_type == "segments", "segment features", segment_type)))
  }

  # Handle cluster filtering
  if(!is.null(clusters) && !is.null(cluster_col)) {
    if(!cluster_col %in% colnames(segments_df)) {
      stop("Cluster column not found in segment data")
    }
    segments_df <- segments_df[segments_df[[cluster_col]] %in% clusters, ]
    if(nrow(segments_df) == 0) {
      stop("No segments found in specified clusters")
    }
  }

  # Get unique labels and clusters
  available_labels <- unique(segments_df$label)
  if(!is.null(cluster_col)) {
    available_clusters <- unique(segments_df[[cluster_col]])
  } else {
    available_clusters <- NULL
  }

  # Validate requested labels
  if(is.null(labels)) {
    labels <- available_labels
  } else {
    invalid_labels <- setdiff(labels, available_labels)
    if(length(invalid_labels) > 0) {
      stop("Invalid labels requested: ", paste(invalid_labels, collapse = ", "))
    }
  }

  # Validate requested clusters
  if(!is.null(clusters) && !is.null(available_clusters)) {
    invalid_clusters <- setdiff(clusters, available_clusters)
    if(length(invalid_clusters) > 0) {
      stop("Invalid clusters requested: ", paste(invalid_clusters, collapse = ", "))
    }
    clusters <- intersect(clusters, available_clusters)
  } else {
    clusters <- available_clusters
  }

  # Set seed if provided
  if(!is.null(seed)) set.seed(seed)

  # Determine grouping structure
  if(!is.null(cluster_col) && !is.null(clusters)) {
    group_vars <- list(labels = labels, clusters = clusters)
    n_groups <- length(labels) * length(clusters)
  } else {
    group_vars <- list(labels = labels)
    n_groups <- length(labels)
  }

  # Determine number of samples per group
  if(is.null(n_samples)) {
    n_samples <- if(n_groups > 0) {
      min(sapply(group_vars, function(g) {
        if(length(g) > 0) min(table(segments_df[segments_df$label %in% labels, ]$label))
        else nrow(segments_df)
      }))
    } else {
      5  # default sample size
    }
  }

  # Sample rows for each group
  selected_rows <- list()
  if(!is.null(cluster_col) && !is.null(clusters)) {
    for(cl in clusters) {
      for(lbl in labels) {
        group_rows <- which(segments_df$label == lbl & segments_df[[cluster_col]] == cl)
        if(length(group_rows) > 0) {
          selected_rows[[paste(lbl, cl, sep="_")]] <- sample_rows(group_rows, n_samples)
        }
      }
    }
  } else {
    for(lbl in labels) {
      group_rows <- which(segments_df$label == lbl)
      if(length(group_rows) > 0) {
        selected_rows[[lbl]] <- sample_rows(group_rows, n_samples)
      }
    }
  }

  # Combine selected rows
  selected_indices <- unlist(selected_rows)
  sampled_segments <- segments_df[selected_indices, ]

  # Set up multi-panel plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # Calculate layout dimensions
  if (!is.null(clusters)) {
    # Create group labels
    groups <- expand.grid(label = labels, cluster = clusters)
    n_groups <- nrow(groups)

    if (by_column) {
      # Columns: Groups, Rows: Samples
      layout_matrix <- matrix(1:(n_groups * n_samples),
                              nrow = n_samples,
                              ncol = n_groups,
                              byrow = FALSE)
    } else {
      # Rows: Groups, Columns: Samples
      layout_matrix <- matrix(1:(n_groups * n_samples),
                              nrow = n_groups,
                              ncol = n_samples,
                              byrow = TRUE)
    }
  } else {
    # Original label-only layout
    if (by_column) {
      layout_matrix <- matrix(1:(length(labels) * n_samples),
                              nrow = n_samples,
                              ncol = length(labels))
    } else {
      layout_matrix <- matrix(1:(length(labels) * n_samples),
                              nrow = length(labels),
                              ncol = n_samples,
                              byrow = TRUE)
    }
  }
  layout(layout_matrix)

  # Main plotting loop
  if(!is.null(clusters)) {
    # Plot with cluster grouping
    for(cl in clusters) {
      for(lbl in labels) {
        cluster_segments <- sampled_segments[sampled_segments$label == lbl &
                                               sampled_segments[[cluster_col]] == cl, ]
        if(nrow(cluster_segments) > 0) {
          plot_group(cluster_segments, x$base_path, n_samples,
                     fft_window_size, overlap, dark_mode, legend,
                     by_column, lbl, cl)
        }
      }
    }
  } else {
    # Plot without cluster grouping
    for(lbl in labels) {
      label_segments <- sampled_segments[sampled_segments$label == lbl, ]
      if(nrow(label_segments) > 0) {
        plot_group(label_segments, x$base_path, n_samples,
                   fft_window_size, overlap, dark_mode, legend,
                   by_column, lbl)
      }
    }
  }

  invisible(x)
}

#' Sample Rows from a Vector
#'
#' @description
#' An internal helper function to sample a specified number of rows
#' from a given vector of row indices.
#'
#' @param rows A vector of row indices
#' @param n Number of rows to sample
#'
#' @return
#' A vector of sampled row indices
#'
#' @keywords internal
sample_rows <- function(rows, n) {
  if(length(rows) > n) sample(rows, n) else rows[1:min(length(rows), n)]
}

#' Plot Group of Segments
#'
#' @description
#' An internal helper function to plot spectrograms for a group of segments.
#'
#' @param segments A data frame of segment information
#' @param base_path Base directory path for audio files
#' @param n_samples Number of samples to plot
#' @param fft_window_size Size of FFT window
#' @param overlap Overlap between windows
#' @param dark_mode Use dark theme
#' @param legend Show spectrogram legend
#' @param by_column Arrange plots by columns
#' @param label Label for the group of segments
#' @param cluster Optional cluster identifier
#'
#' @return
#' Generates a plot of segment spectrograms as a side effect
#'
#' @keywords internal
plot_group <- function(segments, base_path, n_samples,
                       fft_window_size, overlap, dark_mode, legend,
                       by_column, label, cluster = NULL) {
  for(i in seq_len(nrow(segments))) {
    full_path <- file.path(base_path,
                           segments$day_post_hatch[i],
                           segments$filename[i])

    visualize_song.default(
      x = full_path,
      start_time_in_second = segments$start_time[i],
      end_time_in_second = segments$end_time[i],
      fft_window_size = fft_window_size,
      overlap = overlap,
      dark_mode = dark_mode,
      legend = legend,
      keep.par = TRUE
    )

    # Add labels
    if(i == 1) {
      title_text <- if(!is.null(cluster)) {
        paste(label, "\nCluster", cluster)
      } else {
        label
      }

      if(by_column) {
        title(title_text, cex.main = 0.8)
      } else {
        mtext(title_text, side = 2, line = 4, cex = 0.7)
      }
    }
  }
}
