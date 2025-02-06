
# Visualize Song ----------------------------------------------------------
# Update date : Feb. 4, 2025

#' Visualize Song Data
#'
#' @description
#' A generic function to create spectrograms and visualize acoustic data
#' from WAV files or SAP objects.
#'
#' @param x An object to visualize, either a character path to a WAV file
#'          or a SAP object
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' This generic function supports visualization of song data through two methods:
#' \itemize{
#'   \item Default method for individual WAV files
#'   \item SAP object method for multiple file visualization
#' }
#'
#' @return
#' Generates a spectrogram plot. Returns the input object invisibly.
#'
#' @examples
#' \dontrun{
#' # Visualize a single WAV file
#' visualize_song("path/to/song.wav")
#'
#' # Visualize songs from a SAP object
#' visualize_song(sap_object)
#'
#' # Visualize specific indices from a SAP object
#' visualize_song(sap_object, indices = c(1, 3, 5))
#' }
#'
#' @export
visualize_song <- function(x, ...) {
  UseMethod("visualize_song")
}

#' Visualize Individual Song Spectrogram
#'
#' @description
#' Creates a spectrogram from a single WAV file using FFmpeg's FFT
#' through the av package.
#'
#' @param x Path to a WAV file
#' @param start_time_in_second Numeric start time in seconds (optional)
#' @param end_time_in_second Numeric end time in seconds (optional)
#' @param fft_window_size Size of FFT window (default: 512)
#' @param overlap Overlap between windows (default: 0.5)
#' @param dark_mode Use dark theme for visualization (default: TRUE)
#' @param legend Show spectrogram legend (default: FALSE)
#' @param keep.par Preserve plotting parameters (default: FALSE)
#' @param verbose Print processing messages (default: TRUE)
#' @param ... Additional arguments
#'
#' @details
#' Generates a spectrogram visualization for a single audio file with
#' customizable parameters for time range, FFT settings, and visual style.
#'
#' @return
#' Generates a spectrogram plot of the audio file
#'
#' @importFrom av read_audio_fft
#' @importFrom signal hanning
#' @export
visualize_song.default <- function(x,  # wav file path
                                   start_time_in_second = NULL,
                                   end_time_in_second = NULL,
                                   fft_window_size = 512,
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
    window = signal::hanning(fft_window_size),
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

#' Visualize Songs from a SAP Object
#'
#' @description
#' Creates spectrograms for multiple songs from a SAP object,
#' with options for sampling and customization.
#'
#' @param x A SAP object containing song recordings
#' @param indices Numeric vector of specific indices to visualize
#' @param n_sample Number of random samples to plot (default: 6)
#' @param random Randomly sample songs if TRUE, otherwise sample sequentially
#' @param start_time_in_second Optional start times for each song
#' @param end_time_in_second Optional end times for each song
#' @param fft_window_size Size of FFT window (default: 1024)
#' @param overlap Overlap between windows (default: 0.75)
#' @param keep.par Preserve plotting parameters (default: TRUE)
#' @param verbose Print processing messages (default: FALSE)
#' @param ... Additional arguments
#'
#' @details
#' Generates spectrograms for multiple songs from a SAP object with
#' flexible sampling and visualization options:
#' \itemize{
#'   \item Select specific indices or sample randomly/sequentially
#'   \item Customize time ranges for each song
#'   \item Control spectrogram visualization parameters
#' }
#'
#' @return
#' Generates a multi-panel spectrogram plot. Returns the input SAP object invisibly.
#'
#' @export
visualize_song.Sap <- function(x,  # sap object
                               indices = NULL,
                               n_sample = 6,
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

  # Handle indices or sampling
  if (is.null(indices)) {
    # Validate n_sample
    total_samples <- nrow(x$metadata)
    if (n_sample < 1 || n_sample > total_samples) {
      stop(sprintf("n_sample must be between 1 and %d", total_samples))
    }

    # Select indices based on random parameter
    if (random) {
      indices <- sample(1:total_samples, size = n_sample)
    } else {
      # Create or get the environment to store the last index
      if (!exists("._sap_state", envir = .GlobalEnv)) {
        assign("._sap_state", new.env(), envir = .GlobalEnv)
      }
      sap_state <- get("._sap_state", envir = .GlobalEnv)

      # Create a unique key for this SAP object based on its base path
      state_key <- paste0("idx_", digest::digest(x$base_path))

      # Get or initialize the last index
      if (!exists(state_key, envir = sap_state)) {
        assign(state_key, 0, envir = sap_state)
      }
      last_idx <- get(state_key, envir = sap_state)

      # Calculate next set of indices
      start_idx <- last_idx + 1
      end_idx <- min(start_idx + n_sample - 1, total_samples)

      indices <- start_idx:end_idx

      # Update last plotted index
      assign(state_key, end_idx, envir = sap_state)

      # Reset if we've reached the end
      if (end_idx >= total_samples) {
        assign(state_key, 0, envir = sap_state)
        message("Reached the end of samples. Next call will start from the beginning.")
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
  cat("\nPlotting the following songs:\n")
  for (i in seq_along(indices)) {
    index <- indices[i]
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
    song_path <- file.path(x$base_path,
                           x$metadata$day_post_hatch[index],
                           x$metadata$filename[index])

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
    title_text <- sprintf("Day %s%s",
                          x$metadata$day_post_hatch[index],
                          if ("label" %in% names(x$metadata))
                            sprintf(" (%s)", x$metadata$label[index]) else "")
    title(title_text, line = 0.5)
  }

  invisible(x)
}


# Visualize Song Segments -------------------------------------------------
# Update date : Feb. 4, 2025

#' Visualize Song Segments
#'
#' @description
#' A generic function to visualize audio segments from various sources,
#' supporting both data frames and SAP objects.
#'
#' @param x An object to visualize, either a data frame or a SAP object
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' This generic function supports visualization of song segments through two methods:
#' \itemize{
#'   \item Default method for visualizing segments from a data frame
#'   \item SAP object method for visualizing segments from a Sound Analysis Pro object
#' }
#'
#' @return
#' Generates a multi-panel spectrogram plot. Returns the input object invisibly.
#'
#' @examples
#' \dontrun{
#' # Visualize segments from a data frame
#' song_df <- data.frame(
#'   filename = c("song1.wav", "song2.wav"),
#'   start_time = c(0, 10),
#'   end_time = c(30, 40)
#' )
#' visualize_segments(song_df, wav_file_dir = "path/to/wav/files")
#'
#' # Visualize segments from a SAP object
#' visualize_segments(sap_object, segment_type = "motifs")
#' }
#'
#' @export
visualize_segments <- function(x, ...) {
  UseMethod("visualize_segments")
}

#' Visualize Segments from a Data Frame
#'
#' @description
#' Creates a multi-panel plot of audio segments specified in a data frame.
#'
#' @param x A data frame containing segment information
#' @param wav_file_dir Directory path containing the WAV files
#' @param n_samples Number of samples to display (default: all samples)
#' @param seed Random seed for sample selection
#' @param fft_window_size Size of the FFT window (default: 1024)
#' @param overlap Overlap between windows (default: 0.75)
#' @param dark_mode Use dark theme for visualization (default: TRUE)
#' @param legend Show spectrogram legend (default: FALSE)
#' @param ... Additional arguments
#'
#' @details
#' Generates a multi-panel spectrogram visualization for segments:
#' \itemize{
#'   \item Validates input data frame
#'   \item Samples segments if specified
#'   \item Creates individual spectrograms for each segment
#' }
#'
#' @return
#' Generates a multi-panel spectrogram plot
#'
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

#' Visualize Segments from a SAP Object
#'
#' @description
#' Creates a multi-panel plot of audio segments from a Sound Analysis Pro object.
#'
#' @param x A SAP object containing song recordings
#' @param segment_type Type of segments to visualize
#'        (one of "motifs", "bouts", "segments")
#' @param labels Labels to include in visualization
#' @param clusters Specific clusters to visualize
#' @param n_samples Number of samples per label/cluster
#' @param by_column Arrange labels/clusters by columns (default: TRUE)
#' @param seed Random seed for sample selection
#' @param fft_window_size Size of FFT window (default: 1024)
#' @param overlap Overlap between windows (default: 0.75)
#' @param dark_mode Use dark theme for visualization (default: TRUE)
#' @param legend Show spectrogram legend (default: FALSE)
#' @param ... Additional arguments
#'
#' @details
#' Generates a multi-panel spectrogram visualization for segments from a SAP object:
#' \itemize{
#'   \item Supports visualization by labels and/or clusters
#'   \item Flexible sampling of segments
#'   \item Customizable plot layout
#' }
#'
#' @return
#' Generates a multi-panel spectrogram plot. Returns the input SAP object invisibly.
#'
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
