# Amplitude Envelope module


# Extract Amplitude Envelope ----------------------------------------------
# Update date : Feb. 7, 2025

#' Calculate Amplitude Envelope for Audio Segment
#'
#' @description
#' Calculates the amplitude envelope for a specified segment of an audio file.
#'
#' @param segment_row A row from either find_motif output or segment class in SAP object
#' @param wav_dir Optional path to WAV files directory (default: NULL)
#' @param msmooth Numeric vector of length 2 for envelope smoothing:
#'        \itemize{
#'          \item First value: Window length in number of points
#'          \item Second value: Overlap between windows (percentage)
#'        }
#'        If NULL, no smoothing is applied
#' @param norm Logical, Whether to normalize envelope to range between 0 and 1 (default: FALSE)
#' @param plot Logical, whether to plot the envelope (default: FALSE)
#'
#' @details
#' The function processes audio segments with the following steps:
#' \itemize{
#'   \item Validates input segment data
#'   \item Constructs correct file path
#'   \item Reads specified portion of audio file
#'   \item Calculates amplitude envelope
#'   \item Optionally smooths and normalizes the envelope
#' }
#'
#' @return
#' A numeric vector containing the amplitude envelope values
#'
#' @examples
#' \dontrun{
#' # Basic envelope calculation
#' env <- amp_env(segments[1,], wav_dir = "path/to/wavs")
#'
#' # With smoothing and normalization
#' env <- amp_env(segments[1,],
#'                wav_dir = "path/to/wavs",
#'                msmooth = c(256, 50),
#'                norm = TRUE)
#'
#' # With plot
#' env <- amp_env(segments[1,],
#'                wav_dir = "path/to/wavs",
#'                msmooth = c(256, 50),
#'                plot = TRUE)
#' }
#'
#' @seealso \code{\link{find_motif}} for creating segment data
#'
#' @export
amp_env <- function(segment_row,
                    wav_dir = NULL,
                    msmooth = NULL,
                    norm = FALSE,
                    plot = FALSE) {

  # Check if input is valid
  if (!is.data.frame(segment_row) || nrow(segment_row) != 1) {
    stop("segment_row must be a single row from a data frame")
  }

  # Check for required timing columns
  if (!all(c("start_time", "end_time") %in% names(segment_row))) {
    stop("Input row must contain 'start_time' and 'end_time' columns")
  }

  # Validate msmooth if provided
  if (!is.null(msmooth)) {
    if (!is.numeric(msmooth) || length(msmooth) != 2) {
      stop("msmooth must be a numeric vector of length 2")
    }
  }

  # Construct file path
  sound_path <- construct_wav_path(segment_row, wav_dir = wav_dir)

  # Verify file exists
  if (!file.exists(sound_path)) {
    stop(sprintf("Sound file not found: %s", sound_path))
  }

  # Read wave file using timestamps
  wv <- tuneR::readWave(sound_path,
                        from = segment_row$start_time,
                        to = segment_row$end_time,
                        units = "seconds")

  # Calculate envelope
  env <- seewave::env(wv, msmooth = msmooth, norm = norm, plot = plot)

  return(env)
}


# Plot Heatmap ----------------------------------------------
# Update date : Feb. 7, 2025

#' Plot Heatmap of Amplitude Envelopes
#'
#' @description
#' Creates heatmap visualizations of amplitude envelopes from audio segments, supporting
#' multiple data sources and visualization options.
#'
#' @param x An object to visualize (data frame, SAP object, or matrix)
#' @param wav_dir For default method: Path to WAV files directory
#' @param msmooth Smoothing parameters c(window_length, overlap_percentage)
#' @param color_palette Function generating color palette
#' @param n_colors Number of colors in heatmap (default: 500)
#' @param contrast Contrast factor for visualization (default: 3)
#' @param segment_type For SAP objects: Type of segments ('motifs', 'bouts', 'syllables', 'segments')
#' @param sample_percent For SAP objects: Percentage to sample
#' @param balanced For SAP objects: Balance across labels
#' @param clusters Numeric vector of cluster IDs to filter
#' @param labels Optional vector of labels to include
#' @param cores For SAP objects: Number of processing cores
#' @param seed For SAP objects: Random seed (default: 222)
#' @param ordered For SAP objects: Order by embeddings
#' @param descending For SAP objects: Direction of ordering
#' @param padding_quantile For SAP objects: Quantile for bout padding (default: 0.9)
#' @param window For SAP objects: Numeric vector of length 2 (pre, post) specifying time windows
#'   in seconds around the alignment point for bouts. When NULL (default),
#'   windows are auto-calculated. Ignored for motifs. Example: \code{c(0.5, 1.5)}.
#' @param reference_lines Numeric vector specifying time positions (in seconds)
#'   relative to the alignment point where vertical dashed reference lines should
#'   be drawn. Can include both positive and negative values.
#' @param reference_line_color Character vector specifying the color(s) of the
#'   reference lines. Can be a single color (applied to all lines) or a vector
#'   of colors (one for each line). If fewer colors than lines are provided,
#'   colors will be recycled.
#' @param ylabel Character string specifying the label for the y-axis. (default: "Labels")
#' @param verbose For SAP objects: Print progress messages
#' @param main For matrix method: Plot title
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For data frames:
#' \itemize{
#'   \item Requires columns: filename, start_time, end_time
#'   \item Calculates envelopes for each segment
#'   \item Creates matrix of aligned envelopes
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Supports multiple segment types
#'   \item Optional balanced sampling
#'   \item Parallel processing support
#'   \item Ordering by feature embeddings
#' }
#'
#' For matrices:
#' \itemize{
#'   \item Direct visualization of pre-computed envelopes
#'   \item Label-based organization
#'   \item Visual separation between groups
#' }
#'
#' @return
#' For default method: List containing segments, matrix, and plot
#' For SAP objects: Updated object with amplitude features
#' For matrices: Lattice plot object
#'
#' @examples
#' \dontrun{
#' # Basic usage with data frame
#' plot_heatmap(segments, wav_dir = "path/to/wavs")
#'
#' # SAP object with options
#' plot_heatmap(sap_obj,
#'              segment_type = "motifs",
#'              balanced = TRUE,
#'              ordered = TRUE)
#'
#' # Matrix with specific labels
#' plot_heatmap(amp_matrix,
#'              labels = c("a", "b"),
#'              contrast = 2)
#'
#' # Advanced SAP object usage
#' plot_heatmap(sap_obj,
#'              segment_type = "bouts",
#'              sample_percent = 80,
#'              cores = 4,
#'              ordered = TRUE,
#'              descending = FALSE)
#' }
#'
#' @rdname plot_heatmap
#' @export
plot_heatmap <- function(x, ...) {
  UseMethod("plot_heatmap")
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.default <- function(x,
                                 wav_dir = NULL,
                                 msmooth = c(256,50),
                                 color_palette = NULL,
                                 n_colors = 500,
                                 contrast = 3,
                                 ...) {

  # Check required columns
  required_cols <- c("filename", "start_time", "end_time")
  missing_cols <- required_cols[!required_cols %in% names(x)]
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Check if data is empty
  if (nrow(x) == 0) {
    stop("Input data frame is empty")
  }

  # Handle wav_dir
  if (is.null(wav_dir) && is.null(attr(x, "wav_dir"))) {
    stop("wav_dir must be provided either as argument or attribute")
  }

  # # If wav_dir is provided as argument and not as attribute, add it as attribute
  # if (!is.null(wav_dir) && is.null(attr(x, "wav_dir"))) {
  #   attr(x, "wav_dir") <- normalizePath(wav_dir, mustWork = FALSE)
  # }

  # Generate amplitude envelope matrix using appropriate progress bar function
  if (Sys.info()["sysname"] == "Darwin") {
    amp_list <- pbmcapply::pbmclapply(1:nrow(x),
                                      function(i) amp_env(x[i,],
                                                          wav_dir = wav_dir,
                                                          msmooth = msmooth))
    amp_matrix <- do.call(cbind, amp_list)
  } else {
    amp_matrix <- pbapply::pbsapply(1:nrow(x),
                                    function(i) amp_env(x[i,],
                                                        wav_dir = wav_dir,
                                                        msmooth = msmooth),
                                    simplify = TRUE)
  }

  # Calculate time window from first segment
  time_window <- round(x$end_time[1] - x$start_time[1], 2)

  # Set color palette
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("black", "red", "yellow", "white"))
  }

  # Create heatmap
  heatmap <- lattice::levelplot(amp_matrix,
                                col.regions = color_palette(n_colors),
                                at = seq(min(amp_matrix),
                                         max(amp_matrix)/contrast,
                                         length.out = n_colors),

                                border = "transparent",
                                ylab = "Rendition",
                                xlab = "Time (s)",
                                main = "Amplitude Envelope Heatmap",
                                scales = list(
                                  x = list(
                                    at = seq(0, nrow(amp_matrix), length.out = 5),
                                    labels = sprintf("%.1f", seq(0, time_window, length.out = 5))
                                  )
                                ),
                                pretty = TRUE,
                                aspect = "fill",
                                colorkey = list(
                                  space = "right",
                                  width = 1,
                                  height = 0.75,
                                  labels = list(
                                    at = seq(min(amp_matrix),
                                             max(amp_matrix) / contrast,
                                             length.out = 3),
                                    cex = 1,
                                    col = "black",
                                    labels = sprintf("%.1f", seq(0, 1, length.out = 3))
                                  )
                                ))

  print(heatmap)

  # Return results
  result <- list(
    segments = x,
    amp_matrix = amp_matrix,
    plot = heatmap
  )

  class(result) <- c("segment_heatmap", class(result))
  invisible(result)
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.matrix <- function(x,
                                labels = NULL,
                                color_palette = NULL,
                                n_colors = 500,
                                contrast = 3,
                                reference_lines = NULL,
                                reference_line_color = "black",
                                main = "Amplitude Envelope Heatmap",
                                ylabel = "Labels",
                                ...) {
  # Check if input is matrix
  if (!is.matrix(x)) {
    stop("Input must be a matrix")
  }

  # Check for required attributes
  if (is.null(attr(x, "time_window"))) {
    stop("Matrix must have 'time_window' attribute")
  }

  if (is.null(colnames(x))) {
    stop("Matrix must have column names as labels")
  }

  # Store time_window before subsetting
  time_window <- attr(x, "time_window")
  if (!is.numeric(time_window) || length(time_window) != 1) {
    stop("time_window attribute must be a single numeric value")
  }

  # Handle label subsetting if provided
  if (!is.null(labels)) {
    if (!all(labels %in% colnames(x))) {
      stop("Some provided labels not found in matrix")
    }
    x <- x[, colnames(x) %in% labels]
    attr(x, "time_window") <- time_window  # Preserve just the time_window attribute
  }

  # Create reversed matrix for plotting
  reversed_amp_matrix <- x[, ncol(x):1]

  # Calculate positions using reserved order
  current_labels <- colnames(x)

  # Extract numeric part from dph labels
  if (any(grepl("dph", current_labels))) {
    numeric_labels <- gsub("dph", "", current_labels)
    numeric_values <- as.numeric(numeric_labels)
    current_labels <- numeric_labels
  }

  # reverse the order
  reversed_labels <- rev(unique(current_labels))

  ordered_labels <- factor(current_labels, levels = unique(current_labels))
  samples_per_label <- rev(table(ordered_labels))  # rev() for bottom-to-top plotting

  cumulative_positions <- cumsum(c(0, head(samples_per_label, -1)))
  label_positions <- cumulative_positions + samples_per_label/2
  hline_positions <- cumsum(samples_per_label)[-length(samples_per_label)]

  # Handle reference lines
  if (!is.null(reference_lines)) {
    # Ensure reference_line_color has the same length as reference_lines
    if (length(reference_line_color) == 1) {
      reference_line_color <- rep(reference_line_color, length(reference_lines))
    } else if (length(reference_line_color) != length(reference_lines)) {
      reference_line_color <- rep(reference_line_color, length.out = length(reference_lines))
    }
  }

  # Set color palette
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("black", "red", "yellow", "white"))
  }

  # Create heatmap
  heatmap <- lattice::levelplot(reversed_amp_matrix,
                                col.regions = color_palette(n_colors),
                                at = seq(min(reversed_amp_matrix),
                                         max(reversed_amp_matrix)/contrast,
                                         length.out = n_colors),
                                border = "transparent",
                                ylab = ylabel,
                                xlab = "Time (s)",
                                main = main,
                                scales = list(
                                  x = list(
                                    at = seq(0, nrow(reversed_amp_matrix), length.out = 5),
                                    labels = sprintf("%.1f", seq(0, time_window, length.out = 5))
                                  ),
                                  y = list(
                                    at = label_positions,
                                    labels = reversed_labels
                                  )
                                ),
                                panel = function(x, y, z, ...) {
                                  lattice::panel.levelplot(x, y, z, ...)

                                  # Draw horizontal lines between labels
                                  lattice::panel.abline(h = hline_positions,
                                                        col = "white",
                                                        lwd = 2)

                                  # Draw vertical reference lines if provided
                                  if (!is.null(reference_lines)) {
                                    for (i in seq_along(reference_lines)) {
                                      lattice::panel.abline(v = reference_lines[i],
                                                            col = reference_line_color[i],
                                                            lwd = 3,
                                                            lty = 2)  # dashed lines
                                    }
                                  }
                                },
                                aspect = "fill",
                                colorkey = list(
                                  space = "right",
                                  width = 1,
                                  height = 0.75,
                                  labels = list(
                                    at = seq(min(x),
                                             max(x) / contrast,
                                             length.out = 3),
                                    cex = 1,
                                    col = "black",
                                    labels = sprintf("%.1f", seq(0, 1, length.out = 3))
                                  )
                                ))

  print(heatmap)
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.Sap <- function(x,
                             segment_type = c("motifs", "bouts"),
                             sample_percent = NULL,
                             balanced = FALSE,
                             labels = NULL,
                             clusters = NULL,
                             cores = NULL,
                             seed = 222,
                             msmooth = c(256,50),
                             color_palette = NULL,
                             n_colors = 500,
                             contrast = 3,
                             ordered = FALSE,
                             descending = TRUE,
                             padding_quantile = 0.9,
                             window = NULL,
                             reference_lines = NULL,
                             reference_line_color = "white",
                             ylabel = "Labels",
                             verbose = TRUE,
                             ...) {
  if(verbose) message(sprintf("\n=== Starting Heatmap Plotting ===\n"))

  # Input validation
  segment_type <- match.arg(segment_type)

  # Validate window argument
  if (!is.null(window)) {
    if (segment_type != "bouts") {
      warning("'window' argument is only used when segment_type = 'bouts'. Ignoring...")
      window <- NULL
    } else {
      if (!is.numeric(window) || length(window) != 2) {
        stop("'window' must be a numeric vector of length 2")
      }
      if (any(window <= 0)) {
        stop("'window' values must be positive")
      }
    }
  }

  # Validate reference_lines argument
  if (!is.null(reference_lines)) {
    if (!is.numeric(reference_lines)) {
      stop("'reference_lines' must be a numeric vector")
    }

    # Validate and prepare reference_line_color
    if (length(reference_line_color) == 1) {
      # Single color for all lines
      reference_line_color <- rep(reference_line_color, length(reference_lines))
    } else if (length(reference_line_color) != length(reference_lines)) {
      # Recycle colors if different length
      reference_line_color <- rep(reference_line_color, length.out = length(reference_lines))
      if (verbose) {
        message(sprintf("Note: Recycling %d colors for %d reference lines",
                        length(unique(reference_line_color)), length(reference_lines)))
      }
    }
  }

  # Special handling for motifs
  if (segment_type == "motifs") {
    # Use original motifs if no ordering/cluster filtering needed
    if (!ordered && is.null(clusters)) {
      segments_df <- x[["motifs"]]

      # Use feature embeddings when ordering or clusters requested
    } else {
      if (is.null(x$features$motif$feat.embeds)) {
        stop("Feature embeddings required for ordered/clustered motif plots")
      }

      segments_df <- x$features$motif$feat.embeds |>
        as_segment()

      # Apply UMAP-based ordering if requested
      if (ordered) {
        segments_df <- segments_df |>
          dplyr::arrange(
            day_post_hatch,
            if (descending) dplyr::desc(UMAP2) else UMAP2,
            if (descending) dplyr::desc(UMAP1) else UMAP1
          )
      }
    }

    # For non-motif segment types
  } else {
    segments_df <- x[[segment_type]]
  }

  # Validation
  if (!inherits(segments_df, "segment") || nrow(segments_df) == 0) {
    stop("No segments found in the specified segment type")
  }

  # Select and balance segments
  segments_df <- select_segments(segments_df,
                                 labels = labels,
                                 clusters = clusters,
                                 balanced = balanced,
                                 sample_percent = sample_percent,
                                 seed = seed)

  # Check if segments_df is empty after subsetting
  if (nrow(segments_df) == 0) {
    stop("No segments remaining after subsetting. Check labels/clusters.")
  }

  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }

  # Function to process a single row
  if (segment_type == "motifs") {
    process_row <- function(i) {
      amp_env(segments_df[i,],
              wav_dir = x$base_path,
              msmooth = msmooth)
    }
    time_window <- segments_df$duration[1]
    alignment_point <- segments_df$detection_time[1] - segments_df$start_time[1]

    # Get sampling rate for motifs
    test_env <- amp_env(segments_df[1,], wav_dir = x$base_path, msmooth = msmooth)
    samples_per_second <- length(test_env) / segments_df$duration[1]

  } else if (segment_type == "bouts") {
    required_cols <- c("start_time", "end_time", "align_time")
    if (!all(required_cols %in% colnames(segments_df))) {
      stop("Bouts require: start_time, end_time, align_time columns")
    }

    # Calculate windows using quantiles (existing method)
    segments_df <- segments_df |>
      mutate(
        pre_window = .data$align_time - .data$start_time,
        post_window = .data$end_time - .data$align_time
      ) |>
      arrange(.data$day_post_hatch, .data$duration)

    # Determine max_pre_window and max_post_window
    if (!is.null(window)) {
      # Use provided window values
      max_pre_window <- window[1]
      max_post_window <- window[2]

      if (verbose) {
        message(sprintf("Using provided window: pre=%.3f, post=%.3f seconds",
                        max_pre_window, max_post_window))
      }
    } else {
      # Calculate windows using quantiles (existing method)
      max_pre_window <- quantile(segments_df$pre_window, padding_quantile, na.rm = TRUE)
      max_post_window <- quantile(segments_df$post_window, padding_quantile, na.rm = TRUE)

      if (verbose) {
        message(sprintf("Auto-calculated window: pre=%.3f, post=%.3f seconds (quantile=%.2f)",
                        max_pre_window, max_post_window, padding_quantile))
      }
    }

    # Get sampling rate
    test_row <- segments_df[1, ]
    test_env <- amp_env(test_row, wav_dir = x$base_path, msmooth = msmooth)
    samples_per_second <- length(test_env) / test_row$duration

    # Create processing function
    process_row <- function(i) {
      pad_amp_env(
        segment_row = segments_df[i, ],
        wav_dir = x$base_path,
        msmooth = msmooth,
        max_pre_window = max_pre_window,
        max_post_window = max_post_window,
        samples_per_second = samples_per_second
      )
    }
    time_window <- max_pre_window + max_post_window
    alignment_point <- max_pre_window   # Alignment point for bouts (in seconds from start)
  }

  # Generate amplitude envelope matrix using parallel processing
  if (cores > 1) {
    if (Sys.info()["sysname"] == "Darwin") {#"Linux"
      amp_list <- pbmcapply::pbmclapply(
        seq_len(nrow(segments_df)),
        process_row,
        mc.cores = cores,
        mc.preschedule = FALSE
      )
      amp_matrix <- do.call(cbind, amp_list)
    } else {
      amp_matrix <- pbapply::pbsapply(
        seq_len(nrow(segments_df)),
        process_row,
        cl = cores,
        simplify = TRUE)
    }
  } else {
    amp_matrix <- pbapply::pbsapply(
      seq_len(nrow(segments_df)),
      process_row,
      simplify = TRUE)
  }

  # Check for non-finite values and handle
  if (any(!is.finite(amp_matrix))) {
    warning("amp_matrix contains non-finite values. Replacing with 0.")
    amp_matrix[!is.finite(amp_matrix)] <- 0
  }

  # Convert reference lines from seconds to sample positions if provided
  reference_line_positions <- NULL
  final_reference_colors <- NULL
  if (!is.null(reference_lines) && !is.null(samples_per_second)) {
    # Convert time positions relative to alignment point to absolute sample positions
    reference_line_positions <- (alignment_point + reference_lines) * samples_per_second

    # Filter out positions that would be outside the matrix
    valid_positions <- reference_line_positions >= 1 & reference_line_positions <= nrow(amp_matrix)
    reference_line_positions <- reference_line_positions[valid_positions]
    final_reference_colors <- reference_line_color[valid_positions]

    if (length(reference_line_positions) < length(reference_lines) && verbose) {
      message(sprintf("Note: %d reference line(s) outside the display window were removed",
                      length(reference_lines) - length(reference_line_positions)))
    }

    if (verbose && length(reference_line_positions) > 0) {
      color_info <- paste(paste0(reference_lines[valid_positions], "s (", final_reference_colors, ")"),
                          collapse = ", ")
      message(sprintf("Reference lines at time positions: %s (relative to alignment)", color_info))
    }
  }

  # Store attributes in amp_matrix
  attr(amp_matrix, "msmooth") <- msmooth
  # Store segments_df in attributes for future reference
  attr(amp_matrix, "segments_df") <- segments_df
  attr(amp_matrix, "segment_type") <- segment_type
  attr(amp_matrix, "time_window") <- time_window
  if (!is.null(window) && segment_type == "bouts") {
    attr(amp_matrix, "window") <- window
  }

  # Set column names as labels
  colnames(amp_matrix) <- segments_df$label

  # # Store label order based on order in segments_df
  # label_order <- unique(segments_df$label)
  # attr(amp_matrix, "label_order") <- label_order

  # Update SAP object's features
  feature_type <- sub("s$", "", segment_type) # Remove 's' from end
  x$features[[feature_type]][["amp_env"]] <- amp_matrix

  #  Guidance on accessing the amplitude envelope matrix
  if(verbose) {
    message(sprintf("\nAccess amplitude envelope matrix via: x$features$%s$amp_env", feature_type))
    message(sprintf("Access attributes via: attributes(x$features$%s$amp_env)", feature_type))
  }

  # Create heatmap using matrix method
  plot_heatmap(amp_matrix,
               labels = labels,
               color_palette = color_palette,
               n_colors = n_colors,
               contrast = contrast,
               reference_lines = reference_line_positions,
               reference_line_color = final_reference_colors,
               main = paste("Heatmap of", segment_type),
               ylabel = ylabel,
               ...)

  invisible(x)
}


#' Pad Amplitude Envelope for Alignment
#'
#' @description
#' Internal function to pad and align amplitude envelopes for bout analysis.
#'
#' @param segment_row Single row of segment data
#' @param wav_dir Directory containing WAV files
#' @param msmooth Smoothing parameters
#' @param max_pre_window Maximum pre-window time
#' @param max_post_window Maximum post-window time
#' @param samples_per_second Sampling rate
#' @param pad_value Value to use for padding
#' @param ... Additional arguments
#'
#' @return
#' Padded and aligned amplitude envelope vector
#'
#' @keywords internal
pad_amp_env <- function(segment_row,
                        wav_dir = NULL,
                        msmooth = NULL,
                        max_pre_window = NULL,
                        max_post_window = NULL,
                        samples_per_second = NULL,
                        pad_value = 0,
                        ...) {

  # Extract original amplitude envelope
  orig_env <- amp_env(segment_row = segment_row,
                      wav_dir = wav_dir,
                      msmooth = msmooth,
  )

  # Calculate pre and post windows in samples
  pre_samples <- floor(segment_row$pre_window * samples_per_second)
  post_samples <- floor(segment_row$post_window * samples_per_second)

  # Calculate max allowed samples based on quantile windows
  max_pre_samples <- floor(max_pre_window * samples_per_second)
  max_post_samples <- floor(max_post_window * samples_per_second)

  # Handle pre-window
  if (pre_samples > max_pre_samples) {
    # Truncate
    start_idx <- pre_samples - max_pre_samples + 1
    pre_env <- orig_env[start_idx:pre_samples]
    pre_pad <- numeric(0)  # No padding needed
  } else {
    # Pad
    pre_env <- orig_env[1:pre_samples]
    pre_pad <- rep(pad_value, max_pre_samples - pre_samples)
  }

  # Handle post-window
  if (post_samples > max_post_samples) {
    # Truncate
    post_env <- orig_env[(pre_samples + 1):(pre_samples + max_post_samples)]
    post_pad <- numeric(0)  # No padding needed
  } else {
    # Pad
    post_env <- orig_env[(pre_samples + 1):(pre_samples + post_samples)]
    post_pad <- rep(pad_value, max_post_samples - post_samples)
  }

  # Combine all parts
  padded_env <- c(pre_pad, pre_env, post_env, post_pad)

  return(padded_env)
}



