# Amplitude Envelope module


# Extract Amplitude Envelope ----------------------------------------------
# Update date : Feb. 4, 2025

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
#' @param norm Logical, whether to normalize envelope to range [0,1] (default: FALSE)
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
#' @importFrom tuneR readWave
#' @importFrom seewave env
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
# Update date : Feb. 4, 2025

#' Plot Heatmap of Amplitude Envelopes
#'
#' @description
#' A generic function to create heatmap visualizations of amplitude envelopes
#' from audio segments.
#'
#' @param x An object to visualize: data frame, SAP object, or amplitude matrix
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' This generic function supports heatmap creation through three methods:
#' \itemize{
#'   \item Default method for segment data frames
#'   \item SAP object method for organized song data
#'   \item Matrix method for pre-computed amplitude envelopes
#' }
#'
#' @return
#' A heatmap visualization or updated SAP object
#'
#' @examples
#' \dontrun{
#' # Plot from segment data frame
#' plot_heatmap(segments, wav_dir = "path/to/wavs")
#'
#' # Plot from SAP object
#' plot_heatmap(sap_obj, segment_type = "motifs")
#'
#' # Plot from amplitude matrix
#' plot_heatmap(amp_matrix)
#' }
#'
#' @export
plot_heatmap <- function(x, ...) {
  UseMethod("plot_heatmap")
}

#' Plot Heatmap from Segment Data Frame
#'
#' @description
#' Creates a heatmap visualization from a data frame containing segment information.
#'
#' @param x A data frame with segment information
#' @param wav_dir Path to WAV files directory
#' @param msmooth Numeric vector of length 2 for envelope smoothing:
#'        c(window_length, overlap_percentage)
#' @param color_palette Function generating color palette
#' @param n_colors Number of colors in heatmap
#' @param contrast Contrast factor for visualization
#' @param ... Additional arguments
#'
#' @details
#' Creates a heatmap visualization with the following steps:
#' \itemize{
#'   \item Validates input segment data
#'   \item Calculates amplitude envelopes for each segment
#'   \item Creates matrix of aligned envelopes
#'   \item Generates heatmap visualization
#' }
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item segments: Input segment data frame
#'   \item amp_matrix: Matrix of amplitude envelopes
#'   \item plot: Heatmap plot object
#' }
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pbsapply
#' @importFrom lattice levelplot
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
  return(result)
}

#' Plot Heatmap from Amplitude Matrix
#'
#' @description
#' Creates a heatmap visualization from a pre-computed amplitude envelope matrix.
#'
#' @param x Matrix of amplitude envelopes with time_window attribute
#' @param labels Optional vector of labels to subset the matrix
#' @param color_palette Function generating color palette
#' @param n_colors Number of colors in heatmap
#' @param contrast Contrast factor for visualization
#' @param main Title for the heatmap
#' @param ... Additional arguments
#'
#' @details
#' Creates a heatmap visualization with the following features:
#' \itemize{
#'   \item Label-based organization
#'   \item Customizable color scheme
#'   \item Automatic time scaling
#'   \item Visual separation between labels
#' }
#'
#' @return
#' A lattice plot object containing the heatmap visualization
#'
#' @examples
#' \dontrun{
#' # Plot full matrix
#' plot_heatmap(sap_obj$features$motif$amp_env)
#'
#' # Plot subset of labels
#' plot_heatmap(sap_obj$features$motif$amp_env,
#'              labels = c("a", "b"))
#' }
#'
#' @export
plot_heatmap.matrix <- function(x,
                                labels = NULL,
                                color_palette = NULL,
                                n_colors = 500,
                                contrast = 3,
                                main = "Amplitude Envelope Heatmap",
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

  # Get reversed labels and calculate positions
  reversed_labels <- rev(unique(colnames(x)))
  samples_per_label <- rev(table(colnames(x)))
  cumulative_positions <- cumsum(c(0, head(samples_per_label, -1)))
  label_positions <- cumulative_positions + samples_per_label/2
  hline_positions <- cumsum(samples_per_label)[-length(samples_per_label)]

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
                                ylab = "Labels",
                                xlab = "Time (s)",
                                main = main,
                                scales = list(
                                  x = list(
                                    #at = seq(0, ncol(reversed_amp_matrix), length.out = 5),
                                    at = pretty(1:nrow(reversed_amp_matrix), n = 5),
                                    labels = sprintf("%.1f", seq(0, time_window, length.out = 5))
                                  ),
                                  y = list(
                                    at = label_positions,
                                    labels = reversed_labels
                                  )
                                ),
                                panel = function(x, y, z, ...) {
                                  lattice::panel.levelplot(x, y, z, ...)
                                  lattice::panel.abline(h = hline_positions,
                                                        col = "white",
                                                        lwd = 2)
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

#' Plot Heatmap from SAP Object
#'
#' @description
#' Creates a heatmap visualization from segments stored in a SAP object,
#' with options for sampling and organization.
#'
#' @param x A SAP object containing song recordings
#' @param segment_type Type of segments to analyze ("motifs", "bouts", "syllables", "segments")
#' @param sample_percent Percentage of segments to sample from each label
#' @param balanced Whether to balance groups across labels
#' @param labels Specific labels to include
#' @param cores Number of cores for parallel processing
#' @param seed Random seed for sampling
#' @param msmooth Numeric vector for envelope smoothing: c(window_length, overlap_percentage)
#' @param color_palette Function generating color palette
#' @param n_colors Number of colors in heatmap
#' @param contrast Contrast factor for visualization
#' @param ordered Whether to order segments based on feature embeddings
#' @param descending Direction of ordering
#' @param padding_quantile Quantile for padding in bout analysis
#' @param verbose Whether to print progress messages
#' @param ... Additional arguments
#'
#' @details
#' Creates a heatmap with the following features:
#' \itemize{
#'   \item Supports multiple segment types
#'   \item Optional balanced sampling across labels
#'   \item Parallel processing for large datasets
#'   \item Ordering based on feature embeddings
#'   \item Automatic padding for bout analysis
#' }
#'
#' @return
#' Updated SAP object with amplitude envelope features and visualization
#'
#' @importFrom dplyr group_by summarise arrange ungroup slice_sample n count pull group_modify
#' @importFrom lattice levelplot panel.levelplot panel.abline
#' @export
plot_heatmap.Sap <- function(x,
                             segment_type = c("motifs", "bouts", "syllables", "segments"),
                             sample_percent = NULL,
                             balanced = FALSE,
                             labels = NULL,
                             cores = NULL,
                             seed = 222,
                             msmooth = c(256,50),
                             color_palette = NULL,
                             n_colors = 500,
                             contrast = 3,
                             ordered = FALSE,
                             descending = TRUE,
                             padding_quantile = 0.9,
                             verbose = TRUE,
                             ...) {
  if(verbose) message(sprintf("\n=== Starting Heatmap Plotting ===\n"))

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
  } else if (segment_type == "bouts") {
    required_cols <- c("start_time", "end_time", "align_time")
    if (!all(required_cols %in% colnames(segments_df))) {
      stop("Bouts require: start_time, end_time, align_time columns")
    }

    # Calculate windows using quantiles
    segments_df <- segments_df |>
      mutate(
        pre_window = align_time - start_time,
        post_window = end_time - align_time
      ) |>
      arrange(day_post_hatch, duration)

    max_pre_window <- quantile(segments_df$pre_window, padding_quantile, na.rm = TRUE)
    max_post_window <- quantile(segments_df$post_window, padding_quantile, na.rm = TRUE)

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

  # Store attributes in amp_matrix
  attr(amp_matrix, "msmooth") <- msmooth
  # Store segments_df in attributes for future reference
  attr(amp_matrix, "segments_df") <- segments_df
  attr(amp_matrix, "segment_type") <- segment_type
  attr(amp_matrix, "time_window") <- round(time_window, 2)

  # Set column names as labels
  colnames(amp_matrix) <- segments_df$label

  # Update SAP object's features
  feature_type <- sub("s$", "", segment_type) # Remove 's' from end
  x$features[[feature_type]][["amp_env"]] <- amp_matrix


  # if order is requested
  if (isTRUE(ordered)) {
    # Check if the feature embeddings 'feat.embeds' exist
    if (!is.null(x$features[[feature_type]][["feat.embeds"]])) {
      embeds <- x$features[[feature_type]][["feat.embeds"]]

      # Create a unique identifier for each segment
      segments_df$segment_id <- paste(segments_df$filename,
                                      segments_df$start_time,
                                      segments_df$end_time,
                                      sep = "_")
      embeds$segment_id <- paste(embeds$filename,
                                 embeds$start_time,
                                 embeds$end_time,
                                 sep = "_")

      # Check if all segments in 'segments_df' are in 'embeds'
      if (all(segments_df$segment_id %in% embeds$segment_id)) {
        # Arrange embeddings
        if (isTRUE(descending)) {
          embeds_arranged <- embeds |>
            dplyr::arrange(day_post_hatch, desc(UMAP2), desc(UMAP1))
        } else {
          embeds_arranged <- embeds |>
            dplyr::arrange(day_post_hatch, UMAP2, UMAP1)
        }

        # Get the new order of segment IDs
        new_order_segment_ids <- embeds_arranged$segment_id

        # Create a mapping from new order to columns in amp_matrix
        mapping <- match(new_order_segment_ids, segments_df$segment_id)

        # Reorder amp_matrix columns
        amp_matrix <- amp_matrix[, mapping]

        # Reorder segments_df accordingly
        segments_df <- segments_df[mapping, ]

        # Update column names of amp_matrix
        colnames(amp_matrix) <- segments_df$label

        # Reset attributes after reordering
        attr(amp_matrix, "time_window") <- round(time_window, 2)

      } else {
        warning("Not all segments in 'segments_df' are present in 'feat.embeds'. Using original order.")
      }
    } else {
      warning("feat.embeds not found in features. Using original order.")
    }
  }

  # Create heatmap using matrix method
  plot_heatmap(amp_matrix,
               labels = labels,
               color_palette = color_palette,
               n_colors = n_colors,
               contrast = contrast,
               main = paste("Heatmap of", segment_type))

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


# Plot Traces -------------------------------------------------------------

#' Plot amplitude envelope traces
#'
#' @param x Amplitude envelope matrix with time_window attribute and label colnames
#' @param labels Optional vector of labels to subset the matrix (default: NULL, uses all labels)
#' @param plot_type Type of plot: "individual", "average", or "combined" (default: "combined")
#' @param alpha Transparency for individual traces (default: 0.2)
#' @param ncol Number of columns for facets (default: 1)
#' @param palette Color palette for plotting (default: "Set1")
#'
#' @return A ggplot object
#'
#' @return A ggplot object
#'
#' @importFrom dplyr mutate select filter arrange left_join bind_rows %>%
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs theme_minimal scale_color_brewer
#'             scale_fill_brewer ggtitle theme element_text stat_summary
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' # Direct use of stored matrix with all labels
#' plot_traces(sap_obj$features$motif$amp_env)
#'
#' # With specific labels
#' plot_traces(sap_obj$features$motif$amp_env,
#'            labels = c("a", "b"),
#'            plot_type = "average")
#' }
#'
#' @export
plot_traces <- function(x,
                        labels = NULL,
                        plot_type = c("combined", "individual", "average"),
                        alpha = 0.2,
                        ncol = 1,
                        palette = "Set1") {

  # Input validation
  plot_type <- match.arg(plot_type)

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

  # Create data frame from matrix and make column names unique
  dat <- as.data.frame(x)
  original_names <- colnames(dat)
  colnames(dat) <- make.unique(colnames(dat))

  # Calculate time points
  time_window <- attr(x, "time_window")
  time_step <- time_window / nrow(x)
  time_points <- seq(0, time_window, by = time_step)[1:nrow(x)]

  # Create mapping dataframe with original and unique names
  mapping_df <- data.frame(
    col_name = colnames(dat),
    label = original_names,
    rendition_no = sequence(rle(original_names)$lengths)
  )

  # Create long format data
  res_long <- dat |>
    dplyr::mutate(time = time_points) |>
    tidyr::pivot_longer(cols = -time, names_to = "col_name", values_to = "amplitude") |>
    dplyr::left_join(mapping_df, by = "col_name") |>
    dplyr::select(-col_name) |>
    dplyr::arrange(rendition_no, time)

  # Filter labels if specified
  if (!is.null(labels)) {
    if (!all(labels %in% unique(res_long$label))) {
      stop("Some provided labels not found in matrix")
    }
    res_long <- res_long |>
      dplyr::filter(label %in% labels)
  }

  # Load mean_se function from ggplot2
  mean_se <- function(x) {
    x <- stats::na.omit(x)
    se <- stats::sd(x) / sqrt(length(x))
    mean <- mean(x)
    data.frame(y = mean,
               ymin = mean - se,
               ymax = mean + se)
  }

  # Create plot based on type
  if (plot_type == "individual") {
    p <- ggplot2::ggplot(res_long,
                         ggplot2::aes(x = time, y = amplitude,
                                      group = rendition_no, color = label)) +
      ggplot2::geom_line(alpha = alpha) +
      ggplot2::facet_wrap(~label, ncol = ncol) +
      ggplot2::labs(x = "Time (s)",
                    y = "Amplitude Envelope",
                    color = "Label") +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_brewer(palette = palette) +
      ggplot2::ggtitle("Amplitude Envelope by Label") +
      ggplot2::theme(legend.position = "right",
                     strip.text = ggplot2::element_text(size = 12, face = "bold"))

  } else if (plot_type == "average") {
    p <- ggplot2::ggplot(res_long,
                         ggplot2::aes(x = time, y = amplitude,
                                      color = label, fill = label)) +
      ggplot2::stat_summary(fun.data = mean_se, geom = "ribbon",
                            alpha = 0.3, color = NA) +
      ggplot2::stat_summary(fun = mean, geom = "line", size = 1) +
      ggplot2::labs(x = "Time (s)", y = "Amplitude Envelope") +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_brewer(palette = palette) +
      ggplot2::scale_fill_brewer(palette = palette) +
      ggplot2::ggtitle("Mean Amplitude Envelope with Standard Error")
  }

  else if (plot_type == "combined") {
    # Create combined data for both plots
    res_long <- res_long |>
      dplyr::mutate(plot_type = "Individual Traces")

    res_long_mean <- res_long |>
      dplyr::mutate(plot_type = "Mean ± SE")

    res_combined <- dplyr::bind_rows(res_long, res_long_mean)

    p <- ggplot2::ggplot(res_combined,
                         ggplot2::aes(x = time, y = amplitude,
                                      color = label, fill = label)) +
      ggplot2::geom_line(data = . %>% dplyr::filter(plot_type == "Individual Traces"),
                         ggplot2::aes(group = interaction(label, rendition_no)),
                         alpha = alpha) +
      ggplot2::stat_summary(data = . %>% dplyr::filter(plot_type == "Mean ± SE"),
                            fun.data = mean_se, geom = "ribbon",
                            alpha = 0.3, color = NA) +
      ggplot2::stat_summary(data = . %>% dplyr::filter(plot_type == "Mean ± SE"),
                            fun = mean, geom = "line", size = 1) +
      ggplot2::facet_wrap(~plot_type, ncol = ncol) +
      ggplot2::labs(x = "Time (s)",
                    y = "Amplitude Envelope",
                    color = "Label",
                    fill = "Label") +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_brewer(palette = palette) +
      ggplot2::scale_fill_brewer(palette = palette) +
      ggplot2::ggtitle("Amplitude Envelope") +
      ggplot2::theme(legend.position = "right",
                     strip.text = ggplot2::element_text(size = 12, face = "bold"))
  } else {
    stop("Invalid plot_type. Choose 'individual', 'average', or 'combined'")
  }

  return(p)
}

