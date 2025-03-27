# Pitch Module

# Fundamental Frequency Estimation ----------------------------------------------
# Update date : Mar. 26, 2025

#' Fundamental Frequency Estimation
#'
#' Estimates the fundamental frequency (F0) of an audio segmentusing either cepstral analysis or YIN method
#' (currently under development) for a specified audio segment.
#'
#' @param segment_row A single-row data frame containing segment information
#' @param wav_dir Directory containing wav files. If NULL, attempts to use an attribute from the input
#' @param method Pitch estimation method. Currently supports "cepstrum"
#' @param wl Window length for spectral analysis (default: 512)
#' @param ovlp Overlap percentage between windows (default: 80)
#' @param fmax Maximum frequency to consider (default: 1400 Hz)
#' @param threshold Threshold Amplitude threshold for cepstral method in % (default = 10).
#' @param plot Logical, whether to plot the pitch estimation (default: FALSE)
#'
#' @return A matrix with two columns:
#'   - First column: Time points
#'   - Second column: Fundamental frequency values
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'segment' is a data frame with start_time and end_time columns
#' segment_row <- segment[1, ]
#' F0 <- FF(segment_row, wav_dir = "/path/to/wav/files")
#'
#' # Plot the fundamental frequency
#' FF(segment_row, plot = TRUE)
#' }
#'
#' @seealso
#' \code{\link[seewave]{fund}} for underlying pitch detection method
FF <- function(segment_row,
               wav_dir = NULL,
               method = c("cepstrum", "yin"),
               wl = 512,
               ovlp = 80,
               fmax = 1400,
               threshold = 10,
               plot = FALSE) {

  # Check if input is valid
  if (!is.data.frame(segment_row) || nrow(segment_row) != 1) {
    stop("segment_row must be a single row from a data frame")
  }

  # Check for required timing columns
  if (!all(c("start_time", "end_time") %in% names(segment_row))) {
    stop("Input row must contain 'start_time' and 'end_time' columns")
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

  # Properly match the method argument
  method <- match.arg(method)

  # Process based on method
  if(method == "cepstrum") {
    F0 <- seewave::fund(wv, wl = wl, ovlp = ovlp, fmax = fmax,
                        threshold = threshold, plot = plot)

    # Replace NA values with 0 in frequency column
    if (!is.null(F0) && is.matrix(F0) && ncol(F0) >= 2) {
      F0[,2][is.na(F0[,2])] <- 0
    }
  } else if(method == "yin") {
    stop(
      "The YIN pitch detection method is not yet fully implemented.\n",
      "  * Please use method = 'cepstrum' for current analyses\n")
  }

  return(F0)
}

# Plot Heatmap for Fundamental Frequency ----------------------------------------------
# Update date : Mar. 26, 2025

#' Plot Heatmap of Fundamental Frequency
#'
#' @description
#' Creates heatmap visualizations of fundamental frequency from audio segments, supporting
#' multiple data sources and visualization options.
#'
#' @param x An object to visualize (data frame, SAP object, or matrix)
#' @param wav_dir For default method: Path to WAV files directory
#' @param wl Window length for spectral analysis (default: 512)
#' @param ovlp Overlap percentage between windows (default: 50)
#' @param fmax Maximum frequency to consider (default: 1400 Hz)
#' @param threshold Threshold for pitch detection (default: 10)
#' @param plot_freq_lim Optional two-element vector to limit frequency range in plot
#' @param color_palette Function generating color palette
#' @param n_colors Number of colors in heatmap (default: 500)
#' @param segment_type For SAP objects: Type of segments ('motifs')
#' @param sample_percent For SAP objects: Percentage to sample
#' @param balanced For SAP objects: Balance across labels
#' @param clusters Numeric vector of cluster IDs to filter
#' @param labels Optional vector of labels to include
#' @param cores Number of processing cores
#' @param seed Random seed (default: 222)
#' @param ordered For SAP objects: Order by embeddings
#' @param descending For SAP objects: Direction of ordering
#' @param verbose Print progress messages
#' @param main For matrix method: Plot title
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' Fundamental frequency estimation methods:
#' \itemize{
#'   \item Uses cepstrum-based pitch detection
#'   \item Normalizes time series across renditions
#'   \item Supports optional frequency range limiting
#' }
#'
#' For data frames:
#' \itemize{
#'   \item Requires columns: filename, start_time, end_time
#'   \item Extracts fundamental frequency for each segment
#'   \item Creates matrix of aligned F0 values
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Currently supports only 'motifs' segment type
#'   \item Optional balanced sampling
#'   \item Parallel processing support
#'   \item Ordering by feature embeddings
#' }
#'
#' For matrices:
#' \itemize{
#'   \item Direct visualization of pre-computed F0 data
#'   \item Label-based organization
#'   \item Visual separation between groups
#' }
#'
#' @return
#' For default method: List containing F0 matrix, reference time, and plot
#' For SAP objects: Updated object with fundamental frequency features
#' For matrices: Lattice plot object
#'
#' @examples
#' \dontrun{
#' # Basic usage with data frame
#' plot_FundFreq(segments, wav_dir = "path/to/wavs")
#'
#' # Limit frequency range
#' plot_FundFreq(segments,
#'               wav_dir = "path/to/wavs",
#'               plot_freq_lim = c(0.5, 1.5))
#'
#' # SAP object with options
#' plot_FundFreq(sap_obj,
#'               ordered = TRUE,
#'               descending = FALSE)
#'
#' # Matrix with specific labels and color palette
#' plot_FundFreq(f0_matrix,
#'               labels = c("a", "b"),
#'               color_palette = custom_palette)
#' }
#'
#' @seealso
#' \code{\link{FF}} for fundamental frequency estimation
#'
#' @rdname plot_FundFreq
#' @export
plot_FundFreq <- function(x, ...) {
  UseMethod("plot_FundFreq")
}

#' @rdname plot_FundFreq
#' @export
plot_FundFreq.default <- function(x,
                                  wav_dir = NULL,
                                  wl = 512,
                                  ovlp = 50,
                                  fmax = 1400,
                                  threshold = 10,
                                  plot_freq_lim = NULL,
                                  color_palette = NULL,
                                  n_colors = 500,
                                  cores = NULL,
                                  ...) {

  # Validate input structure
  required_cols <- c("filename", "start_time", "end_time")
  if (!all(required_cols %in% names(x))) {
    stop(paste("Missing required columns:", paste(setdiff(required_cols, names(x)), collapse = ", ")))
  }
  if (nrow(x) == 0) stop("Input data frame is empty")

  # Handle audio directory
  if (is.null(wav_dir) && is.null(attr(x, "wav_dir"))) {
    stop("wav_dir must be provided either as argument or attribute")
  }
  wav_dir <- wav_dir %||% attr(x, "wav_dir")

  # Parallel processing of F0 extraction
  f0_list <- parallel_apply(
    indices = 1:nrow(x),
    FUN = function(i) {
      FF_result <- FF(
        segment_row = x[i,],
        wav_dir = wav_dir,
        wl = wl,
        ovlp = ovlp,
        fmax = fmax,
        threshold = threshold,
        plot = FALSE
      )
      list(time = FF_result[,1], f0 = FF_result[,2])
    },
    cores = cores
  )

  # Extract time and frequency components
  time_list <- lapply(f0_list, function(x) x$time)
  f0_list <- lapply(f0_list, function(x) x$f0)

  # Create reference time axis from first rendition
  ref_time <- time_list[[1]]
  ref_duration <- round(max(ref_time) - min(ref_time), 2)

  # Normalize all F0 vectors to reference time scale
  f0_matrix <- sapply(f0_list, function(f0) {
    if(length(f0) != length(ref_time)) {
      approx(
        x = seq(0, 1, length.out = length(f0)),
        y = f0,
        xout = seq(0, 1, length.out = length(ref_time))
      )$y
    } else {
      f0
    }
  })

  # Modify F0 matrix capping
  if (!is.null(plot_freq_lim)) {
    f0_matrix[f0_matrix < plot_freq_lim[1]] <- plot_freq_lim[1]
    f0_matrix[f0_matrix > plot_freq_lim[2]] <- plot_freq_lim[2]

    # Update range for color breaks
    f0_range <- plot_freq_lim
  } else {
    # Only calculate f0_range if plot_freq_lim is not provided
    f0_range <- round(range(f0_matrix, na.rm = FALSE), 2)
  }

  # Configure color mapping
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("black", "darkblue","blue","white", "yellow","orange", "red"))
  }

  # Create color breaks based on the determined f0_range
  color_breaks <- seq(f0_range[1], f0_range[2], length.out = n_colors)

  # Configure visualization parameters
  time_axis <- list(
    at = seq(1, length(ref_time), length.out = 5),
    labels = sprintf("%.1f", seq(0, ref_duration, length.out = 5))
  )

  # Generate heatmap plot
  heatmap <- lattice::levelplot(
    f0_matrix,
    col.regions = color_palette(n_colors),
    at = color_breaks,
    xlab = "Time (s)",
    ylab = "Rendition",
    main = paste("Fundamental Frequency Heatmap\n",
                 "Window:", wl, "samples | Overlap:", ovlp, "%"),
    scales = list(
      x = list(
        at = time_axis$at,
        labels = time_axis$labels,
        rot = 45
      ),
      y = list(draw = FALSE)
    ),
    aspect = "fill",
    colorkey = list(
      space = "right",
      side =4,
      width = 1.5,           # Width of color bar (default=1)
      height = 0.75,         # Height of color bar
      padding.text = 2,
      labels = list(
        at = seq(f0_range[1], f0_range[2], length.out = 3),  # Align with color_breaks
        labels = round(seq(f0_range[1], f0_range[2], length.out = 3),2),
        cex = 0.8
      ),
      title = "Frequency (kHz)",
      title.gpar = list(
        fontsize = 8,
        fontface = "plain",
        lineheight = 0.8
      )
    ),
    ...
  )


  print(heatmap)

  # Return structure includes original times for reference
  structure(
    list(
      f0_matrix = f0_matrix,
      reference_time = ref_time,
      original_times = time_list,
      plot = heatmap
    ),
    class = c("f0_heatmap", "list")
  )
}

#' @rdname plot_FundFreq
#' @export
plot_FundFreq.matrix <- function(x,
                                 labels = NULL,
                                 plot_freq_lim = NULL,
                                 color_palette = NULL,
                                 n_colors = 500,
                                 main = "Fundamental Frequency Heatmap",
                                 ...) {
  # Check that input is a matrix
  if (!is.matrix(x)) {
    stop("Input must be a matrix")
  }

  # Require a 'time_window' attribute for the x-axis scaling
  if (is.null(attr(x, "time_window"))) {
    stop("Matrix must have a 'time_window' attribute")
  }
  time_window <- attr(x, "time_window")
  if (!is.numeric(time_window) || length(time_window) != 1) {
    stop("The 'time_window' attribute must be a single numeric value")
  }

  # Check that the matrix has column names to use as rendition labels
  if (is.null(colnames(x))) {
    stop("Matrix must have column names (to be used as labels)")
  }

  # If a subset of labels is provided, subset the matrix accordingly
  # and preserve the time_window attribute
  if (!is.null(labels)) {
    if (!all(labels %in% colnames(x))) {
      stop("Some provided labels are not found in the matrix's column names")
    }
    x <- x[, colnames(x) %in% labels, drop = FALSE]
    attr(x, "time_window") <- time_window
  }

  # Modify F0 matrix capping
  if (!is.null(plot_freq_lim)) {
    x[x < plot_freq_lim[1]] <- plot_freq_lim[1]
    x[x > plot_freq_lim[2]] <- plot_freq_lim[2]

    # Update range for color breaks
    f0_range <- plot_freq_lim
  } else {
    # Only calculate f0_range if plot_freq_lim is not provided
    f0_range <- round(range(x, na.rm = TRUE), 2)
  }

  # Reverse the matrix so that the first rendition appears on top
  reversed_f0_matrix <- x[, ncol(x):1, drop = FALSE]

  # Reverse the column labels and compute positions for labeling and horizontal lines
  reversed_labels <- rev(unique(colnames(x)))
  samples_per_label <- rev(as.numeric(table(colnames(x))))
  cumulative_positions <- cumsum(c(0, head(samples_per_label, -1)))
  label_positions <- cumulative_positions + samples_per_label / 2
  hline_positions <- cumsum(samples_per_label)[-length(samples_per_label)]

  # Set a default color_palette if none given
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("black", "darkblue", "blue", "white", "yellow", "orange", "red"))
  }

  # Create color breaks based on the determined f0_range
  color_breaks <- seq(f0_range[1], f0_range[2], length.out = n_colors)

  # Define the x-axis (time) scale
  x_scale <- list(
    at = seq(0, nrow(reversed_f0_matrix), length.out = 5),
    labels = sprintf("%.1f", seq(0, time_window, length.out = 5))
  )

  # Create the heatmap using levelplot
  heatmap <- lattice::levelplot(
    reversed_f0_matrix,
    col.regions = color_palette(n_colors),
    at = color_breaks,
    border = "transparent",
    xlab = "Time (s)",
    ylab = "Rendition",
    main = main,
    scales = list(
      x = x_scale,
      y = list(
        at = label_positions,
        labels = reversed_labels
      )
    ),
    panel = function(x, y, z, ...) {
      lattice::panel.levelplot(x, y, z, ...)
      lattice::panel.abline(h = hline_positions, col = "white", lwd = 2)
    },
    aspect = "fill",
    colorkey = list(
      space = "right",
      width = 1,
      height = 0.75,
      labels = list(
        at = seq(f0_range[1], f0_range[2], length.out = 3),
        labels = sprintf("%.2f", seq(f0_range[1], f0_range[2], length.out = 3)),
        cex = 1,
        col = "black"
      ),
      title = "Frequency (kHz)",
      title.gpar = list(
        fontsize = 8,
        fontface = "plain",
        lineheight = 0.8
      )
    ),
    ...
  )

  print(heatmap)

  invisible(heatmap)
}

#' @rdname plot_FundFreq
#' @export
plot_FundFreq.Sap <- function(x,
                              segment_type = c("motifs", "syllables", "segments"),
                              sample_percent = NULL,
                              balanced = FALSE,
                              labels = NULL,
                              clusters = NULL,
                              cores = NULL,
                              seed = 222,
                              wl = 512,
                              ovlp = 50,
                              fmax = 1400,
                              threshold = 10,
                              plot_freq_lim = NULL,
                              color_palette = NULL,
                              n_colors = 500,
                              ordered = FALSE,
                              descending = TRUE,
                              verbose = TRUE,
                              ...) {
  if(verbose) message(sprintf("\n=== Starting Fundamental Frequency Heatmap Plotting ===\n"))

  # Input validation
  segment_type <- match.arg(segment_type)

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
  } else {
    stop("Currently only 'motifs' segment type is supported for fundamental frequency heatmap")
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

  # Parallel processing of F0 extraction
  f0_list <- parallel_apply(
    indices = 1:nrow(segments_df),
    FUN = function(i) {
      FF_result <- FF(
        segment_row = segments_df[i,],
        wav_dir = x$base_path,
        wl = wl,
        ovlp = ovlp,
        fmax = fmax,
        threshold = threshold,
        plot = FALSE
      )
      list(time = FF_result[,1], f0 = FF_result[,2])
    },
    cores = cores
  )

  # Extract time and frequency components
  time_list <- lapply(f0_list, function(x) x$time)
  f0_list <- lapply(f0_list, function(x) x$f0)

  # Create reference time axis from first rendition
  ref_time <- time_list[[1]]
  time_window <- round(max(ref_time) - min(ref_time), 2)

  # Normalize all F0 vectors to reference time scale
  f0_matrix <- sapply(f0_list, function(f0) {
    if(length(f0) != length(ref_time)) {
      approx(
        x = seq(0, 1, length.out = length(f0)),
        y = f0,
        xout = seq(0, 1, length.out = length(ref_time))
      )$y
    } else {
      f0
    }
  })

  # Store attributes in f0_matrix
  attr(f0_matrix, "wl") <- wl
  attr(f0_matrix, "ovlp") <- ovlp
  attr(f0_matrix, "segments_df") <- segments_df
  attr(f0_matrix, "segment_type") <- segment_type
  attr(f0_matrix, "time_window") <- time_window

  # Set column names as labels
  colnames(f0_matrix) <- segments_df$label

  # Update SAP object's features
  x$features$motif[["fund_freq"]] <- f0_matrix

  if(verbose) {
    message(sprintf("\nAccess fundamental frequency matrix via: x$features$motif$fund_freq"))
    message(sprintf("Access attributes via: attributes(x$features$motif$fund_freq)"))
  }

  # Create heatmap using matrix method
  plot_FundFreq(f0_matrix,
                labels = labels,
                plot_freq_lim = plot_freq_lim,
                color_palette = color_palette,
                n_colors = n_colors,
                main = "Heatmap of Motifs: Fundamental Frequency")

  invisible(x)
}
