# Pitch Module

# Fundamental Frequency Estimation ----------------------------------------------
# Update date : Apr.2, 2025

#' Fundamental Frequency Analysis and Visualization
#'
#' @description
#' Creates fundamental frequency analyses and visualizations from audio segments, supporting
#' multiple data types (data frames, SAP objects) and visualization options.
#' Provides both analytical results and optional single trial or heatmap visualizations.
#'
#' @param x Input object:
#'   \itemize{
#'     \item data frame (default method)
#'     \item SAP object
#'     \item Pre-computed F0 matrix
#'   }
#' @param wav_dir Directory containing WAV files (for data frame methods)
#' @param method Pitch estimation method ("cepstrum" or "yin")
#' @param wl Window length for spectral analysis (default: 512)
#' @param ovlp Overlap percentage between windows (default: 50)
#' @param fmax Maximum frequency to consider (default: 1400 Hz)
#' @param threshold Amplitude threshold for pitch detection in % (default: 10)
#' @param plot Logical, whether to generate visualization(default: TRUE)
#' @param plot_freq_lim Optional vector of length 2 specifying frequency limits for plotting
#' @param color_palette Function generating color palette (default: black-to-red spectrum)
#' @param n_colors Number of colors in heatmap (default: 500)
#' @param main Plot title (for matrix method)
#' @param labels Optional vector of labels for subsetting or grouping
#' @param segment_type For SAP objects: Type of segments (currently only 'motifs')
#' @param sample_percent For SAP objects: Percentage of segments to sample
#' @param balanced For SAP objects: Whether to balance samples across labels
#' @param clusters For SAP objects: Numeric vector of cluster IDs to include
#' @param cores Number of cores for parallel processing
#' @param seed Random seed for sampling (default: 222)
#' @param ordered For SAP objects: Whether to order by feature embeddings
#' @param descending For SAP objects: Direction of embedding-based ordering
#' @param verbose Logical, whether to print progress messages
#' @param ... Additional arguments passed to methods or lattice::levelplot
#'
#' @details
#' The function provides different methods depending on the input type:
#'
#'
#' Default method (Data frame method)
#' \itemize{
#'   \item Performs fundamental frequency analysis on a single segment
#'          or processes multiple segments in parallel
#'   \item Supports both cepstrum and YIN-based pitch detection
#'   \item Normalizes time series across renditions
#'   \item Creates aligned F0 matrix
#'   \item Optional spectrogram overlay (single segment)
#'         or heatmap visualization(multiple segments)
#' }
#'
#' SAP object method:
#' \itemize{
#'   \item Integrates with SAP feature analysis pipeline
#'   \item Supports segment filtering and ordering
#'   \item Stores results in features component
#' }
#'
#' Matrix method:
#' \itemize{
#'   \item Visualizes pre-computed F0 data
#'   \item Supports label-based organization
#'   \item Adds visual separators between groups
#' }
#'
#' @return
#' Returns an object depending on the method used:
#' \itemize{
#'   \item Default method: Matrix of time and frequency values(F0 matrix)
#'         or List with F0 matrix and metadata
#'   \item SAP method: Updated SAP object with F0 features
#'   \item Matrix method: Lattice plot object (invisibly)
#' }
#'
#' @examples
#' \dontrun{
#' # Single segment analysis
#' FF(segment_row, wav_dir = "path/to/wavs", plot = TRUE)
#'
#' # Multiple segment analysis
#' FF(segments_df, wav_dir = "path/to/wavs",
#'    plot = TRUE, plot_freq_lim = c(0.5, 1.5))
#'
#' # SAP object analysis
#' FF(sap_obj, ordered = TRUE, balanced = TRUE,
#'    sample_percent = 80, seed = 123)
#'
#' # Direct matrix visualization
#' FF(f0_matrix, labels = c("a", "b", "c"),
#'    main = "Custom Title")
#' }
#'
#' @seealso
#' \code{\link[lattice]{levelplot}} for underlying plotting function
#' \code{\link[seewave]{fund}} for cepstral analysis
#' \code{\link[librosa]{yin}} for YIN algorithm (when using method = "yin")
#'
#' @rdname Fundamental_Frequency
#' @export
FF <- function(x, ...) {
  UseMethod("FF")
}

#' @rdname Fundamental_Frequency
#' @export
FF.default <- function(x,
                       wav_dir = NULL,
                       wl = 512,
                       ovlp = 50,
                       fmax = 1400,
                       threshold = 10,
                       method = c("cepstrum", "yin"),
                       plot = TRUE,
                       plot_freq_lim = NULL,
                       color_palette = NULL,
                       n_colors = 500,
                       cores = NULL,
                       ...) {
  # Properly match the method argument
  method <- match.arg(method)

  # Validate input
  if (!is.data.frame(x)) {
    stop("Input must be a data frame")
  }
  if (nrow(x) == 0) {
    stop("Input data frame is empty")
  }

  required_cols <- c("filename", "start_time", "end_time")
  if (!all(required_cols %in% names(x))) {
    stop(paste("Missing required columns:", paste(setdiff(required_cols, names(x)), collapse = ", ")))
  }

  # Handle audio directory
  if (is.null(wav_dir) && is.null(attr(x, "wav_dir"))) {
    stop("wav_dir must be provided either as argument or attribute")
  }
  wav_dir <- wav_dir %||% attr(x, "wav_dir")

  if (nrow(x) == 1) {
    # Call the single row implementation
    return(FF_single_row(segment_row = x,
                         wav_dir =  wav_dir,
                         wl = wl,
                         ovlp = ovlp,
                         fmax = fmax,
                         threshold = threshold,
                         method = method,
                         plot = plot))
  } else {
    # Validate python dependencies if method "yin" is chosen
    if(method == "yin") {

      check_python_dependencies()

      tryCatch({
        librosa <- reticulate::import("librosa")
        np <- reticulate::import("numpy")
      }, error = function(e) {
        stop("Python dependencies not found. Install with:\n",
             "reticulate::py_install(c('librosa', 'numpy'))\n",
             "Then restart R and try again.")
      })
    }

    # Parallel processing of F0 extraction
    f0_list <- parallel_apply(
      indices = 1:nrow(x),
      FUN = function(i) {
        FF_result <- FF_single_row(
          segment_row = x[i,],
          wav_dir = wav_dir,
          wl = wl,
          ovlp = ovlp,
          fmax = fmax,
          threshold = threshold,
          method = method,
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

    # Return structure includes original times for reference
    result <- structure(
      list(
        f0_matrix = f0_matrix,
        reference_time = ref_time,
        original_times = time_list
      ),
      class = c("f0_matrix", "list")
    )

    if (plot) {
      # Validate plot parameters
      if (!is.null(plot_freq_lim)) {
        if (length(plot_freq_lim) != 2 || !is.numeric(plot_freq_lim)) {
          stop("plot_freq_lim must be a numeric vector of length 2")
        }
      }

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
      result$plot <- heatmap
    }

    return(result)

  }
}

#' @rdname Fundamental_Frequency
#' @export
FF.Sap <- function(x,
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
                   method = c("cepstrum", "yin"),
                   plot = TRUE,
                   plot_freq_lim = NULL,
                   color_palette = NULL,
                   n_colors = 500,
                   ordered = FALSE,
                   descending = TRUE,
                   verbose = TRUE,
                   ...) {
  if(verbose) message(sprintf("\n=== Starting Fundamental Frequency Analysis ===\n"))

  # Properly match the method argument
  method <- match.arg(method)

  # Validate python dependencies if method "yin" is chosen
  if(method == "yin") {

    check_python_dependencies()

    tryCatch({
      librosa <- reticulate::import("librosa")
      np <- reticulate::import("numpy")
    }, error = function(e) {
      stop("Python dependencies not found. Install with:\n",
           "reticulate::py_install(c('librosa', 'numpy'))\n",
           "Then restart R and try again.")
    })
  }

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
      FF_result <- FF_single_row(
        segment_row = segments_df[i,],
        wav_dir = x$base_path,
        wl = wl,
        ovlp = ovlp,
        fmax = fmax,
        threshold = threshold,
        method = method,
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

  if(plot) {
  # Create heatmap using matrix method
  FF(f0_matrix,
     labels = labels,
     plot_freq_lim = plot_freq_lim,
     color_palette = color_palette,
     n_colors = n_colors,
     main = "Heatmap of Motifs: Fundamental Frequency")
  }

  invisible(x)
}

#' @rdname Fundamental_Frequency
#' @export
FF.matrix <- function(x,
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

  # Validate plot_freq_lim
  if (!is.null(plot_freq_lim) &&
      (!is.numeric(plot_freq_lim) || length(plot_freq_lim) != 2)) {
    stop("plot_freq_lim must be a numeric vector of length 2")
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


#' Internal function to calculate the fundamental frequency of time series data derived from a single audio segment
#'
#' Estimates the fundamental frequency (F0) of an audio segment using either cepstral analysis
#' or YIN method for a specified audio segment.
#'
#' @param segment_row A single-row data frame containing segment information
#' @param wav_dir Directory containing wav files.
#' @param method Pitch estimation method. Either "cepstrum" or "yin" (default: "cepstrum")
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
#' @seealso
#' \code{\link[seewave]{fund}} for cepstral analysis method
#' \code{\link[librosa]{yin}} for YIN algorithm implementation
#'
#' @keywords internal
FF_single_row  <- function(segment_row,
                   wav_dir = NULL,
                   method = c("cepstrum", "yin"),
                   wl = 512,
                   ovlp = 50,
                   fmax = 1400,
                   threshold = 10,
                   plot = FALSE,
               ...) {
  # Validate python dependencies if method "yin" is chosen
  if(method == "yin") {

    check_python_dependencies()

    tryCatch({
      librosa <- reticulate::import("librosa")
      np <- reticulate::import("numpy")
    }, error = function(e) {
      stop("Python dependencies not found. Install with:\n",
           "reticulate::py_install(c('librosa', 'numpy'))\n",
           "Then restart R and try again.")
    })
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
                        threshold = threshold, plot = FALSE)

    # Replace NA values with 0 in frequency column
    if (!is.null(F0) && is.matrix(F0) && ncol(F0) >= 2) {
      F0[,2][is.na(F0[,2])] <- 0
    }
  } else if(method == "yin") {
    # Convert tuneR wave to numpy array
    audio <- np$array(wv@left)
    sr <- wv@samp.rate

    # Calculate hop length (similar to spectrogram calculation)
    hop_length <- as.integer(wl * (1 - ovlp/100))

    # Compute pitch using librosa YIN
    pitch <- librosa$yin(
      y = audio,
      fmin = 100,  # Minimum frequency (adjust as needed)
      fmax = fmax,
      sr = sr,
      hop_length = hop_length
    )

    # Create time vector
    times <- seq(0,
                 by = hop_length/sr,
                 length.out = length(pitch))

    F0 <- cbind(time = times, frequency = as.numeric(pitch/1000))
  }


  # Plotting
  if(plot) {
    # Generate the spectrogram using the visualize_song.default function
    visualize_song.default(x = sound_path,
                           start_time_in_second = segment_row$start_time,
                           end_time_in_second = segment_row$end_time,
                           fft_window_size = wl,
                           overlap = ovlp/100,
                           dark_mode = TRUE,
                           legend = FALSE,
                           keep.par = TRUE,
                           verbose = FALSE)

    F0[,1]  <-  F0[,1] + segment_row$start_time

    if (!"duration" %in% names(segment_row)) {
      segment_row$duration <- segment_row$end_time - segment_row$start_time
    }
    offset <- segment_row$duration * 0.035 #adjust xlim offset

    # Overlay the F0 trace
    par(new = TRUE)
    plot(F0[,1], F0[,2]*1000,
         type = "p",  # Changed from "l" to "p" for points
         col = "white",
         pch = 16,    # Solid dot
         cex = 0.5,   # Size of the dot
         axes = FALSE,
         xlab = "",
         ylab = "",
         xlim = c(segment_row$start_time + offset, segment_row$end_time - offset),
         ylim = c(900, wv@samp.rate/2))  # manual adjust ylim offset

    # Add a title to the plot
    title(main = "Spectrogram w/ Fundamental Frequency Overlay")
  }

  return(F0)
}

# Goodness of pitch ----------------------------------------------
# Update date : Apr.3, 2025

#' Pitch Goodness Analysis
#'
#' @description
#' Calculates and visualizes pitch goodness for audio segments using cepstral analysis.
#' Supports both single segment analysis and batch processing with visualization options.
#'
#' @param x Input object:
#'   \itemize{
#'     \item  Data frame (default method)
#'     \item SAP object
#'     \item Pre-computed goodness matrix
#'   }
#' @param wav_dir Directory containing WAV files (for data frame methods)
#' @param wl Window length for spectral analysis (default: 512)
#' @param ovlp Overlap percentage between windows (default: 50)
#' @param fmax Maximum frequency to consider (default: 1500 Hz)
#' @param plot Logical, whether to generate visualization (default: TRUE)
#' @param plot_lim Optional vector of length 2 specifying goodness limits for plotting
#' @param color_palette Function generating color palette (default: black-to-white spectrum)
#' @param n_colors Number of colors in heatmap (default: 500)
#' @param main Plot title (for matrix method)
#' @param labels Optional vector of labels for subsetting or grouping
#' @param segment_type For SAP objects: Type of segments (currently only 'motifs')
#' @param sample_percent For SAP objects: Percentage of segments to sample
#' @param balanced For SAP objects: Whether to balance samples across labels
#' @param clusters For SAP objects: Numeric vector of cluster IDs to include
#' @param cores Number of cores for parallel processing
#' @param seed Random seed for sampling (default: 222)
#' @param ordered For SAP objects: Whether to order by feature embeddings
#' @param descending For SAP objects: Direction of embedding-based ordering
#' @param verbose Logical, whether to print progress messages(default: TRUE)
#' @param ... Additional arguments passed to methods or lattice::levelplot
#'
#' @param ... Additional arguments passed to methods or plotting functions
#'
#' @details
#' The function provides different methods depending on the input type:
#'
#' Default method (data frame):
#' \itemize{
#'   \item Processes single or multiple segments
#'   \item Uses cepstral analysis for pitch goodness calculation
#'   \item Supports parallel processing for multiple segments
#'   \item Provides optional visualization
#' }
#'
#' The pitch goodness measure:
#' \itemize{
#'   \item Based on cepstral peak prominence
#'   \item Higher values indicate stronger pitch periodicity
#'   \item Computed across time windows
#'   \item Normalized for comparison across segments
#' }
#'
#' @return
#' For single segments (nrow = 1):
#' \itemize{
#'   \item Matrix with columns: time, goodness
#' }
#'
#' For multiple segments (nrow > 1):
#' \itemize{
#'   \item List with components:
#'     \describe{
#'       \item{goodness_matrix}{Matrix of aligned goodness values}
#'       \item{reference_time}{Vector of reference time points}
#'       \item{original_times}{List of original time vectors}
#'       \item{plot}{If plot=TRUE, the generated lattice plot}
#'     }
#' }
#'
#' @examples
#' \dontrun{
#' # Single segment analysis
#' goodness(segment_row, wav_dir = "path/to/wavs", plot = TRUE)
#'
#' # Multiple segment analysis
#' goodness(segments_df, wav_dir = "path/to/wavs",
#'          plot = TRUE, plot_lim = c(0, 1))
#'
#' # SAP object analysis
#' goodness(sap_obj, ordered = TRUE, balanced = TRUE,
#'    sample_percent = 80, seed = 123)
#'
#' # Direct matrix visualization
#' goodness(goodness_matrix, labels = c("a", "b", "c"),
#'    main = "Custom Title")
#'}
#'
#' @seealso
#' \code{\link{FF}} for fundamental frequency analysis
#' \code{\link[seewave]{spectro}} for underlying spectral analysis
#'
#' @export
goodness <- function(x, ...) {
  UseMethod("goodness")
}

#' @rdname goodness
#' @export
goodness.default <- function(x,
                             wav_dir = NULL,
                             wl = 512,
                             ovlp = 50,
                             fmax = 1500,
                             plot = TRUE,
                             plot_lim = NULL,
                             color_palette = NULL,
                             n_colors = 500,
                             cores = NULL,
                             ...) {
  # Input validation
  if (!is.data.frame(x)) {
    stop("Input must be a data frame")
  }
  if (nrow(x) == 0) {
    stop("Input data frame is empty")
  }

  required_cols <- c("filename", "start_time", "end_time")
  if (!all(required_cols %in% names(x))) {
    stop(paste("Missing required columns:",
               paste(setdiff(required_cols, names(x)), collapse = ", ")))
  }

  # Handle audio directory
  if (is.null(wav_dir) && is.null(attr(x, "wav_dir"))) {
    stop("wav_dir must be provided either as argument or attribute")
  }
  wav_dir <- wav_dir %||% attr(x, "wav_dir")

  if (nrow(x) == 1) {
    # Single row case
    return(pitch_goodness(segment_row = x,
                          wav_dir = wav_dir,
                          wl = wl,
                          ovlp = ovlp,
                          fmax = fmax,
                          plot = plot))
  } else {
    # Multiple row case - parallel processing
    goodness_list <- parallel_apply(
      indices = 1:nrow(x),
      FUN = function(i) {
        pg_result <- pitch_goodness(
          segment_row = x[i,],
          wav_dir = wav_dir,
          wl = wl,
          ovlp = ovlp,
          fmax = fmax,
          plot = FALSE
        )
        list(time = pg_result[,1], goodness = pg_result[,2])
      },
      cores = cores
    )

    # Extract time and goodness components
    time_list <- lapply(goodness_list, function(x) x$time)
    goodness_list <- lapply(goodness_list, function(x) x$goodness)

    # Create reference time axis from first rendition
    ref_time <- time_list[[1]]
    ref_duration <- round(max(ref_time) - min(ref_time), 2)

    # Normalize all goodness vectors to reference time scale
    goodness_matrix <- sapply(goodness_list, function(goodness) {
      if(length(goodness) != length(ref_time)) {
        approx(
          x = seq(0, 1, length.out = length(goodness)),
          y = goodness,
          xout = seq(0, 1, length.out = length(ref_time))
        )$y
      } else {
        goodness
      }
    })

    # Create result structure
    result <- structure(
      list(
        goodness_matrix = goodness_matrix,
        reference_time = ref_time,
        original_times = time_list
      ),
      class = c("goodness_matrix", "list")
    )

    if (plot) {
      # Validate plot parameters
      if (!is.null(plot_lim)) {
        if (length(plot_lim) != 2 || !is.numeric(plot_lim)) {
          stop("plot_lim must be a numeric vector of length 2")
        }
        goodness_matrix[goodness_matrix < plot_lim[1]] <- plot_lim[1]
        goodness_matrix[goodness_matrix > plot_lim[2]] <- plot_lim[2]
        goodness_range <- plot_lim
      } else {
        goodness_range <- round(range(goodness_matrix, na.rm = TRUE), 2)
      }

      # Configure color mapping
      if (is.null(color_palette)) {
        color_palette <- colorRampPalette(c("black", "red" , "yellow", "white"))
      }

      # Create color breaks
      color_breaks <- seq(goodness_range[1], goodness_range[2],
                          length.out = n_colors)

      # Configure time axis
      time_axis <- list(
        at = seq(1, length(ref_time), length.out = 5),
        labels = sprintf("%.1f", seq(0, ref_duration, length.out = 5))
      )

      # Generate heatmap plot
      heatmap <- lattice::levelplot(
        goodness_matrix,
        col.regions = color_palette(n_colors),
        at = color_breaks,
        xlab = "Time (s)",
        ylab = "Rendition",
        main = paste("Pitch Goodness Heatmap\n",
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
          side = 4,
          width = 1.5,
          height = 0.75,
          padding.text = 2,
          labels = list(
            at = seq(goodness_range[1], goodness_range[2], length.out = 3),
            labels = round(seq(goodness_range[1], goodness_range[2],
                               length.out = 3), 2),
            cex = 0.8
          ),
          title = "Goodness",
          title.gpar = list(
            fontsize = 8,
            fontface = "plain",
            lineheight = 0.8
          )
        ),
        ...
      )

      print(heatmap)
      result$plot <- heatmap
    }

    return(result)
  }
}

#' @rdname goodness
#' @export
goodness.Sap <- function(x,
                         segment_type = c("motifs", "syllables", "segments"),
                         sample_percent = NULL,
                         balanced = FALSE,
                         labels = NULL,
                         clusters = NULL,
                         cores = NULL,
                         seed = 222,
                         wl = 512,
                         ovlp = 50,
                         fmax = 1500,
                         plot = TRUE,
                         plot_lim = NULL,
                         color_palette = NULL,
                         n_colors = 500,
                         ordered = FALSE,
                         descending = TRUE,
                         verbose = TRUE,
                         ...) {

  if(verbose) message(sprintf("\n=== Starting Pitch Goodness Analysis ===\n"))

  # Input validation
  segment_type <- match.arg(segment_type)

  # Special handling for motifs
  if (segment_type == "motifs") {
    # Use original motifs if no ordering/cluster filtering needed
    if (!ordered && is.null(clusters)) {
      segments_df <- x[["motifs"]]
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
    stop("Currently only 'motifs' segment type is supported for pitch goodness analysis")
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

  # Process segments
  goodness_list <- parallel_apply(
    indices = 1:nrow(segments_df),
    FUN = function(i) {
      pg_result <- pitch_goodness(
        segment_row = segments_df[i,],
        wav_dir = x$base_path,
        wl = wl,
        ovlp = ovlp,
        fmax = fmax,
        plot = FALSE
      )
      list(time = pg_result[,1], goodness = pg_result[,2])
    },
    cores = cores
  )

  # Extract components and create matrix
  time_list <- lapply(goodness_list, function(x) x$time)
  goodness_list <- lapply(goodness_list, function(x) x$goodness)

  ref_time <- time_list[[1]]
  time_window <- round(max(ref_time) - min(ref_time), 2)

  goodness_matrix <- sapply(goodness_list, function(goodness) {
    if(length(goodness) != length(ref_time)) {
      approx(
        x = seq(0, 1, length.out = length(goodness)),
        y = goodness,
        xout = seq(0, 1, length.out = length(ref_time))
      )$y
    } else {
      goodness
    }
  })

  # Store attributes
  attr(goodness_matrix, "wl") <- wl
  attr(goodness_matrix, "ovlp") <- ovlp
  attr(goodness_matrix, "segments_df") <- segments_df
  attr(goodness_matrix, "segment_type") <- segment_type
  attr(goodness_matrix, "time_window") <- time_window

  # Set column names as labels
  colnames(goodness_matrix) <- segments_df$label

  # Update SAP object's features
  x$features$motif[["pitch_goodness"]] <- goodness_matrix

  if(verbose) {
    message(sprintf("\nAccess pitch goodness matrix via: x$features$motif$pitch_goodness"))
    message(sprintf("Access attributes via: attributes(x$features$motif$pitch_goodness)"))
  }

  if(plot) {
    # Create heatmap using matrix method
    goodness(goodness_matrix,
             labels = labels,
             plot_lim = plot_lim,
             color_palette = color_palette,
             n_colors = n_colors,
             main = "Heatmap of Motifs: Pitch Goodness")
  }

  invisible(x)
}

#' @rdname goodness
#' @export
goodness.matrix <- function(x,
                            labels = NULL,
                            plot_lim = NULL,
                            color_palette = NULL,
                            n_colors = 500,
                            main = "Pitch Goodness Heatmap",
                            ...) {
  # Check that input is a matrix
  if (!is.matrix(x)) {
    stop("Input must be a matrix")
  }

  # Validate plot_lim
  if (!is.null(plot_lim) &&
      (!is.numeric(plot_lim) || length(plot_lim) != 2)) {
    stop("plot_lim must be a numeric vector of length 2")
  }

  # Require a 'time_window' attribute for the x-axis scaling
  if (is.null(attr(x, "time_window"))) {
    stop("Matrix must have a 'time_window' attribute")
  }
  time_window <- attr(x, "time_window")
  if (!is.numeric(time_window) || length(time_window) != 1) {
    stop("The 'time_window' attribute must be a single numeric value")
  }

  # Check that the matrix has column names
  if (is.null(colnames(x))) {
    stop("Matrix must have column names (to be used as labels)")
  }

  # Handle label subsetting
  if (!is.null(labels)) {
    if (!all(labels %in% colnames(x))) {
      stop("Some provided labels are not found in the matrix's column names")
    }
    x <- x[, colnames(x) %in% labels, drop = FALSE]
    attr(x, "time_window") <- time_window
  }

  # Handle plot limits
  if (!is.null(plot_lim)) {
    x[x < plot_lim[1]] <- plot_lim[1]
    x[x > plot_lim[2]] <- plot_lim[2]
    goodness_range <- plot_lim
  } else {
    goodness_range <- round(range(x, na.rm = TRUE), 2)
  }

  # Reverse matrix for top-to-bottom display
  reversed_matrix <- x[, ncol(x):1, drop = FALSE]

  # Prepare labels and positions
  reversed_labels <- rev(unique(colnames(x)))
  samples_per_label <- rev(as.numeric(table(colnames(x))))
  cumulative_positions <- cumsum(c(0, head(samples_per_label, -1)))
  label_positions <- cumulative_positions + samples_per_label / 2
  hline_positions <- cumsum(samples_per_label)[-length(samples_per_label)]

  # Set default color palette
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("black", "red", "yellow", "white"))
  }

  # Create color breaks
  color_breaks <- seq(goodness_range[1], goodness_range[2],
                      length.out = n_colors)

  # Configure time axis
  x_scale <- list(
    at = seq(0, nrow(reversed_matrix), length.out = 5),
    labels = sprintf("%.1f", seq(0, time_window, length.out = 5))
  )

  # Create heatmap
  heatmap <- lattice::levelplot(
    reversed_matrix,
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
        at = seq(goodness_range[1], goodness_range[2], length.out = 3),
        labels = sprintf("%.2f", seq(goodness_range[1], goodness_range[2],
                                     length.out = 3)),
        cex = 1,
        col = "black"
      ),
      title = "Goodness",
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



#' Internal function for processing single-row fundamental frequency estimation
#'
#' @keywords internal
pitch_goodness <- function(segment_row,
                           wav_dir = NULL,
                           wl = 512,
                           ovlp = 80,
                           fmax = 1500,
                           plot = TRUE,
                           ...) {

  # Check input validity
  if (!is.data.frame(segment_row)) stop("segment_row must be a data frame")
  if (!all(c("start_time", "end_time") %in% names(segment_row))) {
    stop("Input must contain 'start_time' and 'end_time' columns")
  }

  # Construct file path
  sound_path <- construct_wav_path(segment_row, wav_dir = wav_dir)
  if (!file.exists(sound_path)) stop("File not found: ", sound_path)

  # Read audio segment with timing
  wv <- tuneR::readWave(sound_path,
                        from = segment_row$start_time,
                        to = segment_row$end_time,
                        units = "seconds")

  # Get waveform parameters
  waveform <- wv@left
  sr <- wv@samp.rate

  # Calculate analysis timestamps
  duration <- segment_row$end_time - segment_row$start_time
  hop <- wl * (1 - ovlp/100)
  n_windows <- floor((length(waveform) - wl)/hop) + 1
  times <- seq(segment_row$start_time + wl/(2*sr),  # Center of window
               by = hop/sr,
               length.out = n_windows)

  # Compute cepstral goodness
  # Compute windowed FFT
  spec <- seewave::spectro(waveform, f = sr, wl = wl, ovlp = ovlp,
                           plot = FALSE, complex = TRUE, norm = FALSE, dB = NULL)

  # Compute cepstra for each window
  cepstra <- apply(spec$amp, 2, function(window_fft) {
    log_mag <- log(Mod(window_fft)^2 + .Machine$double.eps)
    cepstrum <- Re(stats::fft(log_mag, inverse = TRUE)/length(log_mag))
    cepstrum[1:ceiling(length(cepstrum)/2)]
  })

  # Get quefrency values
  quefrencies <- seq(0, nrow(cepstra)-1)/sr

  # Calculate cutoff index for fmax
  cutoff_idx <- which.min(abs(quefrencies - 1/fmax))
  upper_idx <- floor(nrow(cepstra)/2)

  # Compute goodness of pitch
  goodness <- apply(cepstra, 2, function(cepstrum_window) {
    valid_range <- cepstrum_window[cutoff_idx:upper_idx]
    if(all(is.na(valid_range))) NA else max(valid_range, na.rm = TRUE)
  })

  # Create result matrix
  res <- cbind(time = times, goodness = goodness)
  colnames(res) <- c("time", "goodness")

  if(plot) {
    # Generate the spectrogram using the visualize_song.default function
    visualize_song.default(x = sound_path,
                           start_time_in_second = segment_row$start_time,
                           end_time_in_second = segment_row$end_time,
                           fft_window_size = wl,
                           overlap = ovlp/100,
                           dark_mode = TRUE,
                           legend = FALSE,
                           keep.par = TRUE,
                           verbose = FALSE)

    offset <- segment_row$duration * 0.035

    # Overlay the pitch goodness trace
    par(new = TRUE)
    plot(res[,"time"], res[,"goodness"],
         type = "l",
         col = "white",
         lwd = 2,
         axes = FALSE,
         xlab = "",
         ylab = "",
         xlim = c(segment_row$start_time + offset , segment_row$end_time - offset))

    # Add a right-side axis for the pitch goodness values
    axis(side = 4, col = "white", col.axis = "white")
    mtext("Pitch Goodness", side = 4, line = 3, col = "white")

    # Add a title to the plot
    title(main = "Spectrogram w/ Goodness of Pitch Overlay")
  }


  return(res)
}



# Refine FF ---------------------------------------------------------------

#' Refine Fundamental Frequency Detection
#'
#' @description
#' Refines fundamental frequency detection by identifying regions of reliable pitch tracking
#' using various methods and applying temporal constraints.
#'
#' @param x A Sap object containing fundamental frequency and pitch goodness data
#' @param segment_type Character, type of segment to analyze: "motifs", "syllables", or "segments"
#' @param reference_label Character, label to use as reference for template creation
#' @param method Character, method for template creation: "quantile" or "hmm"
#' @param minimal_duration Numeric, minimum duration (in ms) for valid segments
#' @param split_dips Logical, whether to split segments at significant dips in goodness
#' @param quantile_threshold Numeric, threshold for quantile method (0-1)
#' @param random_seed Numeric, seed for reproducibility
#' @param plot Logical, whether to plot the results
#' @param plot_freq_lim Numeric vector of length 2, frequency limits for plotting
#' @param color_palette Function, color palette for plotting
#' @param stats Logical, whether to calculate segment statistics
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns the modified Sap object with refined fundamental frequency data
#'
#' @details
#' The function uses either quantile-based thresholding or Hidden Markov Models to identify
#' regions of reliable pitch tracking. It can optionally split segments at local minima in
#' pitch goodness and applies minimum duration constraints.
#'
#' @examples
#' \dontrun{
#' sap3 <- refine_FF(sap1,
#'                   reference_label = "BL",
#'                   method = "quantile",
#'                   plot = TRUE)
#' }
#'
#' @export
refine_FF <- function(x,
                     segment_type = c("motifs", "syllables", "segments"),
                     reference_label,
                     method = c("quantile", "hmm"),
                     minimal_duration = 20, # in ms
                     split_dips = TRUE,
                     quantile_threshold = 0.5,
                     random_seed = 222,
                     plot = TRUE,
                     plot_freq_lim = NULL,
                     color_palette = NULL,
                     stats = FALSE,
                     stats_across_label = FALSE,
                      ...) {

  # Input validation
  segment_type <- match.arg(segment_type)
  method <- match.arg(method)

  # Check if input is a Sap object
  if (!inherits(x, "Sap")) {
    stop("Input must be a Sap object")
  }

  # Update feature type by removing 's' from segment_type
  feature_type <- sub("s$", "", segment_type)

  # Get the relevant matrices
  ff_matrix <- x$features[[feature_type]]$fund_freq
  goodness_matrix <- x$features[[feature_type]]$pitch_goodness

  # Validate matrices
  if (is.null(ff_matrix) || is.null(goodness_matrix)) {
    stop("Both fundamental frequency and pitch goodness matrices must exist")
  }
  if (!identical(dim(ff_matrix), dim(goodness_matrix))) {
    stop("Matrices must have identical dimensions")
  }
  if (!reference_label %in% colnames(ff_matrix)) {
    stop("Reference label not found in data")
  }

  # Calculate time parameters
  time_window <- attr(ff_matrix, "time_window")
  time_step <- time_window / nrow(ff_matrix)
  min_samples <- ceiling(minimal_duration / (time_step * 1000))

  # Get reference label columns and average them
  ref_cols <- which(colnames(ff_matrix) == reference_label)
  ref_goodness <- rowMeans(goodness_matrix[, ref_cols, drop = FALSE])

  # Core segmentation logic
  template <- switch(method,
           "quantile" = {
             quantile_val <- quantile(ref_goodness, quantile_threshold, na.rm = TRUE)
             ref_goodness > quantile_val
           },

           "hmm" = {
             if (!requireNamespace("depmixS4", quietly = TRUE)) {
               stop("Install depmixS4: install.packages('depmixS4')")
             }

             if (!is.null(random_seed)) set.seed(random_seed)

             hmm_data <- data.frame(goodness = ref_goodness)
             hmm_model <- depmixS4::depmix(
               goodness ~ 1,
               data = hmm_data,
               nstates = 2,
               family = gaussian(),
               instart = c(0.5, 0.5),
               trstart = matrix(c(hmm_trans_prob, 1-hmm_trans_prob,
                                  1-hmm_trans_prob, hmm_trans_prob), ncol = 2),
               ntimes = nrow(hmm_data)
             )

             hmm_fit <- tryCatch({
               fit <- depmixS4::fit(hmm_model, verbose = FALSE)
               if(fit@message != "converged") {
                 NULL  # Return NULL if not converged
               } else {
                 fit
               }
             }, error = function(e) NULL)

             if(!is.null(hmm_fit)) {
               states <- depmixS4::posterior(hmm_fit, type = "viterbi")$state
               state_means <- tapply(ref_goodness, states, mean)
               states == which.max(state_means)
             } else {
               warning("HMM failed to converge, using quantile fallback (threshold = 0.5)")
               quantile_val <- quantile(ref_goodness, 0.5, na.rm = TRUE)
               ref_goodness > quantile_val
             }
           }

  )

  # Apply minimal duration filter
  template <- filter_short_segments(template, min_samples)

  # Segmentation splitting at local minima
  if(split_dips) {
    template <- split_at_dips(template, ref_goodness, min_samples)
  }

  # Create filtered FF matrix
  filtered_ff <- ff_matrix
  filtered_ff[!template, ] <- NA
  attributes(filtered_ff) <- attributes(ff_matrix)
  attr(filtered_ff, "temporal_template") <- template

  # Update the Sap object
  x$features[[feature_type]]$filtered_fund_freq <- filtered_ff

  if (stats) {
  # Calculate segment statistics
  segment_list <- calculate_segment_stats(feature_matrix = filtered_ff,
                                           template = template,
                                           time_step = time_step,
                                           stats_across_label = stats_across_label)


  if (stats_across_label) {
    # When comparing across labels
    x$features[[feature_type]]$stats_fund_freq <- segment_list$segment_stats
    x$features[[feature_type]]$fund_freq_across_label <- segment_list$label_comparisons
  } else {
    # When not comparing across labels
    x$features[[feature_type]]$stats_fund_freq <- segment_list
  }
  }

  # Plot
  if (plot) {
    heatmap <- FF.matrix(
      ff_matrix,
      plot_freq_lim = plot_freq_lim,
      color_palette = color_palette,
      main = paste("Fundamental Frequency with", method, "Template")
    )

    # Overlay template
    lattice::trellis.focus("panel", 1, 1, highlight = FALSE)
    template_pos <- which(template)
    if(length(template_pos) > 0) {
      lattice::panel.rect(
        xleft = template_pos - 0.5,
        xright = template_pos + 0.5,
        ybottom = 0.5,
        ytop = ncol(filtered_ff) + 0.5,
        border = NA,
        col = rgb(0, 1, 0, 0.5)  # Cyan overlay
      )
    }
    lattice::trellis.unfocus()
  }

  invisible(x)
}


#' Filter Short Segments from Template
#'
#' @param template Logical vector representing the template
#' @param min_samples Minimum number of samples for a valid segment
#' @return Logical vector with short segments removed
#' @keywords internal
filter_short_segments <- function(template, min_samples) {
  rle_template <- rle(template)
  for(i in seq_along(rle_template$lengths)) {
    if(rle_template$values[i] && rle_template$lengths[i] < min_samples) {
      start <- sum(rle_template$lengths[1:i-1]) + 1
      end <- start + rle_template$lengths[i] - 1
      template[start:end] <- FALSE
    }
  }
  template
}

#' Split Template at Local Minima in Goodness
#'
#' @param template Logical vector representing the template
#' @param goodness Numeric vector of goodness values
#' @param min_duration_samples Minimum number of samples for a valid segment
#' @return Logical vector with segments split at significant dips
#' @keywords internal
split_at_dips <- function(template, goodness, min_duration_samples) {
  runs <- rle(template)
  new_template <- template

  for(i in which(runs$values)) {
    start <- sum(runs$lengths[1:(i-1)]) + 1
    end <- start + runs$lengths[i] - 1
    segment_goodness <- goodness[start:end]

    # Find local minima
    minima <- which(diff(sign(diff(segment_goodness))) == 2) + 1
    if(length(minima) > 0) {
      dip_threshold <- quantile(segment_goodness, 0.25)
      valid_minima <- minima[segment_goodness[minima] < dip_threshold]

       # Apply splits
       for(vm in valid_minima) {
         split_pos <- start + vm - 1
         new_template[split_pos] <- FALSE
       }
    }
  }

  # Return modified template and re-apply duration filter
  filter_short_segments(new_template, min_duration_samples)
}

#' @keywords internal
calculate_segment_stats <- function(feature_matrix,
                                    template,
                                    time_step,
                                    feature = "pitch",
                                    stats_across_label = FALSE,
                                    ...) {

  # Validate inputs
  if (!is.matrix(feature_matrix)) stop("feature_matrix must be a matrix")
  if (!is.logical(template)) stop("template must be a logical vector")
  if (nrow(feature_matrix) != length(template)) {
    stop("feature_matrix rows must match template length")
  }

  # Get labels from column names
  labels <- colnames(feature_matrix)

  # Identify contiguous segments
  rle_template <- rle(template)
  segment_ends <- cumsum(rle_template$lengths)
  segment_starts <- c(1, segment_ends[-length(segment_ends)] + 1)
  valid_segments <- which(rle_template$values)

  # Initialize list to store results per label
  stats_list <- list()

  # Process each label separately
  for (label in labels) {
    # Get columns belonging to this label
    label_cols <- which(colnames(feature_matrix) == label)
    if (length(label_cols) == 0) next

    # Calculate per-rendition statistics within each segment
    rendition_stats <- lapply(label_cols, function(col) {
      seg_stats <- lapply(valid_segments, function(i) {
        seg_range <- segment_starts[i]:segment_ends[i]
        seg_data <- feature_matrix[seg_range, col]
        valid_vals <- seg_data[!is.na(seg_data)]

        if (length(valid_vals) > 0) {
          data.frame(
            segment_id = i,
            start_time = (segment_starts[i] - 1) * time_step,
            end_time = segment_ends[i] * time_step,
            duration = (segment_ends[i] - segment_starts[i] + 1) * time_step,
            mean_val = mean(valid_vals),
            #sd_val = sd(valid_vals),
            min_val = min(valid_vals),
            max_val = max(valid_vals),
            n_samples = length(valid_vals),
            stringsAsFactors = FALSE
          )
        } else NULL
      })
      do.call(rbind, seg_stats)
    })

    # Combine all renditions for this label
    label_df <- do.call(rbind, rendition_stats)

    if (!is.null(label_df) && nrow(label_df) > 0) {
      # Calculate label-level aggregates across renditions
      label_stats <- label_df %>%
        dplyr::group_by(segment_id, start_time, end_time, duration) %>%
        dplyr::summarize(
          rendition_count = dplyr::n(),
          mean_freq = mean(mean_val, na.rm = TRUE),
          sd_freq = sd(mean_val, na.rm = TRUE),
          min_freq = min(min_val, na.rm = TRUE),
          max_freq = max(max_val, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(label = label)

      stats_list[[label]] <- label_stats
    }
  }

  # Combine all labels into final data frame
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL

  # Reorder columns
  stats_df <- stats_df %>%
    dplyr::select(label, segment_id, start_time, end_time, duration,
           rendition_count, everything())

  # Add statistical comparisons across labels if requested
  if (stats_across_label && length(unique(stats_df$label)) > 1) {
    comparison_stats <- compare_labels(label_df)
    stats_df <- list(
      segment_stats = stats_df,
      label_comparisons = comparison_stats
    )
  }

  return(stats_df)
}


# Helper function for label comparisons
compare_labels <- function(label_df) {
  require(dplyr)
  require(broom)

  # Initialize results list
  comparison_results <- list()

  # 1. ANOVA for each segment across labels
  anova_results <- label_df %>%
    group_by(segment_id) %>%
    filter(n_distinct(label) > 1) %>%  # Only compare segments with multiple labels
    do(tidy(aov(mean_freq ~ label, data = .))) %>%
    ungroup()
    #  %>%
    # mutate(significant = p.value < 0.05)

  comparison_results$anova <- anova_results

  # 2. Pairwise t-tests for each segment
  pairwise_results <- label_df %>%
    group_by(segment_id) %>%
    filter(n_distinct(label) > 1) %>%
    do(tidy(pairwise.t.test(.$mean_freq, .$label, p.adjust.method = "holm"))) %>%
    ungroup()

  comparison_results$pairwise <- pairwise_results

  # 3. Overall label differences (across all segments)
  overall_test <- tryCatch({
    tidy(kruskal.test(mean_freq ~ label, data = label_df))
  }, error = function(e) NULL)

  comparison_results$overall <- overall_test

  # 4. Effect size (Cohen's d for pairwise comparisons)
  effect_sizes <- label_df %>%
    group_by(segment_id) %>%
    filter(n_distinct(label) > 1) %>%
    do({
      labels <- unique(.$label)
      combn(labels, 2, simplify = FALSE) %>%
        purrr::map_dfr(function(pair) {
          x <- filter(., label == pair[1])$mean_freq
          y <- filter(., label == pair[2])$mean_freq
          data.frame(
            segment_id = first(.$segment_id),
            comparison = paste(pair, collapse = " vs "),
            cohens_d = abs(mean(x) - mean(y))/sqrt((sd(x)^2 + sd(y)^2)/2)
          )
        })
    }) %>%
    ungroup()

  comparison_results$effect_sizes <- effect_sizes

  return(comparison_results)
}
