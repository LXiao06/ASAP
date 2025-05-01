# Entropy Module

# Spectral Entropy Estimation ----------------------------------------------
# Update date : Apr.7, 2025

#' Calculate Spectral Entropy for Audio Segments
#'
#' @description
#' Calculates spectral entropy (Wiener or Shannon) for audio segments, providing measures
#' of spectral uniformity or complexity in sound signals.
#'
#' @param x Input object:
#'   \itemize{
#'     \item character: path to WAV file (default method)
#'     \item data frame: containing segment information
#'     \item SAP object
#'     \item Pre-computed matrix for spectral entropy
#'   }
#' @param start_time Numeric, start time in seconds (for default method)
#' @param end_time Numeric, end time in seconds (for default method)
#' @param wav_dir Directory containing WAV files (for data frame methods)
#' @param wl Window length for FFT (default: 512)
#' @param wn Window name for FFT (default: "hanning")
#' @param ovlp Overlap percentage between windows (default: 50)
#' @param fftw Logical, use FFTW or not (default: TRUE)
#' @param freq_range Frequency range for analysis c(min, max) in Hz (default: c(500, 15000))
#' @param threshold Amplitude threshold for power spectrum (default: 10)
#' @param method Entropy type ("weiner" or "shannon")
#' @param normalize Logical, whether to normalize entropy values (default: FALSE)
#' @param plot Logical, whether to plot results (default: TRUE)
#' @param plot_entropy_lim Optional limits for entropy plot c(min, max)
#' @param color_palette Custom color palette function
#' @param n_colors Number of colors in palette (default: 500)
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
#' @param ... Additional arguments passed to levelplot
#'
#' @details
#' Wiener entropy measures the uniformity of power distribution:
#' - Non-normalized (default): returns log-scaled values from -Inf to 0
#'   - 0 indicates uniform distribution (white noise)
#'   - Large negative values indicate structured sound (pure tones)
#' - Normalized: returns values from 0 to 1
#'   - 1 indicates uniform distribution
#'   - Values close to 0 indicate structured sound
#'
#' Shannon entropy measures information content:
#' - Non-normalized: returns values >= 0
#'   - Higher values indicate more uniformity/randomness
#'   - Lower values indicate more structure/predictability
#' - Normalized: returns values from 0 to 1
#'   - Scale relative to maximum possible entropy for vector length
#'
#' @return
#' Depending on the method used:
#' \itemize{
#'   \item Default method (WAV file):
#'     \describe{
#'       \item{Matrix}{With columns: time, entropy}
#'     }
#'   \item Data frame method:
#'     \describe{
#'       \item{List}{With components: entropy_matrix, reference_time, original_times, plot (if plot=TRUE)}
#'     }
#'   \item SAP method:
#'     \describe{
#'       \item{SAP object}{Updated with entropy features}
#'     }
#'   \item Matrix method:
#'     \describe{
#'       \item{Lattice plot}{(invisibly)}
#'     }
#' }
#'
#' @examples
#' \dontrun{
#' # Calculate entropy for a single WAV file
#' entropy <- spectral_entropy("path/to/sound.wav",
#'                           start_time = 1,
#'                           end_time = 2,
#'                           method = "weiner",
#'                           normalize = FALSE)
#'
#' # Calculate normalized Shannon entropy for multiple segments using data frame method
#' entropy_list <- spectral_entropy(segments,
#'                                wav_dir = "path/to/wavs",
#'                                method = "shannon",
#'                                normalize = TRUE,
#'                                freq_range = c(1000, 10000))
#'
#' # Method for Sap objects
#' sap <- spectral_entropy(sap,
#'                        segment_type = "motifs",
#'                        method = "weiner",
#'                        normalize = FALSE,
#'                        plot = TRUE)
#'
#' # Access the entropy matrix
#' wiener_entropy <- sap$features$motif$weiner_entropy
#'
#' # Method for entropy matrices
#' # Plot existing entropy matrix with custom settings
#' spectral_entropy(wiener_entropy,
#'                 plot_entropy_lim = c(-3, 0),
#'                 color_palette = colorRampPalette(c("purple", "blue",
#'                                                   "cyan", "yellow",
#'                                                   "orange", "red")),
#'                 main = "Custom Entropy Visualization")
#'
#' # Plot subset of entropy matrix using labels
#' spectral_entropy(wiener_entropy,
#'                 labels = c("a", "b", "c"),
#'                 plot_entropy_lim = c(0, 1),
#'                 n_colors = 1000)
#' }
#'
#' @export
spectral_entropy <- function(x, ...) {
  UseMethod("spectral_entropy")
}

#' @rdname spectral_entropy
#' @export
spectral_entropy.default <- function(x,  # x is wav file path
                                     start_time = NULL,
                                     end_time = NULL,
                                     wl = 512,
                                     wn = "hanning",
                                     ovlp = 50,
                                     fftw = TRUE,
                                     freq_range = c(500, 15000),
                                     threshold = 10,
                                     method = c("weiner", "shannon"),
                                     normalize = FALSE,
                                     plot = TRUE,
                                     ...) {
  # Match method argument
  method <- match.arg(method)

  # Validate file path
  if (!file.exists(x)) {
    stop("File does not exist: ", x)
  }

  # Always read header first to get duration
  header <- tuneR::readWave(x, header = TRUE)
  duration <- header$samples / header$sample.rate

  # Handle start_time
  if (is.null(start_time)) {
    start_time <- 0
  }

  # Handle end_time
  if (is.null(end_time)) {
    end_time <- duration
  } else {
    # Validate user-provided end_time against actual duration
    if (end_time > duration) {
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

  # Create a single-row data frame for sh_single_row
  segment_row <- data.frame(
    filename = basename(x),
    start_time = start_time,
    end_time = end_time,
    stringsAsFactors = FALSE
  )

  # Call sh_single_row with the constructed data
  if (fftw) ensure_pkgs("fftw")

  sh_result <- sh_single_row(
    segment_row = segment_row,
    wav_dir = dirname(x),
    wl = wl,
    wn = wn,
    ovlp = ovlp,
    fftw = fftw,
    freq_range = freq_range,
    threshold = threshold,
    method = method,
    normalize = normalize,
    plot = plot,
    ...
  )

  return(sh_result)
}



#' @rdname spectral_entropy
#' @export
spectral_entropy.data.frame <- function(x,
                                     wav_dir = NULL,
                                     wl = 512,
                                     wn = "hanning",
                                     ovlp = 50,
                                     fftw = TRUE,
                                     freq_range = c(500, 15000),
                                     threshold = 10,
                                     method = c("weiner", "shannon"),
                                     normalize = FALSE,
                                     plot = TRUE,
                                     plot_entropy_lim = NULL,
                                     color_palette = NULL,
                                     n_colors = 500,
                                     cores = NULL,
                                     ...) {
  # Match method argument
  method <- match.arg(method)

  # Chcek fftw
  if (fftw) ensure_pkgs("fftw")

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
    return(sh_single_row(segment_row = x,
                         wav_dir = wav_dir,
                         wl = wl,
                         wn = wn,
                         ovlp = ovlp,
                         fftw = fftw,
                         freq_range = freq_range,
                         threshold = threshold,
                         method = method,
                         normalize = normalize,
                         plot = plot))
  } else {
    # Parallel processing of entropy extraction
    sh_list <- parallel_apply(
      indices = 1:nrow(x),
      FUN = function(i) {
        sh_result <- sh_single_row(
          segment_row = x[i,],
          wav_dir = wav_dir,
          wl = wl,
          wn = wn,
          ovlp = ovlp,
          fftw = fftw,
          freq_range = freq_range,
          threshold = threshold,
          method = method,
          normalize = normalize,
          plot = FALSE
        )
        list(time = sh_result[,1], entropy = sh_result[,2])
      },
      cores = cores
    )

    # Extract time and entropy components
    time_list <- lapply(sh_list, function(x) x$time)
    entropy_list <- lapply(sh_list, function(x) x$entropy)

    # Create reference time axis from first rendition
    ref_time <- time_list[[1]]
    ref_duration <- round(max(ref_time) - min(ref_time), 2)

    # Normalize all entropy vectors to reference time scale
    entropy_matrix <- sapply(entropy_list, function(entropy) {
      if(length(entropy) != length(ref_time)) {
        approx(
          x = seq(0, 1, length.out = length(entropy)),
          y = entropy,
          xout = seq(0, 1, length.out = length(ref_time))
        )$y
      } else {
        entropy
      }
    })

    # Return structure includes original times for reference
    result <- structure(
      list(
        entropy_matrix = entropy_matrix,
        reference_time = ref_time,
        original_times = time_list
      ),
      class = c("entropy_matrix", "list")
    )

    if (plot) {
      # Validate plot parameters
      if (!is.null(plot_entropy_lim)) {
        if (length(plot_entropy_lim) != 2 || !is.numeric(plot_entropy_lim)) {
          stop("plot_entropy_lim must be a numeric vector of length 2")
        }
      }

      # Handle NA/Inf values in entropy matrix
      entropy_matrix[!is.finite(entropy_matrix)] <- NA

      # Modify entropy matrix capping
      if (!is.null(plot_entropy_lim)) {
        entropy_matrix[entropy_matrix < plot_entropy_lim[1]] <- plot_entropy_lim[1]
        entropy_matrix[entropy_matrix > plot_entropy_lim[2]] <- plot_entropy_lim[2]

        # Update range for color breaks
        entropy_range <- plot_entropy_lim
      } else {
        # Only calculate entropy_range if plot_entropy_lim is not provided
        entropy_range <- round(range(entropy_matrix, na.rm = TRUE), 2)
      }

      # Configure color mapping
      if (is.null(color_palette)) {
        if(normalize) {
          # For normalized entropy (0-1)
          main_colors <- colorRampPalette(c("darkblue", "blue", "cyan",
                                            "yellow", "orange", "red"))
        } else {
          # For non-normalized entropy
          main_colors <- colorRampPalette(c("purple", "blue", "cyan",
                                            "yellow", "orange", "red"))
        }
      }

      # Modify the matrix to use a specific value for NA
      entropy_matrix_mod <- entropy_matrix
      na_value <- entropy_range[1] - (entropy_range[2] - entropy_range[1])/10  # Smaller gap for NA
      entropy_matrix_mod[is.na(entropy_matrix_mod)] <- na_value

      # Create color breaks with smaller NA region
      # Allocate ~10% of the breaks to NA
      na_breaks <- seq(na_value, entropy_range[1], length.out = round(n_colors * 0.1))
      value_breaks <- seq(entropy_range[1], entropy_range[2], length.out = n_colors - length(na_breaks) + 1)
      color_breaks <- c(na_breaks[-length(na_breaks)], value_breaks)

      # Generate corresponding colors
      na_colors <- rep("grey20", length(na_breaks) - 1)
      value_colors <- main_colors(length(value_breaks))

      # Configure visualization parameters
      time_axis <- list(
        at = seq(1, length(ref_time), length.out = 5),
        labels = sprintf("%.1f", seq(0, ref_duration, length.out = 5))
      )

      # Generate heatmap plot
      heatmap <- lattice::levelplot(
        entropy_matrix_mod,
        col.regions =  c(na_colors, value_colors),
        at = color_breaks,
        xlab = "Time (s)",
        ylab = "Rendition",
        main = paste("Spectral Entropy Heatmap\n",
                     tools::toTitleCase(method), "Entropy",
                     ifelse(normalize, "(Normalized)", "(Non-normalized)"),
                     "\nWindow:", wl, "samples | Overlap:", ovlp, "%"),
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
            at = c(na_value + (entropy_range[1] - na_value)/2,  # Center of NA region
                   entropy_range[1],
                   mean(entropy_range),
                   entropy_range[2]),
            labels = c("NA",
                       round(c(entropy_range[1],
                               mean(entropy_range),
                               entropy_range[2]), 2)),
            cex = 0.8
          ),
          title = "Entropy",
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

#' @rdname spectral_entropy
#' @export
spectral_entropy.Sap <- function(x,
                                 segment_type = c("motifs", "syllables", "segments"),
                                 sample_percent = NULL,
                                 balanced = FALSE,
                                 labels = NULL,
                                 clusters = NULL,
                                 cores = NULL,
                                 seed = 222,
                                 wl = 512,
                                 wn = "hanning",
                                 ovlp = 50,
                                 fftw = TRUE,
                                 freq_range = c(500, 15000),
                                 threshold = 10,
                                 method = c("weiner", "shannon"),
                                 normalize = FALSE,
                                 plot = TRUE,
                                 plot_entropy_lim = NULL,
                                 color_palette = NULL,
                                 n_colors = 500,
                                 ordered = FALSE,
                                 descending = TRUE,
                                 verbose = TRUE,
                                 ...) {
  if(verbose) message(sprintf("\n=== Starting Spectral Entropy Analysis ===\n"))

  # Match arguments
  method <- match.arg(method)
  segment_type <- match.arg(segment_type)

  # Chcek fftw
  if (fftw) ensure_pkgs("fftw")

  # Special handling for motifs
  if (segment_type == "motifs") {
    if (!ordered && is.null(clusters)) {
      segments_df <- x[["motifs"]]
    } else {
      if (is.null(x$features$motif$feat.embeds)) {
        stop("Feature embeddings required for ordered/clustered motif plots")
      }

      segments_df <- x$features$motif$feat.embeds |>
        as_segment()

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
    stop("Currently only 'motifs' segment type is supported for spectral entropy heatmap")
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

  if (nrow(segments_df) == 0) {
    stop("No segments remaining after subsetting. Check labels/clusters.")
  }

  # Process segments
  sh_list <- parallel_apply(
    indices = 1:nrow(segments_df),
    FUN = function(i) {
      sh_result <- sh_single_row(
        segment_row = segments_df[i,],
        wav_dir = x$base_path,
        wl = wl,
        wn = wn,
        ovlp = ovlp,
        fftw = fftw,
        freq_range = freq_range,
        threshold = threshold,
        method = method,
        normalize = normalize,
        plot = FALSE
      )
      list(time = sh_result[,1], entropy = sh_result[,2])
    },
    cores = cores
  )

  # Extract components
  time_list <- lapply(sh_list, function(x) x$time)
  entropy_list <- lapply(sh_list, function(x) x$entropy)

  # Create reference time axis
  ref_time <- time_list[[1]]
  time_window <- round(max(ref_time) - min(ref_time), 2)

  # Normalize all entropy vectors
  entropy_matrix <- sapply(entropy_list, function(entropy) {
    if(length(entropy) != length(ref_time)) {
      approx(
        x = seq(0, 1, length.out = length(entropy)),
        y = entropy,
        xout = seq(0, 1, length.out = length(ref_time))
      )$y
    } else {
      entropy
    }
  })

  # Store attributes
  attr(entropy_matrix, "wl") <- wl
  attr(entropy_matrix, "ovlp") <- ovlp
  attr(entropy_matrix, "segments_df") <- segments_df
  attr(entropy_matrix, "segment_type") <- segment_type
  attr(entropy_matrix, "time_window") <- time_window
  attr(entropy_matrix, "method") <- method
  attr(entropy_matrix, "normalize") <- normalize

  # Set column names
  colnames(entropy_matrix) <- segments_df$label

  # Update SAP object's features
  x$features$motif[[paste0(method, "_entropy")]] <- entropy_matrix

  if(verbose) {
    message(sprintf("\nAccess entropy matrix via: x$features$motif$%s_entropy", method))
    message(sprintf("Access attributes via: attributes(x$features$motif$%s_entropy)", method))
  }

  if(plot) {
    # Create heatmap using matrix method
    spectral_entropy(entropy_matrix,
                     labels = labels,
                     plot_entropy_lim = plot_entropy_lim,
                     color_palette = color_palette,
                     n_colors = n_colors,
                     main = sprintf("Heatmap of Motifs: %s Entropy (%s)",
                                    tools::toTitleCase(method),
                                    ifelse(normalize, "Normalized", "Non-normalized")))
  }

  invisible(x)
}

#' @rdname spectral_entropy
#' @export
spectral_entropy.matrix <- function(x,
                                    labels = NULL,
                                    plot_entropy_lim = NULL,
                                    color_palette = NULL,
                                    n_colors = 500,
                                    main = "Spectral Entropy Heatmap",
                                    ...) {
  # Check input
  if (!is.matrix(x)) {
    stop("Input must be a matrix")
  }

  # Validate plot_entropy_lim
  if (!is.null(plot_entropy_lim) &&
      (!is.numeric(plot_entropy_lim) || length(plot_entropy_lim) != 2)) {
    stop("plot_entropy_lim must be a numeric vector of length 2")
  }

  # Check time_window attribute
  if (is.null(attr(x, "time_window"))) {
    stop("Matrix must have a 'time_window' attribute")
  }
  time_window <- attr(x, "time_window")
  if (!is.numeric(time_window) || length(time_window) != 1) {
    stop("The 'time_window' attribute must be a single numeric value")
  }

  # Check column names
  if (is.null(colnames(x))) {
    stop("Matrix must have column names (to be used as labels)")
  }

  # Handle labels subsetting
  if (!is.null(labels)) {
    if (!all(labels %in% colnames(x))) {
      stop("Some provided labels are not found in the matrix's column names")
    }
    x <- x[, colnames(x) %in% labels, drop = FALSE]
    attr(x, "time_window") <- time_window
  }

  # Handle NA/Inf values
  x[!is.finite(x)] <- NA

  # Calculate entropy range
  if (!is.null(plot_entropy_lim)) {
    entropy_range <- plot_entropy_lim
  } else {
    entropy_range <- round(range(x, na.rm = TRUE), 2)
  }

  # Reverse matrix for plotting
  reversed_entropy_matrix <- x[, ncol(x):1, drop = FALSE]

  # Prepare labels and positions
  reversed_labels <- rev(unique(colnames(x)))
  samples_per_label <- rev(as.numeric(table(colnames(x))))
  cumulative_positions <- cumsum(c(0, head(samples_per_label, -1)))
  label_positions <- cumulative_positions + samples_per_label / 2
  hline_positions <- cumsum(samples_per_label)[-length(samples_per_label)]

  # Configure color mapping
  if (is.null(color_palette)) {
    method <- attr(x, "method") %||% "weiner"
    normalize <- attr(x, "normalize") %||% FALSE

    if(normalize) {
      main_colors <- colorRampPalette(c("darkblue", "blue", "cyan",
                                        "yellow", "orange", "red"))
    } else {
      main_colors <- colorRampPalette(c("purple", "blue", "cyan",
                                        "yellow", "orange", "red"))
    }
  } else {
    main_colors <- color_palette
  }

  # Modify matrix for NA representation
  entropy_matrix_mod <- reversed_entropy_matrix
  na_value <- entropy_range[1] - (entropy_range[2] - entropy_range[1])/10
  entropy_matrix_mod[is.na(entropy_matrix_mod)] <- na_value

  # Create color breaks
  na_breaks <- seq(na_value, entropy_range[1], length.out = round(n_colors * 0.1))
  value_breaks <- seq(entropy_range[1], entropy_range[2], length.out = n_colors - length(na_breaks) + 1)
  color_breaks <- c(na_breaks[-length(na_breaks)], value_breaks)

  # Generate corresponding colors
  na_colors <- rep("grey20", length(na_breaks) - 1)
  value_colors <- main_colors(length(value_breaks))

  # Create the heatmap
  heatmap <- lattice::levelplot(
    entropy_matrix_mod,
    col.regions = c(na_colors, value_colors),
    at = color_breaks,
    border = "transparent",
    xlab = "Time (s)",
    ylab = "Rendition",
    main = main,
    scales = list(
      x = list(
        at = seq(0, nrow(entropy_matrix_mod), length.out = 5),
        labels = sprintf("%.1f", seq(0, time_window, length.out = 5))
      ),
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
      side = 4,
      width = 1.5,
      height = 0.75,
      padding.text = 2,
      labels = list(
        at = c(na_value + (entropy_range[1] - na_value)/2,
               entropy_range[1],
               mean(entropy_range),
               entropy_range[2]),
        labels = c("NA",
                   round(c(entropy_range[1],
                           mean(entropy_range),
                           entropy_range[2]), 2)),
        cex = 0.8
      ),
      title = "Entropy",
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


#' Calculate continuous spectral entropy for a single audio segment
#'
#' @param segment_row A data frame row containing start_time and end_time
#' @param wav_dir Directory containing wav files
#' @param wl Window length for FFT (default: 512)
#' @param wn Window name for FFT (default: "hanning")
#' @param ovlp Overlap percentage between windows (default: 50)
#' @param fftw Logical, use FFTW or not (default: TRUE)
#' @param freq_range Frequency range for analysis c(min, max) in Hz
#' @param threshold Threshold for power spectrum (default: 10)
#' @param method Entropy method ("weiner" or "shannon")
#' @param normalize Logical, whether to normalize entropy values
#' @param plot Logical, whether to plot results (default: TRUE)
#'
#' @return Matrix with columns: time and spectral entropy values
#'
#' @keywords internal
#' @noRd
sh_single_row <- function(segment_row,
                         wav_dir = NULL,
                         wl = 512,
                         wn = "hanning",
                         ovlp = 50,
                         fftw = TRUE,
                         freq_range = c(500, 15000),
                         threshold = 10,
                         method = c("weiner", "shannon"),
                         normalize = FALSE,
                         plot = TRUE) {

  # Match method argument
  method <- match.arg(method)

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

  # Calculate continuous entropy using csh2 function
  sh <- csh2(wave = wv,
            f = wv@samp.rate,
            channel = 1,
            wl = wl,
            wn = wn,
            ovlp = ovlp,
            fftw = fftw,
            freq_range = freq_range,
            threshold = threshold,
            method = method,
            normalize = normalize)

  # Plotting if requested
  if(plot) {
    # Generate spectrogram
    visualize_song.default(x = sound_path,
                           start_time_in_second = segment_row$start_time,
                           end_time_in_second = segment_row$end_time,
                           fft_window_size = wl,
                           overlap = ovlp/100,
                           dark_mode = TRUE,
                           legend = FALSE,
                           keep.par = TRUE,
                           verbose = FALSE)

    # Adjust time values for plotting
    sh[,1] <- sh[,1] + segment_row$start_time

    if (!"duration" %in% names(segment_row)) {
      segment_row$duration <- segment_row$end_time - segment_row$start_time
    }
    offset <- segment_row$duration * 0.035

    # Determine y-axis limits based on normalization
    if(normalize) {
      ylim <- c(0, 1)
    } else {
      ylim <- range(sh[,2], na.rm = TRUE)
      # Add some padding to y-limits
      ylim_range <- diff(ylim)
      ylim <- ylim + c(-0.05, 0.05) * ylim_range
    }

    # Overlay entropy trace
    par(new = TRUE)
    plot(sh[,1], sh[,2],
         type = "l",
         col = "white",
         axes = FALSE,
         xlab = "",
         ylab = "",
         xlim = c(segment_row$start_time + offset, segment_row$end_time - offset),
         ylim = ylim)

    # Add a right-side axis for the entropy values
    axis(side = 4, col = "white", col.axis = "white")

    # # Add right-side label for entropy
    # mtext(paste(tools::toTitleCase(method), "Entropy",
    #             ifelse(normalize, "(Normalized)", "")),
    #       side = 4, col = "white", line = 2)

    # Add title indicating entropy type and normalization
    title_text <- sprintf("Spectrogram w/ %s Entropy Overlay (%s)",
                          tools::toTitleCase(method),
                          ifelse(normalize, "Normalized", "Non-normalized"))
    title(main = title_text)
  }

  return(sh)
}

#' Modified version of seewave::csh using getH
#'
#' @param wave Wave object or time series
#' @param f Sampling frequency (Hz)
#' @param channel Channel to analyze (1 or 2)
#' @param wl Window length for FFT
#' @param wn Window name for FFT ("hanning", "hamming", etc.)
#' @param ovlp Overlap between windows (%)
#' @param fftw Logical, use FFTW or not
#' @param freq_range Frequency range for analysis c(min, max) in Hz
#' @param threshold Amplitude threshold
#' @param method Entropy method ("weiner" or "shannon")
#' @param normalize Logical, whether to normalize entropy values
#'
#' @return Matrix with columns: time and entropy values
#'
#' @keywords internal
#' @noRd
csh2 <- function(wave,
                 f,
                 channel = 1,
                 wl = 512,
                 wn = "hanning",
                 ovlp = 50,
                 fftw = FALSE,
                 freq_range = NULL,
                 threshold = NULL,
                 method = c("weiner", "shannon"),
                 normalize = FALSE) {

  # Match method argument
  method <- match.arg(method)

  # Input processing
  input <- seewave::inputw(wave = wave, f = f, channel = channel)
  wave <- input$w
  f <- input$f
  rm(input)

  # Apply threshold if specified
  if (!is.null(threshold)) {
    wave <- seewave::afilter(wave = wave, f = f, threshold = threshold, plot = FALSE)
  }

  # Calculate spectrogram
  n <- nrow(wave)
  step <- seq(1, n + 1 - wl, wl - (ovlp * wl/100))
  z <- seewave::stdft(wave = wave, f = f, wl = wl, zp = 0, step = step, wn = wn, fftw = fftw)

  # Handle frequency range
  if (!is.null(freq_range)) {
    # Calculate frequency vector
    freq <- seq(0, f/2, length.out = nrow(z))
    # Find indices within frequency range
    freq_idx <- which(freq >= freq_range[1] & freq <= freq_range[2])
    z <- z[freq_idx,]
  }

  # Calculate entropy for each window using getH
  h <- apply(z, MARGIN = 2,
             FUN = function(x) getH(x,
                                    method = method,
                                    normalize = normalize))

  # Create time vector
  t <- seq(0, n/f, length.out = length(step))

  # Replace Inf/-Inf/NaN with NA
  h[!is.finite(h)] <- NA

  return(cbind(t, h))
}

#' Internal function to calculate entropy
#'
#' @description
#' Returns Weiner or Shannon entropy of an input vector such as the spectrum of a sound.
#' This is a modified version of soundgen::getEntropy function.
#'
#' @param x numeric vector
#' @param method character: "weiner" or "shannon"
#' @param normalize logical: whether to normalize the output to [0,1]
#' @param convertNonPositive numeric: value to replace non-positive values
#'
#' @return
#' Numeric value representing entropy. Returns NA if sum(x) = 0
#'
#' @references
#' - Tchernichovski, O., Nottebohm, F., Ho, C. E., Pesaran, B., & Mitra, P. P. (2000).
#'   A procedure for an automated measurement of song similarity.
#'   Animal Behaviour, 59(6), 1167-1176.
#' - Weiner entropy implementation follows SAP2011 algorithm
#' - Original soundgen package: Anikin, A. (2019). Soundgen: An open-source tool
#'   for synthesizing nonverbal vocalizations. Behavior Research Methods, 51(2), 778-792.
#'
#' @keywords internal
#' @noRd
getH <- function(x,
                 method = c("weiner", "shannon"),
                 normalize = FALSE,
                 convertNonPositive = 1e-10) {
  # match the method argument
  method <- match.arg(method)

  sum_x = sum(x)
  if (sum_x == 0)
    return(NA)
  x[x <= 0] = convertNonPositive

  if (method == "weiner") {
    # SAP2011 log-scaled Wiener entropy (default)
    sum_log = sum(log(x))
    log_sum = log(sum_x/length(x))
    entropy = sum_log/length(x) - log_sum

    # Normalize if requested (between 0 and 1)
    if (normalize) {
      entropy = exp(entropy)  # converts log-scaled to original scale
    }
  }
  else if (method == "shannon") {
    p = x/sum_x
    if (normalize) {
      entropy = -sum(p * log2(p)/log(length(p)) * log(2))
    }
    else {
      entropy = -sum(p * log2(p))
    }
  }
  else {
    stop("Implemented entropy methods: \"shannon\" or \"weiner\"")
  }
  entropy
}


# Refine spectral entropy measurement -------------------------------------
#' Refine Spectral Entropy Using Temporal Template
#'
#' @description
#' Refines spectral entropy measurements by creating a temporal template based on pitch goodness.
#' The template is created using either quantile thresholding or Hidden Markov Model (HMM) segmentation
#' of the pitch goodness values from a reference label.
#'
#' @param x A Sap object containing spectral entropy and pitch goodness measurements
#' @param segment_type Character, type of segments to analyze: "motifs", "syllables", or "segments"
#' @param reference_label Character, label to use as reference for template creation
#' @param matrix Character, type of entropy matrix to refine: "weiner" or "shannon"
#' @param method Character, method for template creation: "quantile" or "hmm"
#' @param minimal_duration Numeric, minimum duration (in ms) for segments (default: 20)
#' @param split_dips Logical, whether to split segments at local minima (default: TRUE)
#' @param quantile_threshold Numeric, threshold for quantile method (default: 0.5)
#' @param random_seed Integer, seed for reproducibility in HMM (default: 222)
#' @param plot Logical, whether to plot results (default: TRUE)
#' @param plot_entropy_lim Numeric vector of length 2, limits for entropy plot
#' @param color_palette Function, custom color palette for plotting
#' @param stats Logical, whether to calculate segment statistics (default: TRUE)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param ... Additional arguments passed to plotting functions
#'
#' @details
#' The function creates a temporal template based on pitch goodness values from a reference label,
#' which is then used to filter the spectral entropy measurements. The template can be created
#' using either a quantile threshold or HMM segmentation. The filtered entropy values are stored
#' in the Sap object and can be visualized as a heatmap with the template overlay.
#'
#' @return Invisibly returns the modified Sap object with refined spectral entropy data
#'
#' @examples
#' \dontrun{
#' # Refine Wiener entropy using quantile method
#' sap <- refine_sh(sap,
#'                  segment_type = "motifs",
#'                  reference_label = "a",
#'                  matrix = "weiner",
#'                  method = "quantile",
#'                  plot = TRUE)
#'
#' # Refine Shannon entropy using HMM method
#' sap <- refine_sh(sap,
#'                  segment_type = "motifs",
#'                  reference_label = "a",
#'                  matrix = "shannon",
#'                  method = "hmm",
#'                  plot = TRUE)
#' }
#'
#' @seealso
#' \code{\link{refine_FF}} for refining fundamental frequency
#' \code{\link{spectral_entropy}} for calculating spectral entropy
#'
#' @references
#' For HMM method:
#' Visser, I., & Speekenbrink, M. (2010). depmixS4: An R Package for Hidden Markov Models.
#' Journal of Statistical Software, 36(7), 1-21.
#'
#' @export
refine_sh <- function(x,
                      segment_type = c("motifs", "syllables", "segments"),
                      reference_label,
                      matrix = c("weiner", "shannon"),
                      method = c("quantile", "hmm"),
                      minimal_duration = 20, # in ms
                      split_dips = TRUE,
                      quantile_threshold = 0.5,
                      random_seed = 222,
                      plot = TRUE,
                      plot_entropy_lim = NULL,
                      color_palette = NULL,
                      stats = TRUE,
                      verbose = TRUE,
                      ...) {

  # Input validation
  segment_type <- match.arg(segment_type)
  matrix <- match.arg(matrix)
  method <- match.arg(method)

  # Check for depmixS4 if hmm method is selected
  if (method == "hmm" && !requireNamespace("depmixS4", quietly = TRUE)) {
    stop("Package 'depmixS4' is needed for HMM method. Please install it with:\n",
         "install.packages('depmixS4')",
         call. = FALSE)
  }

  # Check if input is a Sap object
  if (!inherits(x, "Sap")) {
    stop("Input must be a Sap object")
  }

  # Update feature type by removing 's' from segment_type
  feature_type <- sub("s$", "", segment_type)

  # Get the relevant matrices
  entropy_matrix <- x$features[[feature_type]][[paste0(matrix, "_entropy")]]
  goodness_matrix <- x$features[[feature_type]]$pitch_goodness

  # Validate matrices
  if (is.null(entropy_matrix) || is.null(goodness_matrix)) {
    stop("Both spectral entropy and pitch goodness matrices must exist")
  }
  if (!identical(dim(entropy_matrix), dim(goodness_matrix))) {
    stop("Matrices must have identical dimensions")
  }
  if (!reference_label %in% colnames(entropy_matrix)) {
    stop("Reference label not found in data")
  }

  # Calculate time parameters
  time_window <- attr(entropy_matrix, "time_window")
  time_step <- time_window / nrow(entropy_matrix)
  min_samples <- ceiling(minimal_duration / (time_step * 1000))

  # Get reference label columns and average them
  ref_cols <- which(colnames(entropy_matrix) == reference_label)
  ref_goodness <- rowMeans(goodness_matrix[, ref_cols, drop = FALSE])

  # Core segmentation logic
  template <- switch(method,
                     "quantile" = {
                       quantile_val <- quantile(ref_goodness, quantile_threshold, na.rm = TRUE)
                       ref_goodness > quantile_val
                     },

                     "hmm" = {
                       if (!is.null(random_seed)) set.seed(random_seed)

                       hmm_data <- data.frame(goodness = ref_goodness)
                       hmm_model <- depmixS4::depmix(
                         goodness ~ 1,
                         data = hmm_data,
                         nstates = 2,
                         family = gaussian(),
                         instart = c(0.5, 0.5),
                         trstart = matrix(c(0.9, 0.1, 0.1, 0.9), ncol = 2),
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

  # Create filtered entropy matrix
  filtered_entropy <- entropy_matrix
  filtered_entropy[!template, ] <- NA
  attributes(filtered_entropy) <- attributes(entropy_matrix)
  attr(filtered_entropy, "temporal_template") <- template

  # Update the Sap object
  x$features[[feature_type]][[paste0("filtered_", matrix, "_entropy")]] <- filtered_entropy

  if (stats) {
    # Calculate segment statistics
    segment_stat <- calculate_segment_stats(feature_matrix = filtered_entropy,
                                            template = template,
                                            time_step = time_step)

    x$features[[feature_type]][[paste0("stats_", matrix, "_entropy")]] <- segment_stat
  }

  # Plot
  if (plot) {
    # Save original margins
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    # Adjust margins
    par(mar = c(5, 4, 4, 4.5))

    heatmap <- spectral_entropy.matrix(
      entropy_matrix,
      plot_entropy_lim = plot_entropy_lim,
      color_palette = color_palette,
      main = paste("Spectral Entropy (", matrix, ") with", method, "Template")
    )

    # Overlay template
    lattice::trellis.focus("panel", 1, 1, highlight = FALSE)
    template_pos <- which(template)
    if(length(template_pos) > 0) {
      lattice::panel.rect(
        xleft = template_pos - 0.5,
        xright = template_pos + 0.5,
        ybottom = 0.5,
        ytop = ncol(filtered_entropy) + 0.5,
        border = NA,
        col = rgb(1, 1, 1, 0.6)  # white with 60% opacity #rgb(0, 1, 0, 0.5)  # Cyan overlay
      )
    }
    lattice::trellis.unfocus()
  }

  # Print access information
  if(verbose) {
    cat("\nSpectral Entropy refinement completed successfully!")
    cat(sprintf("\nAccess filtered entropy matrix via: sap$features$%s$filtered_%s_entropy",
                feature_type, matrix))
    cat(sprintf("\nAccess template via: attributes(sap$features$%s$filtered_%s_entropy)$temporal_template",
                feature_type, matrix))

    if(stats) {
      cat(sprintf("\nAccess segment statistics via: sap$features$%s$stats_%s_entropy",
                  feature_type, matrix))
    }
  }

  invisible(x)
}
