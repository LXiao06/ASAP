# Analyze Spectral Features -------------------------------------------------------
# Update date : Mar. 4, 2026

#' Analyze Spectral Features of Audio Segments
#'
#' @description
#' Calculates comprehensive spectral features from audio segments with support for
#' parallel processing and organized data structures.
#'
#' @param x An object to analyze, either a data frame or SAP object
#' @param wav_dir Directory containing WAV files
#' @param cores Number of cores for parallel processing
#' @param wl Window length for spectral analysis (default: 512)
#' @param ovlp Overlap percentage (0-100) (default: 50)
#' @param wn Window name ("hanning", "hamming", etc.)
#' @param fftw Logical, use FFTW or not (default: TRUE)
#' @param freq_range Frequency range c(min, max) in kHz
#' @param threshold Threshold for frequency tracking (default: 15)
#' @param fsmooth Frequency smoothing parameter (default: 0.1)
#' @param fast Whether to skip peak frequency calculation (default: TRUE)
#' @param amp_normalize Waveform amplitude normalization before spectral extraction:
#'   one of "none", "peak", or "rms" (default: "none")
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables', 'bouts', 'segments')
#' @param indices For SAP objects: Optional row indices of the selected segment_type to process.
#' @param sample_percent For SAP objects: Percentage of segments to sample
#' @param balanced For SAP objects: Whether to balance groups across labels
#' @param labels For SAP objects: Specific labels to include
#' @param seed For SAP objects: Random seed for sampling (default: 222)
#' @param export_csv If TRUE and \code{indices} exactly match the latest export
#'   \code{exported_indices} for the corresponding segment type, merges stored
#'   export metadata with spectral features and writes a merged CSV.
#' @param csv_filename Output CSV filename. Must be provided if \code{export_csv = TRUE}.
#' @param output_dir Directory to save the CSV. If NULL, defaults to the previous
#'   export's output directory, or the working directory if unavailable.
#' @param time_match_digits Integer digits for matching \code{start_time}/\code{end_time}
#'   when merging metadata and spectral features (default: 3).
#' @param verbose For SAP objects: Whether to print progress messages
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For data frames:
#' \itemize{
#'   \item Requires columns: filename, start_time, end_time
#'   \item Processes segments in parallel
#'   \item Calculates comprehensive spectral features
#'   \item Returns combined results
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Supports multiple segment types
#'   \item Optional balanced sampling
#'   \item Stores results in features slot
#'   \item Preserves segment metadata
#' }
#'
#' @return
#' For default method: A data frame containing spectral features for all segments
#' For SAP objects: Updated SAP object with spectral features stored in features slot
#'
#' @examples
#' \dontrun{
#' # Analyze segments from data frame
#' features <- analyze_spectral(segments,
#'   wav_dir = "path/to/wavs",
#'   cores = 4,
#'   freq_range = c(1, 10)
#' )
#'
#' # Basic analysis from SAP object
#' sap_obj <- analyze_spectral(sap_object,
#'   segment_type = "motifs",
#'   cores = 4
#' )
#'
#' # Balanced sampling with specific labels
#' sap_obj <- analyze_spectral(sap_object,
#'   segment_type = "syllables",
#'   sample_percent = 80,
#'   balanced = TRUE,
#'   labels = c("a", "b")
#' )
#'
#' # Custom spectral parameters
#' sap_obj <- analyze_spectral(sap_object,
#'   segment_type = "motifs",
#'   wl = 1024,
#'   ovlp = 75,
#'   freq_range = c(2, 8)
#' )
#'
#' # With waveform RMS normalization before spectral extraction
#' sap_obj <- analyze_spectral(sap_object,
#'   segment_type = "motifs",
#'   amp_normalize = "rms"
#' )
#' }
#'
#' @rdname analyze_spectral
#' @export
analyze_spectral <- function(x, ...) {
  UseMethod("analyze_spectral")
}

#' @rdname analyze_spectral
#' @export
analyze_spectral.default <- function(x,
                                     wav_dir = NULL,
                                     cores = NULL,
                                     wl = 512,
                                     ovlp = 50,
                                     wn = "hanning",
                                     fftw = TRUE,
                                     freq_range = NULL,
                                     threshold = 15,
                                     fsmooth = 0.1,
                                     fast = TRUE,
                                     amp_normalize = c("none", "peak", "rms"),
                                     ...) {
  # Check required columns
  required_cols <- c("filename", "start_time", "end_time")
  missing_cols <- required_cols[!required_cols %in% names(x)]
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Check if data frame is empty
  if (nrow(x) == 0) {
    stop("Input data frame is empty")
  }

  # Handle wav_dir
  if (is.null(wav_dir) && is.null(attr(x, "wav_dir"))) {
    stop("wav_dir must be provided either as argument or attribute")
  }

  amp_normalize <- .parse_amp_normalize(amp_normalize)

  # Function to process a single row
  process_row <- function(i) {
    spectral_analysis(x[i, ],
      wav_dir = wav_dir,
      wl = wl,
      ovlp = ovlp,
      wn = wn,
      fftw = fftw,
      freq_range = freq_range,
      threshold = threshold,
      fsmooth = fsmooth,
      fast = fast,
      amp_normalize = amp_normalize,
      ...
    )
  }

  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
    cores <- max(1, cores)
  }

  cat(sprintf(
    "\nProcessing %d audio segments using %d cores.\n",
    nrow(x), cores
  ))

  # Choose parallel processing method based on system and cores
  if (fftw) ensure_pkgs("fftw")

  results <- parallel_apply(
    seq_len(nrow(x)),
    process_row,
    cores = cores
  )

  # Combine results
  results <- do.call(rbind, results[!sapply(results, is.null)])

  if (nrow(results) == 0) {
    warning("No valid results produced")
    return(NULL)
  }

  return(results)
}

#' @rdname analyze_spectral
#' @export
analyze_spectral.Sap <- function(x,
                                 segment_type = c("motifs", "syllables", "bouts", "segments"),
                                 indices = NULL,
                                 sample_percent = NULL,
                                 balanced = FALSE,
                                 labels = NULL,
                                 seed = 222,
                                 export_csv = FALSE,
                                 csv_filename = NULL,
                                 output_dir = NULL,
                                 time_match_digits = 3,
                                 cores = NULL,
                                 wl = 512,
                                 ovlp = 50,
                                 wn = "hanning",
                                 fftw = TRUE,
                                 freq_range = NULL,
                                 threshold = 15,
                                 fsmooth = 0.1,
                                 fast = TRUE,
                                 amp_normalize = c("none", "peak", "rms"),
                                 verbose = TRUE,
                                 ...) {
  if (verbose) message(sprintf("\n=== Starting Spectral Features Analysis ===\n"))

  # Input validation
  amp_normalize <- .parse_amp_normalize(amp_normalize)
  segment_type <- match.arg(segment_type)
  segments_df <- x[[segment_type]]

  if (!inherits(segments_df, "segment") || nrow(segments_df) == 0) {
    stop("No segments found in the specified segment type")
  }

  # Subset specific indices if provided
  if (!is.null(indices)) {
    if (!is.numeric(indices) || anyNA(indices)) {
      stop("indices must be a numeric vector without NA values")
    }
    if (any(indices < 1) || any(indices > nrow(segments_df))) {
      stop(sprintf("indices must be between 1 and %d", nrow(segments_df)))
    }
    segments_df <- segments_df[indices, , drop = FALSE]
  } else {
    indices <- seq_len(nrow(segments_df))
  }

  segments_df$.source_index <- as.integer(indices)

  # Select and balance segments
  segments_df <- select_segments(segments_df,
    labels = labels,
    balanced = balanced,
    sample_percent = sample_percent,
    seed = seed
  )

  # Process segments using default method
  if (fftw) ensure_pkgs("fftw")

  result <- analyze_spectral.default(segments_df,
    wav_dir = x$base_path,
    cores = cores,
    wl = wl,
    ovlp = ovlp,
    wn = wn,
    fftw = fftw,
    freq_range = freq_range,
    threshold = threshold,
    fsmooth = fsmooth,
    fast = fast,
    amp_normalize = amp_normalize,
    ...
  )

  # Add attributes to result
  if (".source_index" %in% names(segments_df)) {
    attr(result, "segment_indices") <- as.integer(unique(segments_df$.source_index))
  } else {
    attr(result, "segment_indices") <- which(x[[segment_type]]$filename %in% segments_df$filename &
      x[[segment_type]]$start_time %in% segments_df$start_time &
      x[[segment_type]]$end_time %in% segments_df$end_time)
  }
  attr(result, "segment_type") <- segment_type

  # Update SAP object's features
  feature_type <- sub("s$", "", segment_type) # Remove 's' from end
  x$features[[feature_type]][["spectral_feature"]] <- result

  # Add message about data access
  message(sprintf("Spectral features have been stored in x$features$%s$spectral_feature", feature_type))

  # Optional merge+export with latest export metadata
  if (isTRUE(export_csv)) {
    if (is.null(csv_filename)) {
      stop("csv_filename must be provided when export_csv is TRUE")
    }

    if (verbose) {
      message(sprintf("Exporting %s spectral CSV from selected segments...", segment_type))
    }

    # Extract the correct export list based on segment type
    export_list_name <- switch(segment_type,
      "motifs" = "motif_clip_exports",
      "bouts" = "bout_clip_exports",
      "syllables" = "syllable_clip_exports",
      "segments" = "segment_clip_exports"
    )

    exports <- if (!is.null(export_list_name)) x$misc[[export_list_name]] else NULL

    if (!is.null(exports) && length(exports) > 0) {
      last_export <- exports[[length(exports)]]
      exp_idx <- last_export$exported_indices

      if (!is.null(indices) && !is.null(exp_idx) && length(exp_idx) > 0) {
        same_idx <- setequal(as.integer(indices), as.integer(exp_idx))
        if (same_idx) {
          metadata_df <- last_export$metadata_df
          if (is.null(metadata_df) || !is.data.frame(metadata_df) || nrow(metadata_df) == 0) {
            warning(sprintf("Latest %s export metadata is empty; skipping merge.", segment_type))
          } else {
            # Keep clip-export normalization for provenance
            if ("amp_normalize" %in% names(metadata_df) &&
                !"clip_amp_normalize" %in% names(metadata_df)) {
              metadata_df$clip_amp_normalize <- metadata_df$amp_normalize
            }
            metadata_df$amp_normalize <- amp_normalize

            # Resolve output directory
            if (!is.null(output_dir)) {
              out_dir <- output_dir
            } else if (!is.null(last_export$output_dir)) {
              out_dir <- last_export$output_dir
            } else {
              out_dir <- getwd()
            }

            if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
            out_csv <- file.path(out_dir, csv_filename)

            merged_df <- .merge_metadata_with_spectral(
              metadata_df = metadata_df,
              spectral_df = result,
              output_file = out_csv,
              time_digits = time_match_digits,
              verbose = FALSE
            )

            x$features[[feature_type]]$merged_export_spectral <- merged_df

            if (verbose) {
              message(sprintf(
                "Exported %s spectral CSV: %d/%d rows to %s",
                segment_type,
                nrow(merged_df),
                nrow(metadata_df),
                normalizePath(out_csv, mustWork = TRUE)
              ))
            }
          }
        } else {
          warning("indices do not match latest exported_indices; skipping merge.")
        }
      }
    } else {
      warning(sprintf("No previous exports found for segment_type '%s'; skipping merge.", segment_type))
    }
  }

  # Return updated SAP object invisibly
  invisible(x)
}

#' Internal Spectral Analysis Function
#'
#' @description
#' Internal function to perform spectral analysis on a single audio segment.
#'
#' @param x Single row data frame with segment information
#' @param wav_dir Path to WAV files directory
#' @param wl Window length for analysis
#' @param ovlp Overlap percentage
#' @param wn Window name
#' @param freq_range Frequency range
#' @param threshold Detection threshold
#' @param fsmooth Smoothing parameter
#' @param fast Skip peak frequency calculation
#' @param amp_normalize Waveform amplitude normalization before spectral extraction:
#'   one of "none", "peak", or "rms" (default: "none")
#' @param ... Additional arguments
#'
#' @return
#' Data frame with spectral features
#'
#' @keywords internal
spectral_analysis <- function(x,
                              wav_dir = NULL,
                              wl = 512,
                              ovlp = 50,
                              wn = "hanning",
                              fftw = TRUE,
                              freq_range = NULL,
                              threshold = 15,
                              fsmooth = 0.1,
                              fast = TRUE,
                              amp_normalize = c("none", "peak", "rms"),
                              ...) {
  # Check if input is a single row data frame
  if (!is.data.frame(x) || nrow(x) != 1) {
    stop("Input must be a single row data frame")
  }

  # Check required columns
  required_cols <- c("filename", "start_time", "end_time")
  missing_cols <- required_cols[!required_cols %in% names(x)]
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # # If wav_dir provided, add as attribute
  # if (!is.null(wav_dir) && is.null(attr(x, "wav_dir"))) {
  #   attr(x, "wav_dir") <- normalizePath(wav_dir, mustWork = FALSE)
  # }

  # Get file path
  sound_path <- construct_wav_path(x, wav_dir = wav_dir)

  # Read wave file
  wave <- tuneR::readWave(sound_path,
    from = x$start_time,
    to = x$end_time,
    units = "seconds"
  )

  amp_normalize <- .parse_amp_normalize(amp_normalize)
  if (amp_normalize != "none") {
    wave <- .normalize_wave_amplitude(
      wv = wave,
      method = amp_normalize,
      target_rms = 0.1,
      eps = 1e-8
    )
  }

  # Validate audio data
  if (length(wave@left) < 7) {
    warning("Too few samples in audio segment")
    return(NULL)
  }

  # Set frequency bounds
  if (is.null(freq_range)) {
    freq_range <- c(0, floor(wave@samp.rate / 2000))
  } else {
    if (length(freq_range) != 2 || !is.numeric(freq_range)) {
      stop("freq_range must be NULL or numeric vector of length 2")
    }
    if (freq_range[1] < 0 || freq_range[1] >= freq_range[2]) {
      stop("Invalid freq_range values")
    }
    # Adjust if frequency is too high for sampling rate
    if (freq_range[2] > floor(wave@samp.rate / 2000)) {
      freq_range[2] <- floor(wave@samp.rate / 2000)
    }
  }

  # Adjust window lengths
  wl_adjusted <- ifelse(wl >= length(wave@left), length(wave@left) - 1, wl)
  wl_track <- ifelse(wl >= length(wave@left), length(wave@left) - 2, wl)

  # Calculate spectrogram
  songspec <- suppressWarnings(
    seewave::spec(wave,
      f = wave@samp.rate,
      plot = FALSE,
      wl = wl_adjusted,
      wn = wn,
      flim = freq_range
    )
  )

  if (is.null(songspec)) {
    return(NULL)
  }

  # Calculate basic spectral properties
  analysis <- specprop_wrblr_int(
    spec = songspec,
    f = wave@samp.rate,
    flim = freq_range,
    plot = FALSE
  )

  # Calculate time-based features
  m <- sspectro(wave,
    f = wave@samp.rate, wl = wl_adjusted,
    ovlp = ovlp, wn = wn
  )
  if (!is.matrix(m)) m <- as.matrix(m)

  # Calculate frequency limits in matrix indices.
  # fl values are real-valued; floor them to integers and clamp to [0, nrow(m)-1]
  # so that (fl[1]:fl[2]) + 1 stays within [1, nrow(m)].
  # Without clamping, freq_range[2] == Nyquist gives fl[2] == nrow(m), which
  # makes the +1 index one row past the end of the matrix.
  fl <- floor(freq_range * nrow(m) * 2000 / wave@samp.rate)
  fl[1] <- max(0L, fl[1])
  fl[2] <- min(nrow(m) - 1L, fl[2])
  # Extract relevant frequency range
  m <- m[(fl[1]:fl[2]) + 1L, , drop = FALSE]
  # Ensure matrix format if reduced to vector
  if (is.vector(m)) {
    if (length(m) == 0) {
      warning("No data in specified frequency range")
      return(NULL)
    }
    m <- t(as.matrix(m))
  }

  # Calculate temporal features
  time <- seq(0, length(wave) / wave@samp.rate, length.out = ncol(m))
  t.cont <- apply(m, MARGIN = 2, FUN = sum)
  t.cont <- t.cont / sum(t.cont)
  t.cont.cum <- cumsum(t.cont)
  t.quartiles <- sapply(
    c(0.25, 0.5, 0.75),
    function(x) time[length(t.cont.cum[t.cont.cum <= x]) + 1]
  )

  # Calculate peak time with checks
  peakt <- if (all(t.cont == 0)) {
    NA # If all values are zero
  } else {
    max_indices <- which(t.cont == max(t.cont))
    if (length(max_indices) > 1) {
      # If multiple peaks, take the median position
      time[round(median(max_indices))]
    } else {
      time[max_indices]
    }
  }

  # Calculate frequency tracking with smoothing
  freq_track <- track_harmonic(
    wave = wave,
    f = wave@samp.rate,
    wl = wl_track,
    ovlp = ovlp,
    threshold = threshold,
    bandpass = freq_range * 1000,
    fsmooth = fsmooth,
    plot = FALSE,
    fftw = fftw,
    dfrq = TRUE,
    adjust.wl = TRUE
  )[, 2]

  # Calculate frequency metrics
  freq_metrics <- calculate_freq_metrics(
    freq_track = freq_track,
    freq_range = freq_range,
    start_time = x$start_time,
    end_time = x$end_time,
    wave = wave,
    songspec = songspec,
    wl = wl_adjusted,
    wn = wn,
    fsmooth = fsmooth,
    threshold = threshold,
    ovlp = ovlp,
    fast = fast
  )

  # Combine results
  results <- data.frame(
    filename = x$filename
  )

  # Add day_post_hatch and label columns if they exist in x
  if ("day_post_hatch" %in% names(x)) {
    results$day_post_hatch <- x$day_post_hatch
  }
  if ("label" %in% names(x)) {
    results$label <- x$label
  }

  # Add remaining columns
  results <- cbind(results, data.frame(
    start_time = x$start_time,
    end_time = x$end_time,
    duration = round(x$end_time - x$start_time, digits = 1),
    meanfreq = analysis$mean / 1000,
    sd = analysis$sd / 1000,
    freq.median = analysis$median / 1000,
    freq.Q25 = analysis$Q25 / 1000,
    freq.Q75 = analysis$Q75 / 1000,
    freq.IQR = analysis$IQR / 1000,
    time.median = t.quartiles[2],
    time.Q25 = t.quartiles[1],
    time.Q75 = t.quartiles[3],
    time.IQR = t.quartiles[3] - t.quartiles[1],
    peakt = time[which.max(t.cont)],
    skew = analysis$skewness,
    kurt = analysis$kurtosis,
    sp.ent = analysis$sh,
    time.ent = seewave::th(cbind(time, t.cont)),
    entropy = analysis$sh * seewave::th(cbind(time, t.cont)),
    sfm = analysis$sfm
  ))

  # Add frequency metrics
  results <- cbind(results, freq_metrics)
  row.names(results) <- NULL

  return(results)
}

#' Internal Functions from warbleR Package
#'
#' @name spectral_analysis
#' @description
#' Internal functions adapted from the warbleR package for spectral analysis.
#' These functions are modified versions of the original warbleR code.
#'
#' @details
#' Original source: warbleR package
#' Citation: Araya-Salas, M. and Smith-Vidaurre, G. (2017),
#' warbleR: an r package to streamline analysis of animal acoustic signals.
#' Methods Ecol Evol. 8, 184-191.
#'
#' @keywords internal
NULL

#' Calculate Frequency Metrics
#' @noRd
#' @keywords internal
calculate_freq_metrics <- function(freq_track,
                                   freq_range,
                                   start_time,
                                   end_time,
                                   wave,
                                   songspec,
                                   wl,
                                   wn,
                                   fsmooth,
                                   threshold,
                                   ovlp,
                                   fast = TRUE) {
  # Clean frequency track data
  y <- freq_track[!is.na(freq_track)]
  y <- y[y >= freq_range[1] & y <= freq_range[2] & y != 0]

  # Calculate basic frequency metrics
  if (length(y) > 0) {
    meandom <- mean(y, na.rm = TRUE)
    mindom <- min(y, na.rm = TRUE)
    maxdom <- max(y, na.rm = TRUE)
    dfrange <- maxdom - mindom
    startdom <- y[1]
    enddom <- y[length(y)]

    # Calculate modulation index
    modindx <- if (length(y) > 1 & dfrange != 0) {
      sum(abs(diff(y))) / dfrange
    } else {
      1
    }

    # Calculate frequency slope
    dfslope <- (enddom - startdom) / (end_time - start_time)
  } else {
    meandom <- mindom <- maxdom <- dfrange <- startdom <- enddom <- modindx <- dfslope <- NA
  }

  # Calculate mean peak frequency
  frng <- frd_wrblr_int(
    wave = wave,
    wl = wl,
    fsmooth = fsmooth,
    threshold = threshold,
    wn = wn,
    bp = freq_range,
    ovlp = ovlp
  )
  meanpeakf <- if (length(frng$meanpeakf) > 0) frng$meanpeakf else NA

  # Create base results data frame
  results <- data.frame(
    meandom = meandom,
    mindom = mindom,
    maxdom = maxdom,
    dfrange = dfrange,
    modindx = modindx,
    startdom = startdom,
    enddom = enddom,
    dfslope = dfslope,
    meanpeakf = meanpeakf
  )


  # Calculate peak frequency if fast is FALSE
  if (!fast) {
    peakf <- try(
      seewave::fpeaks(songspec,
        f = wave@samp.rate,
        wl = wl,
        nmax = 3,
        plot = FALSE
      )[1, 1],
      silent = TRUE
    )

    # Add peakf only if calculation was successful
    if (!inherits(peakf, "try-error") && length(peakf) > 0) {
      results$peakf <- peakf
    }
  }

  return(results)
}

#' Calculate Spectral Properties (from warbleR)
#' @noRd
#' @keywords internal
specprop_wrblr_int <- function(spec, f = NULL, flim = NULL, ...) {
  fhz <- f

  if (is.null(f)) {
    if (is.vector(spec)) {
      stop("'f' is missing")
    } else if (is.matrix(spec)) {
      f <- spec[nrow(spec), 1] * 2000 * nrow(spec) / (nrow(spec) - 1)
    }
  }

  if (is.matrix(spec)) {
    freq <- spec[, 1]
    freq <- freq * 1000
    spec <- spec[, 2]
  }
  L <- length(spec)
  wl <- L * 2
  if (any(spec < 0)) {
    stop("The frequency spectrum to be analysed should not be in dB")
  }
  if (!is.null(flim)) {
    if (flim[1] < 0 || flim[2] > fhz / 2) {
      stop("'flim' should range between 0 and f/2")
    }
  } else {
    flim <- c(0, (f / 2 - f / wl) / 1000)
  }
  g <- (1000 * wl / 2) / (f / 2 - f / wl)
  spec <- spec[(flim[1] * g):(flim[2] * g)]
  spec <- spec[!is.na(spec)]
  L <- length(spec)
  amp <- spec / sum(spec)
  cumamp <- cumsum(amp)
  freq <- seq(
    from = flim[1] * 1000,
    to = flim[2] * 1000,
    length.out = L
  )
  mean <- sum(amp * freq)
  sd <- sqrt(sum(amp * ((freq - mean)^2)))
  sem <- sd / sqrt(L)
  median <- freq[length(cumamp[cumamp <= 0.5]) + 1]
  mode <- freq[which.max(amp)]
  Q25 <- freq[length(cumamp[cumamp <= 0.25]) + 1]
  Q75 <- freq[length(cumamp[cumamp <= 0.75]) + 1]
  IQR <- Q75 - Q25
  cent <- sum(freq * amp)
  z <- amp - mean(amp)
  w <- sd(amp)
  skew <- (sum(z^3) / (L - 1)) / w^3
  kurt <- (sum(z^4) / (L - 1)) / w^4
  sfm <- seewave::sfm(amp)
  sh <- seewave::sh(amp)
  prec <- f / wl

  results <- list(
    mean = mean,
    sd = sd,
    median = median,
    sem = sem,
    mode = mode,
    Q25 = Q25,
    Q75 = Q75,
    IQR = IQR,
    cent = cent,
    skewness = skew,
    kurtosis = kurt,
    sfm = sfm,
    sh = sh,
    prec = prec
  )

  return(results)
}

#' Calculate Frequency Range Detection (from warbleR)
#' @noRd
#' @keywords internal
frd_wrblr_int <- function(wave,
                          wl = 512,
                          fsmooth = 0.1,
                          threshold = 10,
                          wn = "hanning",
                          bp = NULL,
                          ovlp = 50,
                          dB.threshold = NULL) {
  # get sampling rate
  f <- wave@samp.rate

  if (wl >= length(wave@left)) {
    wl <- length(wave@left) - 1
  }
  if (wl %% 2 != 0) {
    wl <- wl - 1
  }

  # mean spectrum
  if (is.null(dB.threshold)) {
    spc <-
      seewave::meanspec(
        wave,
        plot = FALSE,
        wl = wl,
        f = f,
        wn = wn,
        ovlp = ovlp
      )

    # get frequency windows length for smoothing
    step <- wave@samp.rate / wl / 1000

    # number of samples
    n <- nrow(spc)

    # smoothing parameter
    if (!is.null(fsmooth)) {
      fsmooth <- fsmooth / step

      FWL <- fsmooth - 1

      # smooth
      z <-
        apply(as.matrix(1:(n - FWL)), 1, function(y) {
          sum(spc[y:(y + FWL), 2])
        })
      zf <-
        seq(min(spc[, 1]), max(spc[, 1]), length.out = length(z))
    } else {
      z <- spc[, 2]
      zf <- spc[, 1]
    }

    # remove range outside bp
    if (!is.null(bp)) {
      # if there are complete freq bins within freq range
      if (any(zf > bp[1] & zf < bp[2])) {
        fbins <- which(zf > bp[1] & zf < bp[2])
      } else {
        # select the one that contains the freq range
        fbins <-
          which.max(ifelse(zf - bp[1] > 0, NA, zf - bp[1])):which.max(ifelse(zf - bp[2] > 0, NA, zf - bp[1]))
      }

      z <- z[fbins]
      zf <- zf[fbins]
    }

    # make minimum amplitude 0
    z <- z - min(z)
    z[z < 0] <- 0

    # normalize amplitude from 0 to 1
    if (length(z) > 1) {
      z <- z / max(z)
    }

    # get freqs crossing threshold
    z1 <- rep(0, length(z))
    z1[z > threshold / 100] <- 1
    z2 <- z1[2:length(z1)] - z1[1:(length(z1) - 1)]

    # add 0 to get same length than z
    z2 <- c(0, z2)

    # determine start and end of amplitude hills
    strt <- zf[z2 == 1]
    nd <- zf[z2 == -1]

    # add NAs when some ends or starts where not found
    if (length(strt) != length(nd)) {
      if (z1[1] == 0) {
        nd <- c(nd, NA)
      } else {
        strt <- c(NA, strt)
      }
    }

    if (length(strt) == 1) {
      if (z1[1] == 1 & z1[length(z1)] == 1 & strt > nd) {
        strt <- c(NA, strt)
        nd <- c(nd, NA)
      }
    }
    # substract half a step to calculate mid point between the 2 freq windows in which the threshold has passed
    nd <- nd - (step / 2)
    strt <- strt - (step / 2)

    meanpeakf <- zf[which.max(z)] + (step / 2)
  } else {
    spc <-
      meanspec(
        wave,
        plot = FALSE,
        wl = wl,
        f = f,
        wn = wn,
        ovlp = ovlp,
        dB = "max0",
        dBref = 2 * 10e-5
      )

    # get frequency windows length for smoothing
    step <- wave@samp.rate / wl / 1000

    # number of samples
    n <- nrow(spc)

    # smoothing parameter
    if (!is.null(fsmooth)) {
      fsmooth <- fsmooth / step

      FWL <- fsmooth - 1

      # smooth
      z <-
        apply(as.matrix(1:(n - FWL)), 1, function(y) {
          sum(spc[y:(y + FWL), 2])
        })
      zf <-
        seq(min(spc[, 1]), max(spc[, 1]), length.out = length(z))

      z <-
        (max(spc[, 2]) - min(spc[, 2])) / (max(z) - min(z)) * (z - max(z)) + max(spc[, 2])
    } else {
      z <- spc[, 2]
      zf <- spc[, 1]
    }

    if (!is.null(bp)) {
      # remove range outsde bp
      z <- z[zf > bp[1] & zf < bp[2]]
      zf <- zf[zf > bp[1] & zf < bp[2]]
    }

    z1 <- rep(0, length(z))
    z1[z > max(z) - dB.threshold] <- 1
    z2 <- z1[2:length(z1)] - z1[1:(length(z1) - 1)]

    # add 0 to get same length than z
    z2 <- c(0, z2)

    # determine start and end of amplitude hills
    strt <- zf[z2 == 1]
    nd <- zf[z2 == -1]

    # add NAs when some ends or starts where not found
    if (length(strt) != length(nd)) {
      if (z1[1] == 0) {
        nd <- c(nd, NA)
      } else {
        strt <- c(NA, strt)
      }
    }

    if (length(strt) == 1) {
      if (z1[1] == 1 & z1[length(z1)] == 1 & strt > nd) {
        strt <- c(NA, strt)
        nd <- c(nd, NA)
      }
    }

    # step <- mean(zf[-1] - zf[1:(length(zf) - 1)])

    # substract half a step to calculate mid point between the 2 freq windows in which the threshold has passed
    nd <- nd - (step / 2)
    strt <- strt - (step / 2)
    meanpeakf <- zf[which.max(z)] + (step / 2)
  }

  # fix range
  # only start lower than peakf
  strt <- strt[strt <= meanpeakf]

  # only ends higher than peakf
  nd <- nd[nd >= meanpeakf]

  # get freq range
  min.strt <-
    ifelse(length(strt) == 1, strt, strt[which.min(meanpeakf - strt)])
  max.nd <-
    ifelse(length(nd) == 1, nd, nd[which.min(nd - meanpeakf)])

  if (!any(is.na(c(min.strt, max.nd)))) {
    if (min.strt > max.nd) {
      min.strt <- NA
      max.nd <- NA
    }
  }

  # force nd and strt the same length adding NAs
  if (length(nd) > length(strt)) {
    strt <- c(strt, rep(NA, length(nd) - length(strt)))
  }
  if (length(strt) > length(nd)) {
    nd <- c(nd, rep(NA, length(strt) - length(nd)))
  }

  # save everything in a list
  rl <-
    list(
      frange = data.frame(bottom.freq = min.strt, top.freq = max.nd),
      af.mat = cbind(z, zf),
      meanpeakf = meanpeakf,
      detections = cbind(start.freq = strt, end.freq = nd),
      threshold = ifelse(is.null(dB.threshold), threshold, max(z) - dB.threshold),
      type = ifelse(is.null(dB.threshold), "percentage", "dB")
    )

  # return rl list
  return(rl)
}

#' Track Harmonics in Audio (from warbleR)
#' @noRd
#' @keywords internal
track_harmonic <- function(wave, f, wl = 512, wn = "hanning", ovlp = 0, fftw = FALSE,
                           at = NULL, tlim = NULL, threshold = 10, bandpass = NULL,
                           clip = NULL, plot = TRUE, xlab = "Times (s)", ylab = "Frequency (kHz)",
                           ylim = c(0, f / 2000), adjust.wl = FALSE, dfrq = FALSE, ...) {
  if (length(inputw(wave = wave, f = f)$w) < wl) {
    if (adjust.wl) {
      wl <- length(wave)
    } else {
      stop("number of samples lower than 'wl' (i.e. no enough samples) \n check 'adjust.wl' argument")
    }
  }

  if (!is.null(at) && ovlp != 0) {
    stop("The 'ovlp' argument cannot bue used in conjunction with the arguement 'at'.")
  }
  if (!is.null(clip)) {
    if (clip <= 0 | clip >= 1) {
      stop("'clip' value has to be superior to 0 and inferior to 1")
    }
  }

  input <- inputw(wave = wave, f = f)
  wave <- input$w
  f <- input$f
  rm(input)
  if (!is.null(tlim)) {
    wave <- seewave::cutw(wave, f = f, from = tlim[1], to = tlim[2])
  }

  if (!is.null(threshold)) {
    wave <- seewave::afilter(
      wave = wave, f = f, threshold = threshold,
      plot = FALSE
    )
  }

  n <- nrow(wave)
  if (!is.null(at)) {
    step <- at * f
    N <- length(step)
    if (step[1] <= 0) {
      step[1] <- 1
    }
    if (step[N] + (wl / 2) >= n) {
      step[N] <- n - wl
    }
    x <- c(0, at, n / f)
  } else {
    step <- seq(1, n - wl, wl - (ovlp * wl / 100))
    N <- length(step)
    x <- seq(0, n / f, length.out = N)
  }

  step <- round(step)
  y1 <- seewave::stdft(
    wave = wave, f = f, wl = wl, zp = 0, step = step,
    wn = wn
  )
  if (!is.null(bandpass)) {
    if (length(bandpass) != 2) {
      stop("'The argument 'bandpass' should be a numeric vector of length 2'")
    }
    if (bandpass[1] > bandpass[2]) {
      stop("The first element of 'bandpass' has to be inferior to the second element, i.e. bandpass[1] < bandpass[2]")
    }
    if (bandpass[1] == bandpass[2]) {
      stop("The limits of the bandpass have to be different")
    }
  }

  # lowlimit <- round((wl * bandpass[1])/f)
  # upperlimit <- round((wl * bandpass[2])/f)
  # y1[-(lowlimit:upperlimit), ] <- 0

  # freq values for each freq window   (using mid point of each window)
  freq.val <- ((1:nrow(y1) * f / wl) - (f / (wl * 2)))

  y1[freq.val < bandpass[1] | freq.val > bandpass[2]] <- 0


  if (dfrq) {
    maxi <- apply(y1, MARGIN = 2, FUN = max)
    y2 <- apply(y1, MARGIN = 2, FUN = which.max)
  } else {
    # find peaks close to first dom freq
    maxi <- NULL
    y2 <- NULL

    for (i in seq_len(ncol(y1)))
    {
      # standardize z between 0-1
      z <- y1[, i] / max(y1[, i])
      z <- ifelse(z > threshold / 100, z, 0)

      # choose the maximum amplitude for the firt time window
      if (i == 1) {
        maxi[i] <- max(z)
        y2[i] <- which.max(z)
      } else {
        ensure_pkgs("pracma")
        pks <- pracma::findpeaks(z, npeaks = 5, sortstr = TRUE)[, 1:3]
        if (is.vector(pks)) pks <- matrix(pks, ncol = 3)
        pks[, 3] <- abs(pks[, 2] - y2[i - 1])
        maxi[i] <- pks[which.min(pks[, 3]), 1]
        y2[i] <- pks[which.min(pks[, 3]), 2]
      }
    }
  }

  y2[which(maxi == 0)] <- NA
  if (!is.null(clip)) {
    maxi <- apply(y1, MARGIN = 2, FUN = max)
    y2[which(maxi < clip)] <- NA
  }
  # y <- (f * y2)/(1000 * wl) - f/(1000 * wl)
  y <- freq.val[y2]

  if (!is.null(at)) {
    y <- c(NA, y, NA)
  }

  y <- y / 1000

  if (plot) {
    plot(
      x = x, y = y, xaxs = "i", xlab = xlab, yaxs = "i",
      ylab = ylab, ylim = ylim, ...
    )
    invisible(cbind(x, y))
  } else {
    return(cbind(x, y))
  }
}

#' Internal Functions from seewave Package
#' @noRd
#'
#' @note Adapted from seewave::sspectro by Jerome Sueur et al.
#' @keywords internal
sspectro <- function(wave, f, wl = 512, ovlp = 0, wn = "hanning",
                     norm = TRUE, correction = "none") {
  # You'll need to implement or import these functions
  input <- inputw(wave = wave, f = f)
  wave <- input$w
  f <- input$f
  rm(input)
  n <- nrow(wave)
  step <- seq(1, n + 1 - wl, wl - (ovlp * wl / 100))
  W <- ftwindow(wl = wl, wn = wn, correction = correction)
  z <- apply(as.matrix(step), 1, function(x) {
    Mod(stats::fft(wave[x:(wl +
      x - 1), ] * W))
  })
  z <- z[2:(1 + wl / 2), ]
  if (norm) {
    z <- z / max(z)
  }
  return(z)
}

#' Merge Segment Metadata with Spectral Features
#'
#' @description
#' Internal helper to merge exported segment metadata with newly extracted
#' spectral features before writing to a CSV.
#'
#' @keywords internal
.merge_metadata_with_spectral <- function(metadata_df,
                                                spectral_df,
                                                output_file,
                                                time_digits = 6,
                                                verbose = TRUE) {
  # metadata_df must have all canonical columns (it is the left / authoritative side).
  required_meta <- c("filename", "start_time", "end_time", "duration", "day_post_hatch", "label")
  if (!all(required_meta %in% names(metadata_df))) {
    stop("metadata_df must contain: filename, start_time, end_time, duration, day_post_hatch, label")
  }

  # Coerce spectral_df to a plain data.frame to guard against tibble / matrix
  # edge cases where names() might behave unexpectedly after parallel rbind.
  spectral_df <- as.data.frame(spectral_df)

  # spectral_df needs only filename + time columns for key-building.
  # duration is intentionally excluded: spectral_analysis() rounds it to
  # 1 decimal while metadata_df stores it at full precision, causing key
  # mismatches. filename + start_time + end_time alone uniquely identify a row.
  required_spec <- c("filename", "start_time", "end_time")
  missing_spec <- setdiff(required_spec, names(spectral_df))
  if (length(missing_spec) > 0) {
    stop("spectral_df must contain: ", paste(missing_spec, collapse = ", "))
  }

  # Build the merge key using only filename, start_time, end_time.
  # day_post_hatch / label are omitted (they come from metadata_df's left-join)
  # and duration is omitted (precision mismatch between the two data frames).
  make_key <- function(df) {
    filename <- as.character(df$filename)
    start_key <- sprintf(paste0("%.", time_digits, "f"), round(as.numeric(df$start_time), time_digits))
    end_key <- sprintf(paste0("%.", time_digits, "f"), round(as.numeric(df$end_time), time_digits))
    paste(filename, start_key, end_key, sep = "||")
  }

  metadata_df$.merge_key <- make_key(metadata_df)
  spectral_df$.merge_key <- make_key(spectral_df)

  # Keep metadata columns as canonical output columns; add only spectral-specific columns.
  keep_spec <- setdiff(names(spectral_df), intersect(names(spectral_df), names(metadata_df)))
  keep_spec <- union(".merge_key", keep_spec)
  spectral_reduced <- spectral_df[, keep_spec, drop = FALSE]

  merged <- merge(
    metadata_df,
    spectral_reduced,
    by = ".merge_key",
    all.x = TRUE,
    sort = FALSE
  )
  merged$.merge_key <- NULL

  utils::write.csv(merged, output_file, row.names = FALSE)

  if (verbose) {
    message(sprintf("Merged CSV path: %s", normalizePath(output_file, mustWork = TRUE)))
  }

  merged
}

