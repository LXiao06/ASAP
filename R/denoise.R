# Audio Denoising Module ---------------------------------------------------
# Update date : Mar. 15, 2026

# Denoise ------------------------------------------------------------------

#' Denoise Audio Files
#'
#' @description
#' Removes stationary background noise from audio recordings using
#' spectral-domain techniques. Two methods are available:
#' \itemize{
#'   \item \strong{Spectral Median Subtraction} (\code{"spectral_median"}):
#'     Estimates the per-frequency noise floor as the quantile magnitude across
#'     all time frames, subtracts it from every frame (half-wave rectified),
#'     then reconstructs the waveform via inverse STFT.  Fast and reliable gold
#'     standard for flat, stationary broadband noise.
#'   \item \strong{Spectral Gating} (\code{"spectral_gate"}):
#'     Audacity-style noise reduction.  Builds a per-frequency noise profile
#'     (floor + spread), then applies a smooth sigmoid gate so that
#'     time-frequency bins significantly above the noise floor pass through
#'     nearly unmodified, while bins close to the noise floor are gently
#'     attenuated.  A configurable floor (\code{gate_floor}) prevents any bin
#'     from being completely silenced, avoiding the spectral "twist" artefacts
#'     that hard gating can introduce in bird-song recordings.
#' }
#'
#' @param x An object to process: a WAV file path (character) or a SAP object.
#' @param method Character, the denoising method to use:
#'   \code{"spectral_median"} (default) or \code{"spectral_gate"}.
#' @param output_dir Character, directory for denoised WAV files.
#'   If \code{NULL} (default), a \code{denoised/} subdirectory is created beside
#'   each source file (default method) or under \code{x$base_path} (SAP method).
#' @param overwrite Logical, whether to overwrite existing denoised files
#'   (default: \code{FALSE}).
#' @param wl Integer, STFT window length in samples (default: 256).
#' @param ovlp Integer, STFT overlap in percent (default: 50).
#' @param wn Character, window function passed to \code{seewave::ftwindow}
#'   (default: \code{"hanning"}).
#' @param plot Logical, for default method only. When \code{TRUE} (default),
#'   plots the spectrogram of the denoised output audio using the same
#'   visualization style as \code{\link{segment}}.
#' @param view_window Numeric vector of length 2 (seconds) for the plot
#'   window in the denoised spectrogram. \code{NULL} (default) shows the
#'   full file. Example: \code{c(1, 4)} plots 1s to 4s.
#' @param freq_range Numeric vector of length 2 giving the frequency range
#'   (in kHz) to denoise, e.g. \code{c(0, 10)} (default).  Only frequency bins
#'   within this range are passed through the denoising algorithm; bins outside
#'   the range are left untouched.  Set to \code{NULL} to denoise the entire
#'   spectrum.  Restricting the range speeds up processing and avoids altering
#'   frequency bands that contain no noise.
#'
#' @param noise_quantile Numeric in (0, 1] or \code{NULL}.  The role differs
#'   by method (auto-selected when \code{NULL}):
#'   \itemize{
#'     \item \strong{spectral_median}: quantile of the per-frequency magnitude
#'       distribution across \emph{all} time frames used as the noise floor
#'       estimate.  Auto-default: \code{0.5} (median).
#'     \item \strong{spectral_gate}: fraction of time frames (ranked by total
#'       frame energy, lowest first) selected as \emph{quiet/noise-only}
#'       frames.  The noise floor and spread are computed exclusively from
#'       those frames, so song syllables never contaminate the noise profile.
#'       Auto-default: \code{0.25} (lowest-energy 25\% of frames).
#'   }
#'   Lower values are more conservative; higher values capture more of the
#'   noise distribution.
#'
#' @param gain Numeric \eqn{\ge 0}, over-subtraction factor for
#'   \strong{spectral_median} (default: 1.0).  Values > 1 increase noise
#'   removal at the risk of artefacts.
#'
#' @param gate_threshold Numeric \eqn{\ge 0}, for \strong{spectral_gate}:
#'   number of noise-spread standard deviations above the noise floor at which
#'   the sigmoid gate reaches 50 \% transmission (default: 1.5).
#'   Lower = less aggressive (more signal preserved);
#'   higher = more aggressive (more noise removed).
#' @param gate_smoothing Integer, for \strong{spectral_gate}: half-width (in
#'   frequency bins) of the box-car smoothing applied to the gate mask
#'   (default: 3).  Smoothing reduces sharp spectral edges that can distort
#'   song structure.  Set to 0 to disable.
#' @param gate_floor Numeric in [0, 1), for \strong{spectral_gate}: minimum
#'   gate value applied to every bin (default: 0.1).  A non-zero floor ensures
#'   no frequency bin is completely silenced, preserving tonal continuity and
#'   preventing the "hollow" or "twisted" sound artefacts of hard gating.
#'
#' @param day For SAP objects: Numeric vector of days post-hatch to process.
#'   \code{NULL} (default) = all days.
#' @param indices For SAP objects: Numeric vector of row indices in
#'   \code{x$metadata}.  \code{NULL} = all rows within selected days.
#' @param cores For SAP objects: Number of parallel cores
#'   (default: \code{parallel::detectCores() - 1}).
#' @param update_base_path Logical, for SAP objects only.  When \code{TRUE},
#'   \code{x$base_path} is replaced with the denoised output directory after
#'   processing completes (default: \code{FALSE}).  This makes all subsequent
#'   pipeline steps (\code{detect_template}, \code{export_clips},
#'   \code{extract_envelope}, etc.) automatically read from the denoised files
#'   without any further configuration, since the denoised directory mirrors
#'   the original \code{day_post_hatch/filename} structure exactly.
#' @param verbose Logical, print progress messages (default: \code{TRUE}).
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' \strong{Spectral Median Subtraction} steps:
#' \enumerate{
#'   \item Compute STFT magnitude and phase.
#'   \item Estimate the per-frequency noise floor as the \code{noise_quantile}
#'         quantile of magnitude across all time frames.
#'   \item Subtract \code{gain * noise_floor} from every frame; zero out
#'         negatives (half-wave rectification).
#'   \item Reconstruct the waveform from cleaned magnitude + original phase
#'         via inverse STFT (overlap-add).
#' }
#'
#' \strong{Spectral Gating} steps:
#' \enumerate{
#'   \item Compute STFT magnitude and phase.
#'   \item Rank all time frames by their total energy
#'         (\eqn{E_t = \sum_f |S(t,f)|^2}).  Select the lowest-energy
#'         \code{noise_quantile} fraction as \emph{quiet frames}.  Because
#'         bird-song syllables are brief and high-energy, they are excluded
#'         from this set, so the noise profile is estimated from
#'         background-only frames only.
#'   \item Compute the per-frequency noise floor (\eqn{\mu_f}) as the mean
#'         and spread (\eqn{\sigma_f}) as the SD of those quiet frames.
#'   \item Build a soft sigmoid gate for every time-frequency bin:
#'         \deqn{G(t,f) = \max\!\left(g_{\min},\;
#'         \sigma\!\left(\frac{|S(t,f)| - \mu_f -
#'         \theta\sigma_f}{\sigma_f/4}\right)\right)}
#'         where \eqn{\theta} = \code{gate_threshold} and
#'         \eqn{g_{\min}} = \code{gate_floor}.
#'   \item Optionally smooth the gate mask along the frequency axis
#'         (\code{gate_smoothing}) to reduce spectral edge artefacts.
#'   \item Multiply the gate mask by the original magnitude and reconstruct
#'         the waveform.
#' }
#'
#' @return
#' \describe{
#'   \item{Default method}{Character path to the denoised WAV file (invisibly).}
#'   \item{SAP method}{Updated SAP object with \code{x$denoised_path} set to
#'     the output directory, mirroring the \code{day_post_hatch/filename}
#'     structure of the originals.}
#' }
#'
#' @examples
#' \dontrun{
#' # --- Single WAV file -------------------------------------------------------
#' # Spectral median (default — fast, reliable)
#' clean_path <- denoise("path/to/recording.wav")
#'
#' # Spectral gate — cleaner background with auto noise_quantile = 0.25
#' clean_path <- denoise("path/to/recording.wav",
#'   method = "spectral_gate"
#' )
#'
#' # Override gate aggressiveness
#' clean_path <- denoise("path/to/recording.wav",
#'   method         = "spectral_gate",
#'   gate_threshold = 2.0, # more aggressive
#'   gate_floor     = 0.05, # tighter floor
#'   gate_smoothing = 5L # extra smoothing
#' )
#'
#' # Spectral median — auto noise_quantile = 0.5
#' clean_path <- denoise("path/to/recording.wav",
#'   method = "spectral_median",
#'   gain   = 1.3
#' )
#'
#' # --- SAP object (batch) ---------------------------------------------------
#' sap_obj <- denoise(sap_obj, method = "spectral_median", cores = 4)
#'
#' sap_obj <- denoise(sap_obj,
#'   method         = "spectral_gate",
#'   day            = c(70, 75, 80),
#'   gate_threshold = 1.0,
#'   gate_floor     = 0.1,
#'   cores          = 8
#' )
#' }
#'
#' @seealso \code{\link{detect_template}} for template matching after denoising.
#'
#' @rdname denoise
#' @export
denoise <- function(x, ...) {
  UseMethod("denoise")
}


# Default method (single WAV file) -----------------------------------------

#' @rdname denoise
#' @export
denoise.default <- function(x,
                            method = c(
                              "spectral_median",
                              "spectral_gate"
                            ),
                            output_dir = NULL,
                            overwrite = FALSE,
                            wl = 256L,
                            ovlp = 50L,
                            wn = "hanning",
                            plot = TRUE,
                            view_window = NULL,
                            freq_range = c(0, 10),
                            # -- shared noise-floor param --
                            noise_quantile = NULL,
                            # -- spectral_median params --
                            gain = 1.0,
                            # -- spectral_gate params --
                            gate_threshold = 1.5,
                            gate_smoothing = 3L,
                            gate_floor = 0.1,
                            verbose = TRUE,
                            ...) {
  method <- match.arg(method)

  # Auto-select noise_quantile based on method when not explicitly provided
  if (is.null(noise_quantile)) {
    noise_quantile <- if (method == "spectral_gate") 0.25 else 0.5
  }

  if (!file.exists(x)) stop("File does not exist: ", x)

  # Output directory
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(x), "denoised")
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  out_path <- file.path(output_dir, basename(x))

  if (file.exists(out_path) && !overwrite) {
    if (verbose) cat(sprintf("Skipping (exists): %s\n", basename(out_path)))
    return(invisible(normalizePath(out_path)))
  }

  # Read audio
  wave <- tuneR::readWave(x)
  sr <- wave@samp.rate
  bit_depth <- wave@bit
  samples <- as.numeric(wave@left)

  # STFT
  hop <- as.integer(wl * (1 - ovlp / 100))
  win <- seewave::ftwindow(wl = wl, wn = wn)
  stft <- .forward_stft(samples, wl = wl, hop = hop, win = win)
  mag <- Mod(stft)
  phase <- Arg(stft)

  # Frequency range restriction
  n_bins <- nrow(mag)
  if (!is.null(freq_range)) {
    bin_lo <- max(1L, floor(freq_range[1] * 1000 * wl / sr) + 1L)
    bin_hi <- min(n_bins, ceiling(freq_range[2] * 1000 * wl / sr) + 1L)
    target_rows <- bin_lo:bin_hi
  } else {
    target_rows <- seq_len(n_bins)
  }

  # Denoise only the target frequency bins
  mag_sub <- mag[target_rows, , drop = FALSE]
  clean_sub <- switch(method,
    spectral_median = .spectral_median_subtract(
      mag_sub,
      noise_quantile = noise_quantile,
      gain           = gain
    ),
    spectral_gate = .spectral_gate(
      mag_sub,
      noise_quantile = noise_quantile,
      threshold      = gate_threshold,
      smoothing      = gate_smoothing,
      floor          = gate_floor
    )
  )

  # Recombine: replace only the target rows in the full magnitude
  clean_mag <- mag
  clean_mag[target_rows, ] <- clean_sub

  # Reconstruct
  clean_complex <- clean_mag * exp(1i * phase)
  clean_samples <- .inverse_stft(clean_complex,
    wl = wl, hop = hop, win = win,
    original_length = length(samples)
  )

  if (plot) {
    .plot_denoised_spectrogram(clean_samples,
      sr = sr,
      wl = wl,
      ovlp = ovlp,
      wn = wn,
      method = method,
      view_window = view_window
    )
  }

  max_amp <- 2^(bit_depth - 1) - 1
  clean_int <- as.integer(round(pmax(pmin(clean_samples, 1), -1) * max_amp))

  # Write
  tuneR::writeWave(
    tuneR::Wave(left = clean_int, samp.rate = sr, bit = bit_depth),
    filename = out_path
  )

  if (verbose) {
    cat(sprintf(
      "Denoised (%s): %s -> %s\n",
      method, basename(x), out_path
    ))
  }
  invisible(normalizePath(out_path))
}


# SAP method ---------------------------------------------------------------

#' @rdname denoise
#' @export
denoise.Sap <- function(x,
                        method = c(
                          "spectral_median",
                          "spectral_gate"
                        ),
                        output_dir = NULL,
                        overwrite = FALSE,
                        wl = 256L,
                        ovlp = 50L,
                        wn = "hanning",
                        plot = FALSE,
                        view_window = NULL,
                        freq_range = c(0, 10),
                        noise_quantile = NULL,
                        gain = 1.0,
                        gate_threshold = 1.5,
                        gate_smoothing = 3L,
                        gate_floor = 0.1,
                        day = NULL,
                        indices = NULL,
                        cores = NULL,
                        update_base_path = FALSE,
                        verbose = TRUE,
                        ...) {
  method <- match.arg(method)

  # Auto-select noise_quantile per method when not explicitly set
  if (is.null(noise_quantile)) {
    noise_quantile <- if (method == "spectral_gate") 0.25 else 0.5
  }

  if (verbose) message(sprintf("\n=== Starting Audio Denoising (%s) ===\n", method))
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  # Output root directory
  if (is.null(output_dir)) output_dir <- file.path(x$base_path, "denoised")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Set cores
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)

  # Determine days to process (mirrors detect_template.Sap)
  if (!is.null(day)) {
    process_metadata <- x$metadata[x$metadata$day_post_hatch %in% day, ]
    days_to_process <- day
    if (nrow(process_metadata) == 0) stop("No files found for specified day(s)")
  } else {
    process_metadata <- x$metadata
    days_to_process <- unique(process_metadata$day_post_hatch)
  }

  # Grand totals for final summary
  total_success <- 0L
  total_skip <- 0L
  total_fail <- 0L

  # ---- Process each day separately ----------------------------------------
  for (current_day in days_to_process) {
    day_metadata <- process_metadata[
      process_metadata$day_post_hatch == current_day,
    ]

    # Apply index filter within this day (if requested)
    if (!is.null(indices)) {
      valid_indices <- indices[indices <= nrow(day_metadata)]
      if (length(valid_indices) == 0) {
        cat(sprintf("\nNo valid indices for day %s. Skipping.\n", current_day))
        next
      }
      day_metadata <- day_metadata[valid_indices, ]
    }

    unique_files <- which(!duplicated(day_metadata$filename))
    n_day <- length(unique_files)

    cat(sprintf(
      "\nProcessing %d files for day %s using %d cores.\n",
      n_day, current_day, cores
    ))

    # Ensure per-day output sub-directory exists
    day_out_dir <- file.path(output_dir, current_day)
    dir.create(day_out_dir, recursive = TRUE, showWarnings = FALSE)

    # Per-file worker for this day — returns "success" / "skipped" / "failed"
    process_file <- function(i) {
      wav_in <- file.path(
        x$base_path,
        day_metadata$day_post_hatch[i],
        day_metadata$filename[i]
      )
      out_path <- file.path(day_out_dir, day_metadata$filename[i])

      if (file.exists(out_path) && !overwrite) {
        return("skipped")
      }

      tryCatch(
        {
          denoise.default(
            x              = wav_in,
            method         = method,
            output_dir     = day_out_dir,
            overwrite      = overwrite,
            wl             = wl,
            ovlp           = ovlp,
            wn             = wn,
            plot           = plot,
            view_window    = view_window,
            freq_range     = freq_range,
            noise_quantile = noise_quantile,
            gain           = gain,
            gate_threshold = gate_threshold,
            gate_smoothing = gate_smoothing,
            gate_floor     = gate_floor,
            verbose        = FALSE
          )
          "success"
        },
        error = function(e) {
          warning(sprintf(
            "Error denoising %s (day %s): %s",
            day_metadata$filename[i], current_day, e$message
          ))
          "failed"
        }
      )
    }

    day_results <- unlist(parallel_apply(unique_files, process_file, cores))

    n_success <- sum(day_results == "success", na.rm = TRUE)
    n_skip <- sum(day_results == "skipped", na.rm = TRUE)
    n_fail <- sum(day_results == "failed", na.rm = TRUE)

    cat(sprintf(
      "Day %s complete: %d denoised, %d skipped (already exists), %d failed.\n",
      current_day, n_success, n_skip, n_fail
    ))

    total_success <- total_success + n_success
    total_skip <- total_skip + n_skip
    total_fail <- total_fail + n_fail
  }

  # ---- Overall summary -----------------------------------------------------
  if (verbose) {
    message(sprintf(
      "\nDenoising complete — %d denoised, %d skipped, %d failed across all days.",
      total_success, total_skip, total_fail
    ))
    message(sprintf("Output directory: %s", normalizePath(output_dir)))
  }

  x$denoised_path <- normalizePath(output_dir)

  if (update_base_path) {
    x$base_path <- normalizePath(output_dir)
    if (verbose) {
      message(sprintf(
        "base_path updated to denoised directory: %s\nAll downstream functions will now read denoised files.",
        x$base_path
      ))
    }
  }

  invisible(x)
}


# Internal helpers ---------------------------------------------------------

# ---- Spectrogram plotting (segment-style) ----

#' @keywords internal
.plot_denoised_spectrogram <- function(samples, sr, wl, ovlp, wn, method,
                                       view_window = NULL) {
  samples_plot <- samples
  time_offset <- 0

  # Segment-style defaults for visualization
  plot_wl <- 256L
  plot_ovlp <- 80L

  if (!is.null(view_window)) {
    if (!is.numeric(view_window) || length(view_window) != 2L) {
      stop("view_window must be a numeric vector of length 2 (seconds).")
    }
    if (any(is.na(view_window)) || any(view_window < 0)) {
      stop("view_window must be non-negative and not NA.")
    }
    if (view_window[1] >= view_window[2]) {
      stop("view_window must be increasing (start < end).")
    }

    duration <- length(samples) / sr
    start_sec <- max(0, view_window[1])
    end_sec <- min(duration, view_window[2])
    if (end_sec <= start_sec) {
      stop("view_window is outside audio duration.")
    }

    start_idx <- max(1L, floor(start_sec * sr) + 1L)
    end_idx <- min(length(samples), ceiling(end_sec * sr))
    samples_plot <- samples[start_idx:end_idx]
    time_offset <- start_sec
  }

  # Compute spectrogram (same backend as segment)
  stft_result <- seewave::spectro(
    samples_plot,
    f = sr,
    wl = plot_wl,
    ovlp = plot_ovlp,
    wn = wn,
    plot = FALSE,
    osc = TRUE,
    cont = TRUE
  )

  times <- stft_result$time + time_offset
  frequencies <- stft_result$freq
  sp <- matrix(stft_result[["amp"]],
    nrow = length(frequencies),
    ncol = length(times)
  )

  # Normalize using segment-style defaults
  spec_norm <- normalize_spec(sp, max_level_db = 40, ref_level_db = 20)
  spec_norm <- spec_norm - matrix(
    rep(apply(spec_norm, 1, median), ncol(spec_norm)),
    nrow = nrow(spec_norm)
  )
  spec_norm[spec_norm < 0] <- 0

  mou_palette <- grDevices::colorRampPalette(c("black", "red", "yellow", "white"))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mar = c(4, 4, 2, 1), bg = "black")
  plot.new()
  plot.window(xlim = range(times), ylim = range(frequencies))

  image(times, frequencies, t(spec_norm),
    col = mou_palette(100),
    add = TRUE,
    useRaster = FALSE,
    interpolate = FALSE
  )

  axis(1, col = "white", col.axis = "white")
  axis(2, col = "white", col.axis = "white")
  title(xlab = "Time (s)", col.lab = "white")
  title(ylab = "Frequency (kHz)", col.lab = "white")
  if (is.null(view_window)) {
    title(main = sprintf("Denoised Spectrogram (%s)", method), col.main = "white")
  } else {
    title(
      main = sprintf(
        "Denoised Spectrogram (%s) [%.2f–%.2fs]",
        method, time_offset, time_offset + max(stft_result$time)
      ),
      col.main = "white"
    )
  }
  box(col = "white")
}

# ---- STFT / iSTFT ----

#' Forward Short-Time Fourier Transform
#' @keywords internal
.forward_stft <- function(samples, wl, hop, win) {
  n_samples <- length(samples)
  n_frames <- floor((n_samples - wl) / hop) + 1L
  n_bins <- floor(wl / 2L) + 1L

  stft <- matrix(complex(real = 0, imaginary = 0),
    nrow = n_bins, ncol = n_frames
  )
  for (t in seq_len(n_frames)) {
    start <- (t - 1L) * hop + 1L
    ft <- stats::fft(samples[start:(start + wl - 1L)] * win)
    stft[, t] <- ft[seq_len(n_bins)]
  }
  stft
}

#' Inverse Short-Time Fourier Transform (overlap-add)
#' @keywords internal
.inverse_stft <- function(stft, wl, hop, win, original_length) {
  n_bins <- nrow(stft)
  n_frames <- ncol(stft)
  out_len <- (n_frames - 1L) * hop + wl
  output <- numeric(out_len)
  win_sum <- numeric(out_len)

  for (t in seq_len(n_frames)) {
    half <- stft[, t]
    full_spec <- if (wl %% 2L == 0L) {
      c(half, Conj(rev(half[-c(1L, n_bins)])))
    } else {
      c(half, Conj(rev(half[-1L])))
    }
    frame <- Re(stats::fft(full_spec, inverse = TRUE)) / wl * win
    idx <- seq.int((t - 1L) * hop + 1L, length.out = wl)
    output[idx] <- output[idx] + frame
    win_sum[idx] <- win_sum[idx] + win^2
  }

  nz <- win_sum > 1e-8
  output[nz] <- output[nz] / win_sum[nz]
  if (length(output) > original_length) output <- output[seq_len(original_length)]
  peak <- max(abs(output))
  if (peak > 0) output <- output / peak
  output
}


# ---- Method 1 : Spectral Median Subtraction ----

#' @keywords internal
.spectral_median_subtract <- function(mag,
                                      noise_quantile = 0.25,
                                      gain = 1.0) {
  noise_floor <- apply(mag, 1L, quantile, probs = noise_quantile, na.rm = TRUE)
  clean_mag <- mag - gain * noise_floor # broadcast along columns
  clean_mag[clean_mag < 0] <- 0
  clean_mag
}


# ---- Method 2 : Spectral Gating ----

#' Spectral Gating (soft sigmoid gate with configurable floor)
#'
#' @description
#' Builds a per-frequency noise profile (floor \eqn{\mu_f} and spread
#' \eqn{\sigma_f}) then computes a gate value for every time-frequency bin:
#' \deqn{G(t,f) = \max\!\left(g_{\min},\;
#'   \frac{1}{1 + e^{-(|S(t,f)| - \mu_f - \theta\sigma_f)\,/\,(\sigma_f/4)}}
#' \right)}
#' The floor \eqn{g_{\min}} (\code{floor}) prevents complete silencing and
#' is the key parameter for preserving tonal continuity of bird song.
#'
#' @param mag      Numeric matrix (freq x time) — magnitude spectrogram.
#' @param noise_quantile Quantile for noise-floor estimation.
#' @param threshold Number of \eqn{\sigma} above noise floor where gate reaches
#'   50 \% transmission (default: 1.5).
#' @param smoothing Half-width (bins) of frequency-axis box-car mask smoother;
#'   0 = disabled (default: 3).
#' @param floor     Minimum gate value in [0, 1).
#' @return Cleaned magnitude matrix.
#' @keywords internal
.spectral_gate <- function(mag,
                           noise_quantile = 0.25,
                           threshold = 1.5,
                           smoothing = 3L,
                           floor = 0.1) {
  # --- Step 1: Identify quiet (background-only) frames ----------------------
  # Rank frames by total power; take the lowest noise_quantile fraction.
  # Song syllables are brief and high-energy so they are excluded here,
  # giving a clean noise-only sample for profile estimation.
  frame_energy <- colSums(mag^2)
  energy_threshold <- quantile(frame_energy,
    probs = noise_quantile,
    na.rm = TRUE
  )
  quiet_frames <- mag[, frame_energy <= energy_threshold, drop = FALSE]

  if (ncol(quiet_frames) < 2L) {
    warning(
      "spectral_gate: fewer than 2 quiet frames found; ",
      "consider raising noise_quantile."
    )
    quiet_frames <- mag # fall back to all frames
  }

  # --- Step 2: Per-frequency noise profile from quiet frames only -----------
  noise_floor <- rowMeans(quiet_frames, na.rm = TRUE)
  noise_sd <- apply(quiet_frames, 1L, stats::sd, na.rm = TRUE)
  noise_sd <- pmax(noise_sd, 1e-8) # guard against zero-variance bins

  # --- Step 3: Soft sigmoid gate --------------------------------------------
  gate_centre <- noise_floor + threshold * noise_sd
  steepness <- noise_sd / 4 # ~1 SD transition width

  z <- sweep(mag, 1L, gate_centre, FUN = "-")
  z <- sweep(z, 1L, steepness, FUN = "/")
  gate <- 1 / (1 + exp(-z))

  # Apply floor — no bin is ever fully silenced
  gate <- pmax(gate, floor)

  # --- Step 4: Optional frequency-axis smoothing of the gate mask -----------
  if (smoothing > 0L) {
    hw <- as.integer(smoothing)
    kern <- rep(1, 2L * hw + 1L) / (2L * hw + 1L)
    gate <- apply(gate, 2L, function(col) {
      out <- stats::filter(col, kern, circular = FALSE, sides = 2)
      out[is.na(out)] <- col[is.na(out)] # restore boundary NAs
      pmax(out, floor)
    })
  }

  mag * gate
}
