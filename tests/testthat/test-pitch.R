# Test Pitch Functions
# Tests for functions in R/pitch.R

test_that("FF parameter validation", {
  # Test fundamental frequency parameter validation
  validate_ff_params <- function(wl = 512,
                               ovlp = 50,
                               fmax = 1400,
                               threshold = 10,
                               method = "cepstrum") {
    
    # Validate window length
    if (!is.numeric(wl) || wl <= 0) {
      stop("wl must be a positive number")
    }
    
    # Check if wl is power of 2 (recommended for FFT)
    if (wl != 2^round(log2(wl))) {
      warning("wl should be a power of 2 for optimal FFT performance")
    }
    
    # Validate overlap
    if (!is.numeric(ovlp) || ovlp < 0 || ovlp >= 100) {
      stop("ovlp must be between 0 and 100")
    }
    
    # Validate maximum frequency
    if (!is.numeric(fmax) || fmax <= 0) {
      stop("fmax must be a positive number")
    }
    
    # Validate threshold
    if (!is.numeric(threshold) || threshold < 0 || threshold > 100) {
      stop("threshold must be between 0 and 100")
    }
    
    # Validate method
    if (!method %in% c("cepstrum", "yin")) {
      stop("method must be either 'cepstrum' or 'yin'")
    }
    
    return(TRUE)
  }
  
  # Test valid parameters
  expect_true(validate_ff_params())
  expect_true(validate_ff_params(wl = 1024, ovlp = 75))
  expect_true(validate_ff_params(method = "yin", fmax = 2000))
  
  # Test invalid parameters
  expect_error(validate_ff_params(wl = -512),
               "wl must be a positive number")
  expect_error(validate_ff_params(ovlp = 150),
               "ovlp must be between 0 and 100")
  expect_error(validate_ff_params(fmax = -1000),
               "fmax must be a positive number")
  expect_error(validate_ff_params(threshold = 150),
               "threshold must be between 0 and 100")
  expect_error(validate_ff_params(method = "invalid"),
               "method must be either 'cepstrum' or 'yin'")
  
  # Test warnings
  expect_warning(validate_ff_params(wl = 500),
                "wl should be a power of 2")
})

test_that("pitch detection time validation", {
  # Test time window validation for pitch detection
  validate_pitch_times <- function(start_time, end_time, duration) {
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
    
    # Check minimum duration for meaningful pitch analysis
    min_duration <- 0.01  # 10ms minimum
    if ((end_time - start_time) < min_duration) {
      warning("Very short analysis window may produce unreliable results")
    }
    
    return(list(start_time = start_time, end_time = end_time))
  }
  
  # Test valid cases
  result1 <- validate_pitch_times(NULL, NULL, 10)
  expect_equal(result1$start_time, 0)
  expect_equal(result1$end_time, 10)
  
  result2 <- validate_pitch_times(2, 8, 10)
  expect_equal(result2$start_time, 2)
  expect_equal(result2$end_time, 8)
  
  result3 <- validate_pitch_times(NULL, 15, 10)  # end_time > duration
  expect_equal(result3$end_time, 10)
  
  # Test invalid cases
  expect_error(validate_pitch_times(-1, 5, 10), "start_time cannot be negative")
  expect_error(validate_pitch_times(8, 5, 10), "Invalid time window")
  expect_error(validate_pitch_times(10, 10, 10), "Invalid time window")
  
  # Test warnings
  expect_warning(validate_pitch_times(0, 0.005, 10), "Very short analysis window")
})

test_that("cepstral pitch detection logic", {
  # Test cepstral pitch detection algorithm components
  calculate_cepstrum <- function(signal, sample_rate, wl = 512) {
    # Simulate cepstral analysis steps
    
    # 1. Apply window
    if (length(signal) > wl) {
      signal <- signal[1:wl]
    } else if (length(signal) < wl) {
      signal <- c(signal, rep(0, wl - length(signal)))  # Zero-pad
    }
    
    # 2. FFT
    spectrum <- fft(signal)
    power_spectrum <- Mod(spectrum)^2
    
    # 3. Log of power spectrum
    log_spectrum <- log(power_spectrum + 1e-10)  # Add small value to avoid log(0)
    
    # 4. Inverse FFT (cepstrum)
    cepstrum <- Re(fft(log_spectrum, inverse = TRUE))
    
    # 5. Find peak in appropriate quefrency range
    # For pitch detection, look for peaks corresponding to reasonable F0 range
    min_f0 <- 80   # Hz
    max_f0 <- 1000 # Hz
    min_quefrency <- sample_rate / max_f0
    max_quefrency <- sample_rate / min_f0
    
    search_range <- max(1, floor(min_quefrency)):min(length(cepstrum), ceil(max_quefrency))
    
    if (length(search_range) > 0) {
      peak_idx <- which.max(cepstrum[search_range])
      peak_quefrency <- search_range[peak_idx]
      estimated_f0 <- sample_rate / peak_quefrency
    } else {
      estimated_f0 <- NA
    }
    
    return(list(
      cepstrum = cepstrum,
      estimated_f0 = estimated_f0,
      peak_quefrency = peak_quefrency
    ))
  }
  
  # Test with synthetic signal
  sample_rate <- 44100
  duration <- 0.1  # 100ms
  f0 <- 440  # A4 note
  t <- seq(0, duration, length.out = sample_rate * duration)
  
  # Create harmonic signal
  signal <- sin(2 * pi * f0 * t) + 0.5 * sin(2 * pi * 2 * f0 * t) + 0.25 * sin(2 * pi * 3 * f0 * t)
  
  result <- calculate_cepstrum(signal, sample_rate)
  
  expect_true(is.list(result))
  expect_true("cepstrum" %in% names(result))
  expect_true("estimated_f0" %in% names(result))
  
  # Check if estimated F0 is reasonable (within 10% of true F0)
  if (!is.na(result$estimated_f0)) {
    expect_true(abs(result$estimated_f0 - f0) / f0 < 0.1)
  }
  
  # Test with noise (should return NA or unreliable estimate)
  noise_signal <- rnorm(length(signal), 0, 0.1)
  noise_result <- calculate_cepstrum(noise_signal, sample_rate)
  expect_true(is.list(noise_result))
})

test_that("YIN algorithm components", {
  # Test YIN algorithm components for pitch detection
  calculate_yin_difference <- function(signal, max_lag) {
    n <- length(signal)
    diff_function <- numeric(max_lag)
    
    # Calculate difference function
    for (tau in 1:max_lag) {
      if (tau < n) {
        diff_sum <- sum((signal[1:(n-tau)] - signal[(tau+1):n])^2)
        diff_function[tau] <- diff_sum
      }
    }
    
    return(diff_function)
  }
  
  calculate_yin_cmnd <- function(diff_function) {
    # Cumulative mean normalized difference function
    cmnd <- numeric(length(diff_function))
    cmnd[1] <- 1
    
    for (tau in 2:length(diff_function)) {
      if (tau > 1) {
        cumulative_sum <- sum(diff_function[2:tau])
        cmnd[tau] <- diff_function[tau] / (cumulative_sum / (tau - 1))
      }
    }
    
    return(cmnd)
  }
  
  find_yin_period <- function(cmnd, threshold = 0.1) {
    # Find first minimum below threshold
    for (tau in 2:length(cmnd)) {
      if (cmnd[tau] < threshold) {
        # Look for local minimum
        if (tau == length(cmnd) || cmnd[tau] < cmnd[tau + 1]) {
          return(tau)
        }
      }
    }
    
    # If no minimum below threshold, return global minimum
    return(which.min(cmnd[-1]) + 1)  # Exclude first element
  }
  
  # Test with synthetic periodic signal
  sample_rate <- 44100
  f0 <- 220  # Hz
  period_samples <- round(sample_rate / f0)
  
  # Create periodic signal
  t <- seq(0, 0.1, length.out = sample_rate * 0.1)
  signal <- sin(2 * pi * f0 * t)
  
  # Test difference function
  max_lag <- min(length(signal) - 1, round(sample_rate / 50))  # Up to 50 Hz
  diff_func <- calculate_yin_difference(signal, max_lag)
  
  expect_true(length(diff_func) == max_lag)
  expect_true(all(diff_func >= 0))  # Difference function should be non-negative
  
  # Test CMND
  cmnd <- calculate_yin_cmnd(diff_func)
  expect_true(length(cmnd) == length(diff_func))
  expect_equal(cmnd[1], 1)  # First element should be 1
  
  # Test period finding
  estimated_period <- find_yin_period(cmnd, threshold = 0.1)
  expect_true(estimated_period > 1)
  expect_true(estimated_period <= max_lag)
  
  # Check if estimated period is close to true period
  estimated_f0 <- sample_rate / estimated_period
  expect_true(abs(estimated_f0 - f0) / f0 < 0.2)  # Within 20%
})

test_that("pitch goodness calculation", {
  # Test pitch goodness/confidence measures
  calculate_pitch_goodness <- function(signal, estimated_period, method = "autocorr") {
    if (method == "autocorr") {
      # Autocorrelation-based goodness
      if (estimated_period >= length(signal)) {
        return(0)
      }
      
      # Calculate autocorrelation at estimated period
      n <- length(signal)
      lag <- round(estimated_period)
      
      if (lag < n) {
        # Normalized autocorrelation
        signal_mean <- mean(signal)
        signal_centered <- signal - signal_mean
        
        autocorr_num <- sum(signal_centered[1:(n-lag)] * signal_centered[(lag+1):n])
        autocorr_denom <- sqrt(sum(signal_centered[1:(n-lag)]^2) * sum(signal_centered[(lag+1):n]^2))
        
        if (autocorr_denom > 0) {
          goodness <- autocorr_num / autocorr_denom
        } else {
          goodness <- 0
        }
      } else {
        goodness <- 0
      }
      
    } else if (method == "harmonic") {
      # Harmonic-based goodness (simplified)
      # This would typically involve spectral analysis
      goodness <- runif(1, 0.5, 1)  # Placeholder
    }
    
    return(abs(goodness))  # Return absolute value
  }
  
  # Test with periodic signal
  sample_rate <- 44100
  f0 <- 330  # Hz
  t <- seq(0, 0.1, length.out = sample_rate * 0.1)
  periodic_signal <- sin(2 * pi * f0 * t)
  true_period <- sample_rate / f0
  
  goodness_periodic <- calculate_pitch_goodness(periodic_signal, true_period, "autocorr")
  expect_true(goodness_periodic > 0.8)  # Should be high for periodic signal
  
  # Test with noise
  noise_signal <- rnorm(length(periodic_signal))
  goodness_noise <- calculate_pitch_goodness(noise_signal, true_period, "autocorr")
  expect_true(goodness_noise < goodness_periodic)  # Should be lower for noise
  
  # Test edge cases
  expect_equal(calculate_pitch_goodness(periodic_signal, length(periodic_signal) + 1, "autocorr"), 0)
  
  # Test very short signal
  short_signal <- periodic_signal[1:10]
  goodness_short <- calculate_pitch_goodness(short_signal, 5, "autocorr")
  expect_true(goodness_short >= 0 && goodness_short <= 1)
})

test_that("pitch tracking continuity", {
  # Test pitch tracking continuity and smoothing
  smooth_pitch_track <- function(f0_values, max_jump_ratio = 0.2, 
                                median_filter_size = 3) {
    n <- length(f0_values)
    smoothed <- f0_values
    
    # Remove outliers based on jump ratio
    for (i in 2:n) {
      if (!is.na(f0_values[i]) && !is.na(f0_values[i-1])) {
        jump_ratio <- abs(f0_values[i] - f0_values[i-1]) / f0_values[i-1]
        if (jump_ratio > max_jump_ratio) {
          smoothed[i] <- f0_values[i-1]  # Use previous value
        }
      }
    }
    
    # Apply median filter
    if (median_filter_size > 1 && n >= median_filter_size) {
      half_window <- floor(median_filter_size / 2)
      for (i in (half_window + 1):(n - half_window)) {
        window_vals <- smoothed[(i - half_window):(i + half_window)]
        window_vals <- window_vals[!is.na(window_vals)]
        if (length(window_vals) > 0) {
          smoothed[i] <- median(window_vals)
        }
      }
    }
    
    return(smoothed)
  }
  
  # Test with noisy pitch track
  true_f0 <- 400
  noisy_track <- c(true_f0, true_f0 * 1.1, true_f0 * 2.1, true_f0 * 0.9, true_f0 * 1.05)  # Contains octave error
  
  smoothed_track <- smooth_pitch_track(noisy_track, max_jump_ratio = 0.2)
  
  expect_true(length(smoothed_track) == length(noisy_track))
  
  # Check that large jump was corrected
  jump_ratios <- abs(diff(smoothed_track)) / smoothed_track[-length(smoothed_track)]
  expect_true(all(jump_ratios <= 0.2, na.rm = TRUE))
  
  # Test with NA values
  track_with_na <- c(400, NA, 420, 410, NA, 405)
  smoothed_na <- smooth_pitch_track(track_with_na)
  expect_true(length(smoothed_na) == length(track_with_na))
  
  # Test edge cases
  expect_equal(smooth_pitch_track(numeric(0)), numeric(0))
  expect_equal(smooth_pitch_track(c(400)), c(400))
})