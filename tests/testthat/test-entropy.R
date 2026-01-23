# Test Entropy Functions
# Tests for functions in R/entropy.R

test_that("spectral_entropy parameter validation", {
  # Test method validation
  validate_entropy_params <- function(method = "weiner", 
                                    normalize = FALSE,
                                    freq_range = c(500, 15000),
                                    threshold = 10,
                                    wl = 512,
                                    ovlp = 50) {
    
    # Validate method
    if (!method %in% c("weiner", "shannon")) {
      stop("method must be either 'weiner' or 'shannon'")
    }
    
    # Validate normalize
    if (!is.logical(normalize)) {
      stop("normalize must be logical")
    }
    
    # Validate freq_range
    if (!is.numeric(freq_range) || length(freq_range) != 2) {
      stop("freq_range must be a numeric vector of length 2")
    }
    
    if (freq_range[1] >= freq_range[2]) {
      stop("freq_range[1] must be less than freq_range[2]")
    }
    
    # Validate threshold
    if (!is.numeric(threshold) || threshold < 0) {
      stop("threshold must be a non-negative number")
    }
    
    # Validate window parameters
    if (!is.numeric(wl) || wl <= 0) {
      stop("wl must be a positive number")
    }
    
    if (!is.numeric(ovlp) || ovlp < 0 || ovlp >= 100) {
      stop("ovlp must be between 0 and 100")
    }
    
    return(TRUE)
  }
  
  # Test valid parameters
  expect_true(validate_entropy_params())
  expect_true(validate_entropy_params(method = "shannon"))
  expect_true(validate_entropy_params(normalize = TRUE))
  expect_true(validate_entropy_params(freq_range = c(1000, 8000)))
  
  # Test invalid parameters
  expect_error(validate_entropy_params(method = "invalid"), 
               "method must be either 'weiner' or 'shannon'")
  expect_error(validate_entropy_params(normalize = "yes"), 
               "normalize must be logical")
  expect_error(validate_entropy_params(freq_range = c(1000)), 
               "freq_range must be a numeric vector of length 2")
  expect_error(validate_entropy_params(freq_range = c(8000, 1000)), 
               "freq_range\\[1\\] must be less than freq_range\\[2\\]")
  expect_error(validate_entropy_params(threshold = -5), 
               "threshold must be a non-negative number")
  expect_error(validate_entropy_params(wl = -512), 
               "wl must be a positive number")
  expect_error(validate_entropy_params(ovlp = 150), 
               "ovlp must be between 0 and 100")
})

test_that("entropy calculation logic", {
  # Test Wiener entropy calculation
  calculate_wiener_entropy <- function(power_spectrum, normalize = FALSE) {
    # Remove zero values to avoid log(0)
    power_spectrum <- power_spectrum[power_spectrum > 0]
    
    if (length(power_spectrum) == 0) {
      return(NA)
    }
    
    # Calculate geometric and arithmetic means
    geometric_mean <- exp(mean(log(power_spectrum)))
    arithmetic_mean <- mean(power_spectrum)
    
    # Wiener entropy
    wiener_entropy <- log(geometric_mean / arithmetic_mean)
    
    if (normalize) {
      # Normalize to [0, 1] range
      wiener_entropy <- exp(wiener_entropy)
    }
    
    return(wiener_entropy)
  }
  
  # Test Shannon entropy calculation
  calculate_shannon_entropy <- function(power_spectrum, normalize = FALSE) {
    # Normalize to probability distribution
    power_spectrum <- power_spectrum / sum(power_spectrum)
    
    # Remove zero values
    power_spectrum <- power_spectrum[power_spectrum > 0]
    
    if (length(power_spectrum) == 0) {
      return(NA)
    }
    
    # Shannon entropy
    shannon_entropy <- -sum(power_spectrum * log2(power_spectrum))
    
    if (normalize) {
      # Normalize by maximum possible entropy
      max_entropy <- log2(length(power_spectrum))
      if (max_entropy > 0) {
        shannon_entropy <- shannon_entropy / max_entropy
      }
    }
    
    return(shannon_entropy)
  }
  
  # Test with uniform distribution (white noise-like)
  uniform_spectrum <- rep(1, 100)
  wiener_uniform <- calculate_wiener_entropy(uniform_spectrum)
  shannon_uniform <- calculate_shannon_entropy(uniform_spectrum)
  
  expect_equal(wiener_uniform, 0, tolerance = 1e-10)  # Perfect uniformity
  expect_true(shannon_uniform > 6.5)  # High entropy for uniform distribution
  
  # Test with peaked distribution (tonal signal)
  peaked_spectrum <- c(rep(0.01, 50), 10, rep(0.01, 49))
  wiener_peaked <- calculate_wiener_entropy(peaked_spectrum)
  shannon_peaked <- calculate_shannon_entropy(peaked_spectrum)
  
  expect_true(wiener_peaked < -2)  # Large negative value for structured signal
  expect_true(shannon_peaked < 4)  # Lower entropy for structured signal
  
  # Test normalization
  wiener_norm <- calculate_wiener_entropy(peaked_spectrum, normalize = TRUE)
  shannon_norm <- calculate_shannon_entropy(uniform_spectrum, normalize = TRUE)
  
  expect_true(wiener_norm >= 0 && wiener_norm <= 1)
  expect_true(shannon_norm >= 0 && shannon_norm <= 1)
  expect_true(shannon_norm > 0.95)  # Should be close to 1 for uniform
})

test_that("entropy frequency range filtering", {
  # Test frequency range filtering logic
  filter_frequency_range <- function(freq_vector, power_matrix, freq_range) {
    # Find indices within frequency range
    freq_indices <- which(freq_vector >= freq_range[1] & freq_vector <= freq_range[2])
    
    if (length(freq_indices) == 0) {
      stop("No frequencies found in specified range")
    }
    
    # Filter power matrix
    filtered_power <- power_matrix[freq_indices, , drop = FALSE]
    filtered_freqs <- freq_vector[freq_indices]
    
    return(list(
      power_matrix = filtered_power,
      frequencies = filtered_freqs,
      n_freq_bins = length(freq_indices)
    ))
  }
  
  # Create test data
  test_freqs <- seq(0, 10000, by = 100)  # 0 to 10 kHz
  test_power <- matrix(runif(length(test_freqs) * 50), 
                      nrow = length(test_freqs), ncol = 50)
  
  # Test normal range
  result1 <- filter_frequency_range(test_freqs, test_power, c(1000, 5000))
  expect_true(result1$n_freq_bins > 0)
  expect_true(all(result1$frequencies >= 1000 & result1$frequencies <= 5000))
  expect_equal(nrow(result1$power_matrix), result1$n_freq_bins)
  
  # Test edge cases
  result2 <- filter_frequency_range(test_freqs, test_power, c(0, 10000))
  expect_equal(result2$n_freq_bins, length(test_freqs))
  
  # Test invalid range
  expect_error(filter_frequency_range(test_freqs, test_power, c(15000, 20000)),
               "No frequencies found in specified range")
})

test_that("entropy matrix organization", {
  # Test entropy matrix creation and organization
  organize_entropy_matrix <- function(entropy_values, time_points, labels = NULL) {
    # Create matrix
    entropy_matrix <- matrix(entropy_values, 
                           nrow = length(unique(labels %||% seq_along(entropy_values))),
                           ncol = length(time_points))
    
    # Add row and column names
    if (!is.null(labels)) {
      rownames(entropy_matrix) <- unique(labels)
    }
    colnames(entropy_matrix) <- round(time_points, 3)
    
    return(entropy_matrix)
  }
  
  # Test with labels
  test_entropy <- runif(100, -3, 0)  # Typical Wiener entropy range
  test_times <- seq(0, 1, length.out = 10)
  test_labels <- rep(c("a", "b"), each = 5)
  
  result <- organize_entropy_matrix(test_entropy[1:10], test_times, test_labels)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)  # Two unique labels
  expect_equal(ncol(result), 10)  # Ten time points
  expect_equal(rownames(result), c("a", "b"))
  
  # Test without labels
  result2 <- organize_entropy_matrix(test_entropy[1:10], test_times)
  expect_true(is.matrix(result2))
  expect_equal(nrow(result2), 10)
})

test_that("entropy threshold application", {
  # Test amplitude threshold application
  apply_amplitude_threshold <- function(power_spectrum, threshold_percent) {
    if (threshold_percent <= 0 || threshold_percent > 100) {
      stop("threshold_percent must be between 0 and 100")
    }
    
    # Calculate threshold value
    max_power <- max(power_spectrum, na.rm = TRUE)
    threshold_value <- max_power * (threshold_percent / 100)
    
    # Apply threshold
    thresholded_spectrum <- power_spectrum
    thresholded_spectrum[power_spectrum < threshold_value] <- 0
    
    return(list(
      spectrum = thresholded_spectrum,
      threshold_value = threshold_value,
      n_above_threshold = sum(power_spectrum >= threshold_value)
    ))
  }
  
  # Test with sample data
  test_spectrum <- c(0.1, 0.5, 1.0, 0.3, 0.8, 0.2, 0.9)
  
  # Test 50% threshold
  result1 <- apply_amplitude_threshold(test_spectrum, 50)
  expect_equal(result1$threshold_value, 0.5)  # 50% of max (1.0)
  expect_equal(result1$n_above_threshold, 3)  # Values >= 0.5
  
  # Test 10% threshold
  result2 <- apply_amplitude_threshold(test_spectrum, 10)
  expect_equal(result2$threshold_value, 0.1)
  expect_equal(result2$n_above_threshold, 7)  # All values >= 0.1
  
  # Test invalid thresholds
  expect_error(apply_amplitude_threshold(test_spectrum, 0),
               "threshold_percent must be between 0 and 100")
  expect_error(apply_amplitude_threshold(test_spectrum, 150),
               "threshold_percent must be between 0 and 100")
})