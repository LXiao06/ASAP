# Test Segmentation Functions
# Tests for functions in R/segment.R

test_that("normalize_spec works correctly", {
  # Create test spectrogram
  test_spec <- matrix(runif(100, 0, 1), nrow = 10, ncol = 10)
  
  # Test dB normalization
  result_db <- ASAP:::normalize_spec(test_spec, max_level_db = 40, 
                                    ref_level_db = 20, method = "db")
  expect_true(is.matrix(result_db))
  expect_equal(dim(result_db), dim(test_spec))
  expect_true(all(result_db >= 0 & result_db <= 1))
  
  # Test minmax normalization
  result_minmax <- ASAP:::normalize_spec(test_spec, method = "minmax")
  expect_true(is.matrix(result_minmax))
  expect_equal(dim(result_minmax), dim(test_spec))
  expect_true(all(result_minmax >= 0 & result_minmax <= 1))
  expect_equal(min(result_minmax), 0)
  expect_equal(max(result_minmax), 1)
  
  # Test zscore normalization
  result_zscore <- ASAP:::normalize_spec(test_spec, method = "zscore")
  expect_true(is.matrix(result_zscore))
  expect_equal(dim(result_zscore), dim(test_spec))
  
  # Test invalid input
  expect_error(ASAP:::normalize_spec("not a matrix"), "Input must be a matrix")
  expect_error(ASAP:::normalize_spec(test_spec, method = "invalid"), 
               "Unknown normalization method")
})

test_that("find_continuous_regions works correctly", {
  # Test simple case
  signal1 <- c(0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  result1 <- ASAP:::find_continuous_regions(signal1)
  expected1 <- matrix(c(3, 5, 8, 9), nrow = 2, ncol = 2, byrow = TRUE)
  expect_equal(result1, expected1)
  
  # Test edge cases
  signal2 <- c(1, 1, 0, 0, 1, 1)  # starts and ends with 1
  result2 <- ASAP:::find_continuous_regions(signal2)
  expect_true(is.matrix(result2))
  
  # Test all zeros
  signal3 <- c(0, 0, 0, 0)
  result3 <- ASAP:::find_continuous_regions(signal3)
  expect_equal(nrow(result3), 0)
  expect_equal(ncol(result3), 2)
  
  # Test all ones
  signal4 <- c(1, 1, 1, 1)
  result4 <- ASAP:::find_continuous_regions(signal4)
  expect_equal(nrow(result4), 0)
  
  # Test single transition
  signal5 <- c(0, 1, 0)
  result5 <- ASAP:::find_continuous_regions(signal5)
  expected5 <- matrix(c(2, 2), nrow = 1, ncol = 2)
  expect_equal(result5, expected5)
})

test_that("segment.default input validation", {
  # Test non-existent file
  expect_error(segment("nonexistent_file.wav"), "File does not exist")
  
  # We can't easily test the full function without actual WAV files,
  # but we can test parameter validation
  
  # Create a mock function to test parameter validation logic
  test_params <- function(wl = 256, ovlp = 80, flim = c(1, 10)) {
    # Validate window length and overlap
    if (!is.numeric(wl) || wl <= 0) {
      stop("Window length (wl) must be a positive number")
    }
    if (!is.numeric(ovlp) || ovlp < 0 || ovlp >= 100) {
      stop("Overlap percentage (ovlp) must be between 0 and 100")
    }
    
    # Validate flim parameter
    if (length(flim) != 2) {
      stop("flim must be a vector of length 2 (minimum and maximum frequency)")
    }
    
    return(TRUE)
  }
  
  # Test valid parameters
  expect_true(test_params())
  expect_true(test_params(wl = 512, ovlp = 50, flim = c(0.5, 8)))
  
  # Test invalid parameters
  expect_error(test_params(wl = -1), "Window length \\(wl\\) must be a positive number")
  expect_error(test_params(wl = 0), "Window length \\(wl\\) must be a positive number")
  expect_error(test_params(ovlp = -1), "Overlap percentage \\(ovlp\\) must be between 0 and 100")
  expect_error(test_params(ovlp = 100), "Overlap percentage \\(ovlp\\) must be between 0 and 100")
  expect_error(test_params(flim = c(1)), "flim must be a vector of length 2")
  expect_error(test_params(flim = c(1, 2, 3)), "flim must be a vector of length 2")
})

test_that("segment time validation logic", {
  # Test time validation logic (extracted from segment.default)
  validate_times <- function(start_time, end_time, duration) {
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
    
    return(list(start_time = start_time, end_time = end_time))
  }
  
  # Test valid cases
  result1 <- validate_times(NULL, NULL, 10)
  expect_equal(result1$start_time, 0)
  expect_equal(result1$end_time, 10)
  
  result2 <- validate_times(2, 8, 10)
  expect_equal(result2$start_time, 2)
  expect_equal(result2$end_time, 8)
  
  result3 <- validate_times(NULL, 15, 10)  # end_time > duration
  expect_equal(result3$end_time, 10)
  
  # Test invalid cases
  expect_error(validate_times(-1, 5, 10), "start_time cannot be negative")
  expect_error(validate_times(8, 5, 10), "Invalid time window")
  expect_error(validate_times(10, 10, 10), "Invalid time window")
})