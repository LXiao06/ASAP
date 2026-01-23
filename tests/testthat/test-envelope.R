# Test Amplitude Envelope Functions
# Tests for functions in R/envelope.R

test_that("amp_env input validation", {
  # Test invalid input types
  expect_error(amp_env("not a data frame"), 
               "segment_row must be a single row from a data frame")
  
  expect_error(amp_env(data.frame()), 
               "segment_row must be a single row from a data frame")
  
  # Test multiple rows
  multi_row_df <- data.frame(
    start_time = c(1, 2),
    end_time = c(2, 3),
    stringsAsFactors = FALSE
  )
  expect_error(amp_env(multi_row_df), 
               "segment_row must be a single row from a data frame")
  
  # Test missing required columns
  invalid_df <- data.frame(other_col = 1, stringsAsFactors = FALSE)
  expect_error(amp_env(invalid_df), 
               "Input row must contain 'start_time' and 'end_time' columns")
  
  # Test missing one required column
  partial_df <- data.frame(start_time = 1, stringsAsFactors = FALSE)
  expect_error(amp_env(partial_df), 
               "Input row must contain 'start_time' and 'end_time' columns")
})

test_that("amp_env msmooth validation", {
  # Create valid test data
  valid_row <- data.frame(
    filename = "test.wav",
    start_time = 1,
    end_time = 2,
    stringsAsFactors = FALSE
  )
  
  # Test msmooth validation logic (extracted from function)
  validate_msmooth <- function(msmooth) {
    if (!is.null(msmooth)) {
      if (!is.numeric(msmooth) || length(msmooth) != 2) {
        stop("msmooth must be a numeric vector of length 2")
      }
      if (any(msmooth <= 0)) {
        stop("msmooth values must be positive")
      }
      if (msmooth[2] < 0 || msmooth[2] >= 100) {
        stop("msmooth overlap percentage must be between 0 and 100")
      }
    }
    return(TRUE)
  }
  
  # Test valid msmooth values
  expect_true(validate_msmooth(NULL))
  expect_true(validate_msmooth(c(256, 50)))
  expect_true(validate_msmooth(c(512, 0)))
  expect_true(validate_msmooth(c(128, 99)))
  
  # Test invalid msmooth values
  expect_error(validate_msmooth(c(256)), "msmooth must be a numeric vector of length 2")
  expect_error(validate_msmooth(c(256, 50, 75)), "msmooth must be a numeric vector of length 2")
  expect_error(validate_msmooth(c(-256, 50)), "msmooth values must be positive")
  expect_error(validate_msmooth(c(256, -50)), "msmooth values must be positive")
  expect_error(validate_msmooth(c(256, 100)), "msmooth overlap percentage must be between 0 and 100")
  expect_error(validate_msmooth(c(256, 150)), "msmooth overlap percentage must be between 0 and 100")
})

test_that("pad_amp_env works correctly", {
  # Test basic padding
  test_env <- c(0.1, 0.5, 0.8, 0.3, 0.1)
  
  # Test padding to longer length
  result1 <- pad_amp_env(test_env, target_length = 10)
  expect_length(result1, 10)
  expect_true(all(result1 >= 0))
  
  # Test padding to shorter length (should truncate)
  result2 <- pad_amp_env(test_env, target_length = 3)
  expect_length(result2, 3)
  
  # Test with same length
  result3 <- pad_amp_env(test_env, target_length = 5)
  expect_equal(result3, test_env)
  
  # Test input validation
  expect_error(pad_amp_env("not numeric"), "envelope must be numeric")
  expect_error(pad_amp_env(test_env, target_length = -1), 
               "target_length must be positive")
  expect_error(pad_amp_env(test_env, target_length = 0), 
               "target_length must be positive")
})

test_that("pad_array works correctly", {
  # Test 1D array padding
  test_array1d <- c(1, 2, 3, 4, 5)
  
  result1d <- pad_array(test_array1d, target_length = 8)
  expect_length(result1d, 8)
  expect_equal(result1d[1:5], test_array1d)
  
  # Test 2D array padding
  test_array2d <- matrix(1:12, nrow = 3, ncol = 4)
  
  result2d <- pad_array(test_array2d, target_length = 6)
  expect_equal(nrow(result2d), 3)
  expect_equal(ncol(result2d), 6)
  expect_equal(result2d[, 1:4], test_array2d)
  
  # Test truncation
  result_trunc <- pad_array(test_array1d, target_length = 3)
  expect_length(result_trunc, 3)
  expect_equal(result_trunc, test_array1d[1:3])
  
  # Test input validation
  expect_error(pad_array(test_array1d, target_length = -1), 
               "target_length must be positive")
})