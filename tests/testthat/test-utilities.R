# Test Utility Functions
# Tests for functions in R/utillties.R

test_that("construct_wav_path works correctly", {
  # Test with wav_dir argument
  test_row <- data.frame(
    filename = "test.wav",
    day_post_hatch = "day1",
    stringsAsFactors = FALSE
  )
  
  # Create temporary directory structure for testing
  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_audio")
  day_dir <- file.path(test_dir, "day1")
  dir.create(day_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create a dummy wav file
  test_file <- file.path(day_dir, "test.wav")
  file.create(test_file)
  
  # Test path construction
  result_path <- ASAP:::construct_wav_path(test_row, wav_dir = test_dir)
  expect_equal(result_path, test_file)
  
  # Test with wav_dir as attribute
  attr(test_row, "wav_dir") <- test_dir
  result_path2 <- ASAP:::construct_wav_path(test_row)
  expect_equal(result_path2, test_file)
  
  # Clean up
  unlink(test_dir, recursive = TRUE)
})

test_that("construct_wav_path error handling", {
  # Test missing filename column
  test_row <- data.frame(other_col = "value", stringsAsFactors = FALSE)
  expect_error(ASAP:::construct_wav_path(test_row, wav_dir = "dummy"),
               "Input must contain 'filename' column")
  
  # Test missing wav_dir
  test_row <- data.frame(filename = "test.wav", stringsAsFactors = FALSE)
  expect_error(ASAP:::construct_wav_path(test_row),
               "wav_dir must be provided either as argument or attribute")
  
  # Test non-existent file
  test_row <- data.frame(filename = "nonexistent.wav", stringsAsFactors = FALSE)
  expect_error(ASAP:::construct_wav_path(test_row, wav_dir = tempdir()),
               "Sound file not found")
})

test_that("parallel_apply works correctly", {
  # Test basic functionality
  test_indices <- 1:5
  test_function <- function(x) x^2
  
  result <- ASAP:::parallel_apply(test_indices, test_function, cores = 1)
  expected <- list(1, 4, 9, 16, 25)
  expect_equal(result, expected)
  
  # Test with different cores
  result2 <- ASAP:::parallel_apply(test_indices, test_function, cores = 2)
  expect_equal(result2, expected)
})

test_that("sample_rows works correctly", {
  # Create test data
  test_df <- data.frame(
    id = 1:100,
    label = rep(c("A", "B", "C"), length.out = 100),
    value = rnorm(100),
    stringsAsFactors = FALSE
  )
  
  # Test basic sampling
  result <- ASAP:::sample_rows(test_df, sample_percent = 50, seed = 123)
  expect_equal(nrow(result), 50)
  expect_true(all(result$id %in% test_df$id))
  
  # Test balanced sampling
  result_balanced <- ASAP:::sample_rows(test_df, sample_percent = 30, 
                                       balanced = TRUE, seed = 123)
  expect_equal(nrow(result_balanced), 30)
  
  # Check that each label is represented
  label_counts <- table(result_balanced$label)
  expect_true(all(label_counts > 0))
  
  # Test with specific labels
  result_specific <- ASAP:::sample_rows(test_df, sample_percent = 50,
                                       labels = c("A", "B"), seed = 123)
  expect_true(all(result_specific$label %in% c("A", "B")))
  expect_false(any(result_specific$label == "C"))
})

test_that("check_pkg works correctly", {
  # Test with installed package
  expect_true(ASAP:::check_pkg("base"))
  
  # Test with non-existent package
  expect_false(ASAP:::check_pkg("nonexistent_package_12345"))
})

test_that("ensure_pkgs works correctly", {
  # Test with installed package (should not error)
  expect_silent(ASAP:::ensure_pkgs("base"))
  
  # Test with non-existent package (should error)
  expect_error(ASAP:::ensure_pkgs("nonexistent_package_12345"),
               "Required package not available")
})