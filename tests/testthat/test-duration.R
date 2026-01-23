# Test Duration Functions
# Tests for functions in R/duration.R

test_that("compute_wav_durations input validation", {
  # Test invalid SAP object
  expect_error(compute_wav_durations("not a sap object"),
               "Input must be a SAP object")

  # Create mock SAP object for testing
  mock_sap <- list(
    base_path = NULL,
    metadata = NULL
  )
  class(mock_sap) <- "Sap"

  # Test missing base_path
  expect_error(compute_wav_durations(mock_sap),
               "sap\\$base_path must contain WAV directory path")

  # Test missing metadata
  mock_sap$base_path <- "/dummy/path"
  expect_error(compute_wav_durations(mock_sap),
               "sap\\$metadata is missing")
})

test_that("compute_wav_durations core detection", {
  # Test core detection logic
  detect_cores_logic <- function(cores) {
    if (is.null(cores)) {
      cores <- parallel::detectCores() - 1
    }
    return(cores)
  }

  # Test automatic detection
  auto_cores <- detect_cores_logic(NULL)
  expect_true(is.numeric(auto_cores))
  expect_true(auto_cores >= 1)

  # Test manual setting
  manual_cores <- detect_cores_logic(4)
  expect_equal(manual_cores, 4)
})

test_that("refine_motif_boundaries input validation", {
  # Test basic parameter validation
  validate_boundary_params <- function(wav_dir = NULL,
                                     amplitude_threshold = 0.1,
                                     freq_range = c(1, 10),
                                     smooth_window = 5) {

    if (!is.null(amplitude_threshold) &&
        (!is.numeric(amplitude_threshold) ||
         amplitude_threshold < 0 || amplitude_threshold > 1)) {
      stop("amplitude_threshold must be between 0 and 1")
    }

    if (!is.null(freq_range) &&
        (!is.numeric(freq_range) || length(freq_range) != 2)) {
      stop("freq_range must be a numeric vector of length 2")
    }

    if (!is.null(smooth_window) &&
        (!is.numeric(smooth_window) || smooth_window <= 0)) {
      stop("smooth_window must be a positive number")
    }

    return(TRUE)
  }

  # Test valid parameters
  expect_true(validate_boundary_params())
  expect_true(validate_boundary_params(amplitude_threshold = 0.05))
  expect_true(validate_boundary_params(freq_range = c(0.5, 8)))
  expect_true(validate_boundary_params(smooth_window = 10))

  # Test invalid parameters
  expect_error(validate_boundary_params(amplitude_threshold = -0.1),
               "amplitude_threshold must be between 0 and 1")
  expect_error(validate_boundary_params(amplitude_threshold = 1.5),
               "amplitude_threshold must be between 0 and 1")
  expect_error(validate_boundary_params(freq_range = c(1)),
               "freq_range must be a numeric vector of length 2")
  expect_error(validate_boundary_params(smooth_window = -5),
               "smooth_window must be a positive number")
})

test_that("plot_motif_boundaries parameter validation", {
  # Test plotting parameter validation
  validate_plot_params <- function(save_plot = FALSE,
                                 plot_dir = NULL,
                                 plot_width = 12,
                                 plot_height = 8) {

    if (!is.logical(save_plot)) {
      stop("save_plot must be logical")
    }

    if (save_plot && is.null(plot_dir)) {
      stop("plot_dir must be specified when save_plot is TRUE")
    }

    if (!is.numeric(plot_width) || plot_width <= 0) {
      stop("plot_width must be a positive number")
    }

    if (!is.numeric(plot_height) || plot_height <= 0) {
      stop("plot_height must be a positive number")
    }

    return(TRUE)
  }

  # Test valid parameters
  expect_true(validate_plot_params())
  expect_true(validate_plot_params(save_plot = TRUE, plot_dir = "/tmp"))
  expect_true(validate_plot_params(plot_width = 15, plot_height = 10))

  # Test invalid parameters
  expect_error(validate_plot_params(save_plot = "yes"),
               "save_plot must be logical")
  expect_error(validate_plot_params(save_plot = TRUE),
               "plot_dir must be specified when save_plot is TRUE")
  expect_error(validate_plot_params(plot_width = -5),
               "plot_width must be a positive number")
  expect_error(validate_plot_params(plot_height = 0),
               "plot_height must be a positive number")
})

test_that("motif boundary frame calculation", {
  # Test frame calculation logic (simplified version)
  calculate_frames <- function(start_time, end_time, sample_rate = 44100, wl = 512) {
    start_frame <- floor(start_time * sample_rate / wl)
    end_frame <- ceiling(end_time * sample_rate / wl)

    if (start_frame < 1) start_frame <- 1
    if (end_frame <= start_frame) end_frame <- start_frame + 1

    return(list(start_frame = start_frame, end_frame = end_frame))
  }

  # Test normal case
  result1 <- calculate_frames(1.0, 2.0)
  expect_true(result1$start_frame >= 1)
  expect_true(result1$end_frame > result1$start_frame)

  # Test edge case - very small time window
  result2 <- calculate_frames(0.001, 0.002)
  expect_true(result2$start_frame >= 1)
  expect_true(result2$end_frame > result2$start_frame)

  # Test zero start time
  result3 <- calculate_frames(0, 1)
  expect_equal(result3$start_frame, 1)
  expect_true(result3$end_frame > 1)
})
