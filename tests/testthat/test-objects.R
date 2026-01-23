# Test SAP Object Creation and Validation
# Tests for functions in R/objects.R

test_that("parse_filename works correctly", {
  # Test valid filename parsing
  test_filename <- "/path/to/bird1_123.45_01_15_10_30_45.wav"
  result <- ASAP:::parse_filename(test_filename)
  
  expect_type(result, "list")
  expect_named(result, c("bird_id", "recording_date", "recording_time"))
  expect_equal(result$bird_id, "bird1")
  expect_equal(result$recording_date, "01-15")
  expect_s3_class(result$recording_time, "POSIXlt")
  
  # Test invalid filename
  expect_error(ASAP:::parse_filename("invalid_filename.wav"))
})

test_that("new_sap creates valid SAP object structure", {
  sap_obj <- ASAP:::new_sap()
  
  # Check class
  expect_s3_class(sap_obj, "Sap")
  
  # Check required components
  required_components <- c("metadata", "base_path", "motifs", "bouts",
                          "syllables", "segments", "templates", 
                          "features", "misc", "version")
  expect_named(sap_obj, required_components)
  
  # Check component types
  expect_s3_class(sap_obj$metadata, "data.frame")
  expect_type(sap_obj$base_path, "character")
  expect_s3_class(sap_obj$templates, "template_collection")
  expect_s3_class(sap_obj$motifs, "segment")
  expect_s3_class(sap_obj$bouts, "segment")
  expect_s3_class(sap_obj$syllables, "segment")
  expect_s3_class(sap_obj$segments, "segment")
  expect_type(sap_obj$features, "list")
  expect_type(sap_obj$misc, "list")
  expect_type(sap_obj$version, "character")
})

test_that("validate_sap correctly validates SAP objects", {
  # Valid SAP object should pass
  valid_sap <- ASAP:::new_sap()
  expect_true(ASAP:::validate_sap(valid_sap))
  
  # Invalid objects should fail
  expect_error(ASAP:::validate_sap(list()), "Object must be of class 'Sap'")
  
  # Missing components
  incomplete_sap <- valid_sap
  incomplete_sap$metadata <- NULL
  expect_error(ASAP:::validate_sap(incomplete_sap), "Missing required components")
  
  # Wrong component types
  wrong_type_sap <- valid_sap
  wrong_type_sap$metadata <- "not a data frame"
  expect_error(ASAP:::validate_sap(wrong_type_sap), "metadata must be a data frame")
  
  wrong_type_sap2 <- valid_sap
  wrong_type_sap2$base_path <- c("path1", "path2")
  expect_error(ASAP:::validate_sap(wrong_type_sap2), "base_path must be a single character string")
})

test_that("create_sap_metadata input validation works", {
  # Test missing labels
  expect_error(create_sap_metadata("dummy_path"), 
               "argument \"labels\" is missing")
  
  # Test invalid labels
  expect_error(create_sap_metadata("dummy_path", labels = character(0)),
               "labels must be a non-empty character vector")
  
  expect_error(create_sap_metadata("dummy_path", labels = 123),
               "labels must be a non-empty character vector")
})

test_that("create_sap_object input validation works", {
  # Test invalid base_path
  expect_error(create_sap_object(c("path1", "path2"), labels = "test"),
               "base_path must be a single character string")
  
  expect_error(create_sap_object(123, labels = "test"),
               "base_path must be a single character string")
  
  # Test invalid labels
  expect_error(create_sap_object("dummy_path", labels = character(0)),
               "labels must be a non-empty character vector")
  
  # Test non-existent path
  expect_error(create_sap_object("/non/existent/path", labels = "test"),
               "Invalid base_path")
})