# ============================================================================
# Global Configuration for ASAP Shiny App
# ============================================================================
# This file loads packages, sets global options, and initializes the app.
# It runs once when the app starts.
# ============================================================================

# Load required packages --------------------------------------------------
library(shiny)
library(bslib)
library(ASAP)
library(httr2)
library(jsonlite)
library(plotly)
library(DT)
library(shinyWidgets)
library(digest)

# Source utility functions ------------------------------------------------
source("utils/ai_integration.R", local = TRUE)

# Global constants --------------------------------------------------------

# Maximum file upload size (100 MB)
options(shiny.maxRequestSize = 100 * 1024^2)

# Supported audio formats
SUPPORTED_FORMATS <- c(".wav", ".WAV")

# Default parameter values for zebra finch analysis
DEFAULT_PARAMS <- list(
  # Visualization
  freq_range = c(1, 10),  # kHz
  
  # Bout detection
  rms_threshold = 0.1,
  min_duration = 0.7,     # seconds
  bout_freq_range = c(1, 8),  # kHz
  
  # Syllable segmentation
  silence_threshold = 0.01,
  min_syllable_ms = 20,
  max_syllable_ms = 240,
  min_level_db = 10,
  
  # Template creation
  template_freq_min = 1,  # kHz
  template_freq_max = 10, # kHz
  template_threshold = 0.5,
  
  # Template detection
  proximity_window = 1,   # seconds
  
  # Motif extraction
  pre_time = 0.7,         # seconds
  lag_time = 0.5          # seconds
)

# App theme configuration -------------------------------------------------
app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#2C3E50",
  secondary = "#18BC9C",
  success = "#18BC9C",
  info = "#3498DB",
  warning = "#F39C12",
  danger = "#E74C3C",
  base_font = font_google("Inter"),
  heading_font = font_google("Inter"),
  code_font = font_google("Fira Code")
)

# Helper functions --------------------------------------------------------

#' Format file size for display
#'
#' @param size_bytes File size in bytes
#' @return Formatted string (e.g., "2.5 MB")
format_file_size <- function(size_bytes) {
  if (size_bytes < 1024) {
    return(paste(size_bytes, "B"))
  } else if (size_bytes < 1024^2) {
    return(paste(round(size_bytes / 1024, 1), "KB"))
  } else if (size_bytes < 1024^3) {
    return(paste(round(size_bytes / 1024^2, 1), "MB"))
  } else {
    return(paste(round(size_bytes / 1024^3, 1), "GB"))
  }
}

#' Validate WAV file
#'
#' @param file_path Path to WAV file
#' @return List with valid (logical) and message (character)
validate_wav_file <- function(file_path) {
  
  if (is.null(file_path) || !file.exists(file_path)) {
    return(list(valid = FALSE, message = "File not found"))
  }
  
  # Check file extension
  ext <- tools::file_ext(file_path)
  if (!paste0(".", ext) %in% SUPPORTED_FORMATS) {
    return(list(
      valid = FALSE, 
      message = paste("Unsupported format. Please use:", 
                     paste(SUPPORTED_FORMATS, collapse = ", "))
    ))
  }
  
  # Try to read file info
  tryCatch({
    info <- tuneR::readWave(file_path, header = TRUE)
    return(list(
      valid = TRUE,
      message = "Valid WAV file",
      sample_rate = info$sample.rate,
      duration = info$samples / info$sample.rate,
      channels = info$channels
    ))
  }, error = function(e) {
    return(list(
      valid = FALSE,
      message = paste("Error reading WAV file:", conditionMessage(e))
    ))
  })
}

#' Create notification message
#'
#' @param message Message text
#' @param type Type of notification ("default", "message", "warning", "error")
#' @return Shiny notification
show_notification <- function(message, type = "message") {
  # Map common type names to Shiny's expected values
  shiny_type <- switch(type,
    "info" = "message",
    "success" = "message",
    "warning" = "warning",
    "error" = "error",
    type  # default: use as-is
  )
  
  showNotification(
    message,
    type = shiny_type,
    duration = 5
  )
}

# Initialize app state ----------------------------------------------------
cat("ASAP Studio initialized\n")
cat("AI Model:", AI_CONFIG$model, "\n")
cat("Theme:", "Bootstrap 5 (Flatly)\n")
