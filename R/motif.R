
# Find motifs ------------------------------------------------------
# Update date : Feb. 7, 2025

#' Find Motifs in Song Data
#'
#' @description
#' Identifies and extracts motifs from song recordings based on detection times.
#'
#' @param x An object to process, either a data frame or SAP object
#' @param pre_time Time in seconds before detection point
#' @param lag_time Time in seconds after detection point
#' @param wav_dir For default method: Directory containing WAV files
#' @param add_path_attr For default method: Add wav_dir as attribute (default: TRUE)
#' @param template_name For SAP objects: Name of template to process
#' @param verbose Whether to print processing information (default: TRUE)
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For detection data frames:
#' \itemize{
#'   \item Requires columns: filename, time
#'   \item Validates motif boundaries against audio duration
#'   \item Processes each unique audio file
#'   \item Returns combined results with metadata
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Processes template-based detections
#'   \item Organizes results by recording day
#'   \item Validates motif boundaries
#'   \item Updates object with extracted motifs
#' }
#'
#' @return
#' For default method: Data frame containing:
#' \itemize{
#'   \item filename: Source WAV file name
#'   \item detection_time: Original detection time
#'   \item start_time, end_time: Motif boundaries
#'   \item duration: Motif duration
#' }
#'
#' For SAP objects: Updated object with motifs stored in motifs slot
#'
#' @examples
#' \dontrun{
#' # Find motifs from detection data frame
#' motifs <- find_motif(detections,
#'                      pre_time = 0.1,
#'                      lag_time = 0.2,
#'                      wav_dir = "path/to/wavs")
#'
#' # Process with path attribute
#' motifs <- find_motif(detections,
#'                      pre_time = 0.1,
#'                      lag_time = 0.2,
#'                      wav_dir = "path/to/wavs",
#'                      add_path_attr = TRUE)
#'
#' # Find motifs in SAP object
#' sap_obj <- find_motif(sap_object,
#'                       template_name = "template1",
#'                       pre_time = 0.7,
#'                       lag_time = 0.5)
#'
#' # Process with custom timing
#' sap_obj <- find_motif(sap_object,
#'                       template_name = "template2",
#'                       pre_time = 0.5,
#'                       lag_time = 0.3,
#'                       verbose = TRUE)
#' }
#'
#' @seealso \code{\link{detect_template}} for template detection
#'
#' @rdname find_motif
#' @export
find_motif <- function(x, ...) {
  UseMethod("find_motif")
}

#' @rdname find_motif
#' @export
find_motif.default <- function(x,
                               pre_time = NULL,
                               lag_time = NULL,
                               wav_dir = NULL,
                               add_path_attr = TRUE,
                               verbose = TRUE,
                               ...) {
  # Validate inputs
  if (!is.data.frame(x)) {
    stop("Input must be a data frame")
  }

  if (!all(c("filename", "time") %in% names(x))) {
    stop("Input data frame must contain 'filename' and 'time' columns")
  }

  if (is.null(wav_dir)) {
    stop("wav_dir must be provided")
  }

  if (is.null(pre_time) || is.null(lag_time)) {
    stop("pre_time and lag_time must be provided")
  }

  # Initialize results list
  motifs <- list()
  total_excluded <- 0

  # Process each unique file
  unique_files <- unique(x$filename)

  for (file in unique_files) {
    # Get full path to wav file
    wav_file <- file.path(wav_dir, file)

    if (!file.exists(wav_file)) {
      warning(sprintf("File not found: %s", wav_file))
      next
    }

    # Get wav duration
    wav_duration <- seewave::duration(wave = tuneR::readWave(wav_file))

    # Get detections for this file
    file_detections <- x[x$filename == file, ]

    # Calculate boundaries for each detection
    motif_data <- data.frame(
      filename = file_detections$filename,
      detection_time = file_detections$time,
      start_time = file_detections$time - pre_time,
      end_time = file_detections$time + lag_time
    )

    # Filter out boundaries that extend beyond wav file
    valid_motifs <- motif_data[
      motif_data$start_time >= 0 &
        motif_data$end_time <= wav_duration,
    ]

    if (nrow(valid_motifs) > 0) {
      # Add duration
      valid_motifs$duration <- valid_motifs$end_time - valid_motifs$start_time

      # Add to results list
      motifs[[file]] <- valid_motifs
    }

    # Count excluded motifs
    excluded <- nrow(motif_data) - nrow(valid_motifs)
    total_excluded <- total_excluded + excluded

    # Report processing if verbose
    if (verbose) {
      cat(sprintf("Processed %s: %d/%d valid motifs (%d excluded)\n",
                  file,
                  nrow(valid_motifs),
                  nrow(motif_data),
                  excluded))
    }
  }

  # Combine all results
  if (length(motifs) > 0) {
    final_motifs <- do.call(rbind, motifs)
    rownames(final_motifs) <- NULL

    # Add wav_dir as attribute if requested
    if (add_path_attr) {
      attr(final_motifs, "wav_dir") <- normalizePath(wav_dir, mustWork = TRUE)
    }

    if (verbose) {
      message(sprintf("Total valid motifs found: %d (excluded: %d)",
                      nrow(final_motifs), total_excluded))
    }
    return(final_motifs)
  } else {
    warning("No valid motifs found")
    return(NULL)
  }
}

#' @rdname find_motif
#' @export
find_motif.Sap <- function(x,
                           template_name,
                           pre_time = NULL,
                           lag_time = NULL,
                           verbose = TRUE,
                           ...) {
  if(verbose) message(sprintf("\n=== Starting Motif Extraction based on template '%s'===", template_name))

  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  if (is.null(template_name) || !template_name %in% names(x$templates$template_matches)) {
    stop("template_name must be provided and must exist in the template_matches")
  }

  if (is.null(pre_time) || is.null(lag_time)) {
    stop("pre_time and lag_time must be provided")
  }

  # Get template matches
  template_matches <- x$templates$template_matches[[template_name]]

  # Process each day separately
  days_to_process <- unique(template_matches$day_post_hatch)
  all_motifs <- list()
  daily_stats <- list()

  for (current_day in days_to_process) {
    # Get day-specific detections
    day_detections <- template_matches[template_matches$day_post_hatch == current_day, ]

    cat(sprintf("\nProcessing day %s:", current_day))

    # Use default method to find motifs with verbose turned off
    day_motifs <- find_motif.default(
      x = day_detections,
      pre_time = pre_time,
      lag_time = lag_time,
      wav_dir = file.path(x$base_path, current_day),
      add_path_attr = FALSE,
      verbose = FALSE,
      ...
    )

    if (!is.null(day_motifs)) {
      # Add day information if not already present
      if (!"day_post_hatch" %in% names(day_motifs)) {
        day_motifs$day_post_hatch <- current_day
      }
      if (!"label" %in% names(day_motifs)) {
        day_motifs$label <- day_detections$label[1]
      }

      all_motifs[[current_day]] <- day_motifs

      # Calculate excluded motifs
      excluded_motifs <- nrow(day_detections) - nrow(day_motifs)
      cat(sprintf("%d valid motifs, %d excluded",
                  nrow(day_motifs), excluded_motifs))

      # Store daily stats
      daily_stats[[current_day]] <- c(
        valid = nrow(day_motifs),
        excluded = excluded_motifs
      )
    }
  }

  # Combine all results
  if (length(all_motifs) > 0) {
    final_motifs <- do.call(rbind, all_motifs)
    rownames(final_motifs) <- NULL

    # Convert to segment object
    final_motifs <- as_segment(final_motifs)

    # Store results in SAP object
    x$motifs <- final_motifs

    # Update last modified time
    x$misc$last_modified <- Sys.time()

    # Calculate and report total statistics
    total_detections <- nrow(template_matches)
    total_motifs <- nrow(final_motifs)
    total_excluded <- total_detections - total_motifs

    message(sprintf("\nSummary:"))
    message(sprintf("Total detections: %d", total_detections))
    message(sprintf("Valid motifs: %d", total_motifs))
    message(sprintf("Excluded motifs: %d", total_excluded))
    message(sprintf("Motif results have been stored in the SAP object"))
  } else {
    warning("No valid motifs found")
  }

  # Return the modified SAP object invisibly
  invisible(x)
}
