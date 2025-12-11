# NPX Object for Neuropixels Recordings
# Update date: Dec. 10, 2025

# NPX Object Constructor --------------------------------------------------------

#' Create a Neuropixels (NPX) Object from Long Audio Recordings
#'
#' @description
#' Creates an NPX object for managing long audio recordings from Neuropixels
#' experiments, with support for chunked template matching.
#'
#' @param audio_path Character string specifying path to audio file or directory
#' @param labels Character vector of labels for recordings (optional)
#' @param metadata Optional pre-existing metadata data frame
#'
#' @details
#' This function creates an NPX object designed for long recordings that
#' require chunked processing. Unlike SAP objects which handle multiple
#' short recordings, NPX objects are optimized for single long recordings.
#'
#' @return An NPX object containing:
#'   \item{metadata}{Data frame with file and recording metadata}
#'   \item{base_path}{Directory containing the audio file(s)}
#'   \item{templates}{Template collection for template matching}
#'   \item{detections}{List of detection results by template}
#'   \item{misc}{List with creation details and timestamps}
#'
#' @examples
#' \dontrun{
#' # Create NPX object from single long recording
#' npx <- create_npx_object("path/to/long_recording.wav",
#'                          labels = "session1")
#'
#' # Create from directory with multiple recordings
#' npx <- create_npx_object("path/to/recordings/",
#'                          labels = c("pre", "post"))
#' }
#'
#' @seealso
#' \code{\link{detect_template.Npx}} for chunked template matching
#'
#' @export
create_npx_object <- function(audio_path,
                              labels = NULL,
                              metadata = NULL) {
  # Input validation
  if (!is.character(audio_path) || length(audio_path) != 1) {
    stop("audio_path must be a single character string")
  }

  # Normalize path
  audio_path <- normalizePath(audio_path, mustWork = TRUE)

  # Determine if path is file or directory
  is_directory <- file.info(audio_path)$isdir

  if (is_directory) {
    base_path <- audio_path
    # Find all wav files in directory
    wav_files <- list.files(audio_path, pattern = "\\.wav$",
                           full.names = FALSE, recursive = FALSE)
    if (length(wav_files) == 0) {
      stop("No WAV files found in directory: ", audio_path)
    }
  } else {
    # Single file
    base_path <- dirname(audio_path)
    wav_files <- basename(audio_path)
  }

  # Create metadata if not provided
  if (is.null(metadata)) {
    metadata <- data.frame(
      filename = wav_files,
      stringsAsFactors = FALSE
    )

    # Add duration for each file
    metadata$duration <- sapply(wav_files, function(f) {
      tryCatch({
        info <- av::av_media_info(file.path(base_path, f))
        info$duration
      }, error = function(e) NA_real_)
    })

    # Add labels if provided
    if (!is.null(labels)) {
      if (length(labels) == 1) {
        metadata$label <- labels
      } else if (length(labels) == length(wav_files)) {
        metadata$label <- labels
      } else {
        stop("Number of labels must be 1 or match number of files")
      }
    } else {
      metadata$label <- NA_character_
    }

    # Add index
    metadata$index <- seq_len(nrow(metadata))
  }

  # Create NPX object
  npx <- new_npx(
    metadata = metadata,
    base_path = base_path,
    misc = list(
      creation_date = Sys.time(),
      last_modified = Sys.time(),
      source_path = audio_path
    )
  )

  # Validate
  validate_npx(npx)

  # Report
  total_duration <- sum(metadata$duration, na.rm = TRUE)
  message(sprintf("Created NPX object with %d file(s), total duration: %.1f seconds (%.1f minutes)",
                  nrow(metadata),
                  total_duration,
                  total_duration / 60))

  return(npx)
}


#' Internal Constructor for NPX Object
#'
#' @param metadata Data frame with file metadata
#' @param base_path Character string of base directory
#' @param templates Template collection object
#' @param detections List of detection results
#' @param misc List of miscellaneous info
#' @param version Version string
#'
#' @return An NPX object of class "Npx"
#'
#' @keywords internal
new_npx <- function(metadata = data.frame(),
                    base_path = character(),
                    templates = new_template_collection(),
                    detections = list(),
                    misc = list(),
                    version = "1.0.0") {

  # Create structure
  x <- structure(
    list(
      metadata = metadata,
      base_path = base_path,
      templates = templates,
      detections = detections,
      misc = misc,
      version = version
    ),
    class = "Npx"
  )

  return(x)
}


#' Validate NPX Object
#'
#' @param x Object to validate
#'
#' @return TRUE if valid, otherwise throws error
#'
#' @keywords internal
validate_npx <- function(x) {
  # Check class

  if (!inherits(x, "Npx")) {
    stop("Object must be of class 'Npx'")
  }

  # Check required components
  required <- c("metadata", "base_path", "templates", "detections", "misc", "version")
  missing <- required[!required %in% names(x)]

  if (length(missing) > 0) {
    stop("Missing required components: ", paste(missing, collapse = ", "))
  }

  # Check metadata is data frame
  if (!is.data.frame(x$metadata)) {
    stop("metadata must be a data frame")
  }

  # Check base_path
  if (!is.character(x$base_path) || length(x$base_path) != 1) {
    stop("base_path must be a single character string")
  }

  # Check templates
  if (!inherits(x$templates, "template_collection")) {
    stop("templates must be a template_collection object")
  }

  # Check detections is list
  if (!is.list(x$detections)) {
    stop("detections must be a list")
  }

  return(TRUE)
}


#' Print method for NPX objects
#'
#' @param x An NPX object
#' @param ... Additional arguments (ignored)
#'
#' @exportS3Method print Npx
print.Npx <- function(x, ...) {
  cat("NPX Object (Neuropixels Audio)\n")
  cat("==============================\n")
  cat("Files:", nrow(x$metadata), "\n")

  total_dur <- sum(x$metadata$duration, na.rm = TRUE)
  cat(sprintf("Total duration: %.1f seconds (%.1f minutes)\n",
              total_dur, total_dur / 60))

  if (!is.null(x$metadata$label) && !all(is.na(x$metadata$label))) {
    cat("Labels:", paste(unique(na.omit(x$metadata$label)), collapse = ", "), "\n")
  }

  cat("Templates:", length(x$templates$template_list), "\n")

  n_detections <- sum(sapply(x$detections, function(d) {
    if (is.data.frame(d)) nrow(d) else 0
  }))
  cat("Total detections:", n_detections, "\n")

  cat("\nBase path:", x$base_path, "\n")
}


#' Summary method for NPX objects
#'
#' @param object An NPX object
#' @param ... Additional arguments (ignored)
#'
#' @exportS3Method summary Npx
summary.Npx <- function(object, ...) {
  cat("NPX Object Summary\n")
  cat("==================\n\n")

  # File info
  cat("Files:\n")
  print(object$metadata[, c("filename", "duration", "label")], row.names = FALSE)

  # Template info
  if (length(object$templates$template_list) > 0) {
    cat("\nTemplates:\n")
    for (tpl_name in names(object$templates$template_list)) {
      n_det <- 0
      if (!is.null(object$detections[[tpl_name]])) {
        n_det <- nrow(object$detections[[tpl_name]])
      }
      cat(sprintf("  - %s: %d detections\n", tpl_name, n_det))
    }
  } else {
    cat("\nNo templates defined.\n")
  }

  # Creation info
  cat("\nCreated:", format(object$misc$creation_date, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Version:", object$version, "\n")
}


# Split Audio Function ----------------------------------------------------------

#' Split Long Audio File into Chunks
#'
#' @description
#' Splits a long audio file into smaller chunks with optional overlap.
#' Used internally by detect_template.Npx for memory-efficient processing.
#'
#' @param audio_path Path to audio file
#' @param chunk_duration Duration of each chunk in seconds (default: 90)
#' @param overlap Overlap between consecutive chunks in seconds (default: 5)
#' @param output_dir Directory for temporary chunk files (default: tempdir())
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A data frame with columns:
#'   \item{chunk_path}{Path to the chunk audio file}
#'   \item{chunk_index}{1-indexed chunk number}
#'   \item{start_time}{Start time in original file (seconds)}
#'   \item{end_time}{End time in original file (seconds)}
#'   \item{duration}{Actual duration of chunk (seconds)}
#'
#' @examples
#' \dontrun{
#' chunks <- split_audio("long_recording.wav",
#'                       chunk_duration = 90,
#'                       overlap = 5)
#' print(chunks)
#' }
#'
#' @keywords internal
split_audio <- function(audio_path,
                        chunk_duration = 90,
                        overlap = 5,
                        output_dir = tempdir(),
                        verbose = TRUE) {

  # Validate inputs
  if (!file.exists(audio_path)) {
    stop("Audio file does not exist: ", audio_path)
  }

  if (overlap >= chunk_duration) {
    stop("Overlap must be less than chunk_duration")
  }

  # Get audio duration
  audio_info <- av::av_media_info(audio_path)
  total_duration <- audio_info$duration

  if (verbose) {
    cat(sprintf("Audio duration: %.1f seconds (%.1f minutes)\n",
                total_duration, total_duration / 60))
  }

  # If audio is shorter than chunk_duration, return single chunk info
  if (total_duration <= chunk_duration) {
    if (verbose) cat("Audio shorter than chunk size - processing as single chunk\n")
    return(data.frame(
      chunk_path = audio_path,
      chunk_index = 1L,
      start_time = 0,
      end_time = total_duration,
      duration = total_duration,
      stringsAsFactors = FALSE
    ))
  }

  # Calculate chunk boundaries
  step_size <- chunk_duration - overlap
  chunk_starts <- seq(0, total_duration - overlap, by = step_size)

  # Create output directory if needed
  chunk_dir <- file.path(output_dir, paste0("chunks_", basename(audio_path)))
  if (!dir.exists(chunk_dir)) {
    dir.create(chunk_dir, recursive = TRUE)
  }

  if (verbose) {
    cat(sprintf("Splitting into %d chunks (%.0fs each, %.0fs overlap)\n",
                length(chunk_starts), chunk_duration, overlap))
  }

  # Create chunks
  chunks <- data.frame(
    chunk_path = character(length(chunk_starts)),
    chunk_index = integer(length(chunk_starts)),
    start_time = numeric(length(chunk_starts)),
    end_time = numeric(length(chunk_starts)),
    duration = numeric(length(chunk_starts)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(chunk_starts)) {
    start_time <- chunk_starts[i]
    end_time <- min(start_time + chunk_duration, total_duration)
    actual_duration <- end_time - start_time

    # Create chunk filename
    chunk_filename <- sprintf("chunk_%03d_%.1f-%.1f.wav", i, start_time, end_time)
    chunk_path <- file.path(chunk_dir, chunk_filename)

    # Extract chunk using av
    av::av_audio_convert(
      audio_path,
      chunk_path,
      start_time = start_time,
      total_time = actual_duration
    )

    chunks$chunk_path[i] <- chunk_path
    chunks$chunk_index[i] <- i
    chunks$start_time[i] <- start_time
    chunks$end_time[i] <- end_time
    chunks$duration[i] <- actual_duration
  }

  if (verbose) {
    cat(sprintf("Created %d chunks in: %s\n", nrow(chunks), chunk_dir))
  }

  # Store chunk directory as attribute for cleanup
  attr(chunks, "chunk_dir") <- chunk_dir

  return(chunks)
}


# Detect Template for NPX -------------------------------------------------------

#' Detect Templates in NPX Object with Chunked Processing
#'
#' @description
#' Performs template matching on long audio recordings using chunked processing
#' for memory efficiency. Automatically handles timestamp correction and
#' deduplication of detections at chunk boundaries.
#'
#' @param x An NPX object
#' @param template_name Name of template to use (must exist in NPX object)
#' @param chunk_duration Duration of each chunk in seconds (default: 90)
#' @param overlap Overlap between chunks in seconds (default: 5)
#' @param threshold Optional new threshold value for detection
#' @param proximity_window Time window in seconds to filter nearby detections
#' @param cor.method Correlation method ("pearson" or "spearman")
#' @param cores Number of cores for parallel processing
#' @param cleanup Remove temporary chunk files after processing (default: TRUE)
#' @param save_plot Whether to save detection plots (default: FALSE)
#' @param verbose Print progress messages (default: TRUE)
#' @param ... Additional arguments passed to corMatch
#'
#' @return Updated NPX object with detection results stored in detections list
#'
#' @examples
#' \dontrun{
#' # Create NPX object and add template
#' npx <- create_npx_object("long_recording.wav")
#' npx <- create_template(npx, template_name = "motif1", ...)
#'
#' # Run chunked detection
#' npx <- detect_template(npx,
#'                        template_name = "motif1",
#'                        chunk_duration = 90,
#'                        overlap = 5)
#'
#' # Access results
#' detections <- npx$detections[["motif1"]]
#' }
#'
#' @rdname detect_template
#' @export
detect_template.Npx <- function(x,
                                template_name,
                                chunk_duration = 90,
                                overlap = 5,
                                threshold = NULL,
                                proximity_window = NULL,
                                cor.method = "pearson",
                                cores = NULL,
                                cleanup = TRUE,
                                save_plot = FALSE,
                                verbose = TRUE,
                                ...) {

  if (verbose) message("\n=== Starting Chunked Template Detection ===")

  # Validate inputs
  if (!inherits(x, "Npx")) {
    stop("Input must be an NPX object")
  }

  if (is.null(template_name) || !template_name %in% names(x$templates$template_list)) {
    stop("template_name must be provided and must exist in the NPX object")
  }

  # Get template
  template <- x$templates$template_list[[template_name]]

  # Update threshold if specified
  if (!is.null(threshold)) {
    original_threshold <- monitoR::templateCutoff(template)[[1]]
    if (original_threshold != threshold) {
      monitoR::templateCutoff(template) <- setNames(threshold, template_name)
      x$templates$template_list[[template_name]] <- template
      if (verbose) cat(sprintf("Template threshold updated: %.2f -> %.2f\n",
                               original_threshold, threshold))
    }
  }

  # Set up cores
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }

  # Process each file in metadata
  all_detections <- list()

  for (file_idx in seq_len(nrow(x$metadata))) {
    audio_file <- file.path(x$base_path, x$metadata$filename[file_idx])

    if (verbose) {
      cat(sprintf("\nProcessing file %d/%d: %s\n",
                  file_idx, nrow(x$metadata), basename(audio_file)))
    }

    # Split audio into chunks
    chunks <- split_audio(
      audio_path = audio_file,
      chunk_duration = chunk_duration,
      overlap = overlap,
      verbose = verbose
    )

    chunk_dir <- attr(chunks, "chunk_dir")

    # Process each chunk
    if (verbose) cat(sprintf("Running template matching on %d chunks...\n", nrow(chunks)))

    process_chunk <- function(i) {
      tryCatch({
        chunk_detections <- detect_template.default(
          x = chunks$chunk_path[i],
          template = template,
          cor.method = cor.method,
          save_plot = FALSE,
          proximity_window = NULL,  # Handle proximity at the end
          ...
        )

        if (!is.null(chunk_detections) && nrow(chunk_detections) > 0) {
          # Correct timestamps by adding chunk start time
          chunk_detections$time <- chunk_detections$time + chunks$start_time[i]
          chunk_detections$chunk_index <- chunks$chunk_index[i]
          chunk_detections$source_file <- x$metadata$filename[file_idx]
        }

        return(chunk_detections)
      }, error = function(e) {
        warning(sprintf("Error processing chunk %d: %s", i, e$message))
        return(NULL)
      })
    }

    # Process chunks (parallel or sequential)
    chunk_results <- parallel_apply(seq_len(nrow(chunks)), process_chunk, cores)

    # Combine results
    valid_results <- chunk_results[!sapply(chunk_results, is.null)]

    if (length(valid_results) > 0) {
      file_detections <- do.call(rbind, valid_results)

      # Sort by time
      file_detections <- file_detections[order(file_detections$time), ]

      # Apply proximity filtering to remove duplicates from overlapping regions
      if (!is.null(proximity_window) && nrow(file_detections) > 1) {
        file_detections <- filter_proximity_detections(file_detections, proximity_window)
      } else if (overlap > 0 && nrow(file_detections) > 1) {
        # Default: filter detections within overlap window
        file_detections <- filter_proximity_detections(file_detections, overlap)
      }

      # Add label if available
      if (!is.na(x$metadata$label[file_idx])) {
        file_detections$label <- x$metadata$label[file_idx]
      }

      all_detections[[x$metadata$filename[file_idx]]] <- file_detections

      if (verbose) {
        cat(sprintf("Found %d detections in %s\n",
                    nrow(file_detections), basename(audio_file)))
      }
    } else {
      if (verbose) cat("No detections found in this file\n")
    }

    # Cleanup chunk files
    if (cleanup && !is.null(chunk_dir) && dir.exists(chunk_dir)) {
      unlink(chunk_dir, recursive = TRUE)
      if (verbose) cat("Cleaned up temporary chunk files\n")
    }
  }

  # Combine all detections
  if (length(all_detections) > 0) {
    final_detections <- do.call(rbind, all_detections)
    row.names(final_detections) <- NULL

    x$detections[[template_name]] <- final_detections
    x$misc$last_modified <- Sys.time()

    if (verbose) {
      message(sprintf("\nTotal detections for '%s': %d",
                      template_name, nrow(final_detections)))
      message(sprintf("Access via: npx$detections[[\"%s\"]]", template_name))
    }
  } else {
    warning(sprintf("No detections found for template '%s'", template_name))
    x$detections[[template_name]] <- data.frame()
  }

  invisible(x)
}


#' Filter detections by proximity
#'
#' @param detections Data frame of detections with 'time' and 'score' columns
#' @param window Time window in seconds
#'
#' @return Filtered data frame with best detection per window
#'
#' @keywords internal
filter_proximity_detections <- function(detections, window) {
  if (nrow(detections) <= 1) return(detections)

  # Sort by time
  detections <- detections[order(detections$time), ]

  # Group by proximity
  time_diff <- diff(detections$time)
  group_ids <- cumsum(c(TRUE, time_diff > window))

  # Keep highest score per group
  detections$group <- group_ids

  result <- detections |>
    dplyr::group_by(.data$group) |>
    dplyr::slice_max(.data$score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-"group")

  return(as.data.frame(result))
}
