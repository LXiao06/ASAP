
# Find motifs ------------------------------------------------------
# Update date : Feb. 7, 2025

#' Find Motifs in Song Data
#'
#' @description
#' Identifies and extracts motifs from song recordings based on detection times.
#' Supports single or multiple templates with template-specific pre- and lag-times.
#'
#' @param x An object to process, either a data frame of detections or a SAP object.
#' @param pre_time Time in seconds before detection point.
#'   Can be a single numeric value (broadcast to all templates) or a numeric vector
#'   matching \code{template_name} in length (optionally named).
#' @param lag_time Time in seconds after detection point.
#'   Can be a single numeric value (broadcast to all templates) or a numeric vector
#'   matching \code{template_name} in length (optionally named).
#' @param wav_dir For default method: Directory containing WAV files.
#' @param add_path_attr For default method: Add wav_dir as attribute (default: TRUE).
#' @param template_name For SAP objects: Character vector of template name(s) to process.
#'   If \code{NULL} (default), processes all templates found in \code{template_matches}.
#' @param day_post_hatch For SAP objects: Numeric value(s) to use for
#'   \code{day_post_hatch} when it cannot be determined from subfolder names
#'   (e.g. non-numeric folder names like \code{"FD_661_667"}).
#'   Can be a single numeric value (applied to all folders) or a named numeric
#'   vector mapping subfolder names to day-post-hatch values
#'   (e.g. \code{c(FD_661_667 = 61, PD_661_667 = 65, UD = 90)}).
#'   If \code{NULL} (default), the value is parsed from the subfolder name.
#' @param verbose Whether to print processing information (default: TRUE).
#' @param ... Additional arguments passed to specific methods.
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
#'   \item Processes template-based detections for one or multiple templates
#'   \item Applies template-specific \code{pre_time} and \code{lag_time} boundaries
#'   \item Organizes results by recording day / subfolder
#'   \item Validates motif boundaries against recording durations
#'   \item Stores all extracted motifs in a unified \code{x$motifs} segment data frame
#'     with a \code{template_name} column identifying the source template
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
#' For SAP objects: Updated SAP object with motifs stored in the \code{motifs} slot
#' as a validated \code{segment} object.
#'
#' @examples
#' \dontrun{
#' # Find motifs from detection data frame
#' motifs <- find_motif(detections,
#'                      pre_time = 0.1,
#'                      lag_time = 0.2,
#'                      wav_dir = "path/to/wavs")
#'
#' # Find motifs in SAP object with single template
#' sap_obj <- find_motif(sap_object,
#'                       template_name = "template1",
#'                       pre_time = 0.7,
#'                       lag_time = 0.5)
#'
#' # Find motifs with multiple templates and template-specific boundaries
#' sap_obj <- find_motif(sap_object,
#'                       template_name = c("c", "c1", "c2"),
#'                       pre_time = c(0.58, 0.78, 1.68),
#'                       lag_time = c(0.70, 0.87, 1.80))
#'
#' # Broadcast same pre_time and lag_time across multiple templates
#' sap_obj <- find_motif(sap_object,
#'                       template_name = c("c", "c1", "c2"),
#'                       pre_time = 0.5,
#'                       lag_time = 0.5)
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
    detection_time <- round(file_detections$time, 2)

    # Calculate boundaries for each detection
    motif_data <- data.frame(
      filename = file_detections$filename,
      detection_time = detection_time,
      start_time = detection_time - pre_time,
      end_time = detection_time + lag_time
    )

    # Carry over additional columns if present
    extra_cols <- setdiff(
      names(file_detections),
      c("filename", "time", "detection_time", "start_time", "end_time", "duration")
    )
    if (length(extra_cols) > 0) {
      for (col in extra_cols) {
        motif_data[[col]] <- file_detections[[col]]
      }
    }

    # Filter out boundaries that extend beyond wav file
    valid_motifs <- motif_data[
      which(motif_data$start_time >= 0 & motif_data$end_time <= wav_duration), ,
      drop = FALSE
    ]

    if (nrow(valid_motifs) > 0) {
      # Add duration
      valid_motifs$duration <- pre_time + lag_time

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
                           template_name = NULL,
                           pre_time = NULL,
                           lag_time = NULL,
                           day_post_hatch = NULL,
                           verbose = TRUE,
                           ...) {
  if (verbose) message(sprintf("\n=== Starting Motif Extraction ==="))

  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  if (is.null(x$templates$template_matches) || length(x$templates$template_matches) == 0) {
    stop("No template matches found in SAP object. Run detect_template() first.")
  }

  # If template_name is NULL, default to all available templates in template_matches
  if (is.null(template_name)) {
    template_name <- names(x$templates$template_matches)
  }

  if (length(template_name) == 0) {
    stop("template_name must be provided or exist in template_matches")
  }

  # Check that all specified templates exist in template_matches
  missing_tpls <- setdiff(template_name, names(x$templates$template_matches))
  if (length(missing_tpls) > 0) {
    stop("The following template(s) do not exist in template_matches: ",
         paste(missing_tpls, collapse = ", "))
  }

  # Validate and parse pre_time
  if (is.null(pre_time)) {
    stop("pre_time must be provided")
  }
  if (!is.numeric(pre_time)) {
    stop("pre_time must be numeric")
  }
  if (length(pre_time) == 1) {
    pre_times <- setNames(rep(pre_time, length(template_name)), template_name)
  } else if (length(pre_time) == length(template_name)) {
    if (!is.null(names(pre_time)) && all(names(pre_time) %in% template_name)) {
      pre_times <- pre_time[template_name]
    } else {
      pre_times <- setNames(pre_time, template_name)
    }
  } else {
    stop(sprintf(
      "length of pre_time (%d) must be 1 or match length of template_name (%d)",
      length(pre_time), length(template_name)
    ))
  }

  # Validate and parse lag_time
  if (is.null(lag_time)) {
    stop("lag_time must be provided")
  }
  if (!is.numeric(lag_time)) {
    stop("lag_time must be numeric")
  }
  if (length(lag_time) == 1) {
    lag_times <- setNames(rep(lag_time, length(template_name)), template_name)
  } else if (length(lag_time) == length(template_name)) {
    if (!is.null(names(lag_time)) && all(names(lag_time) %in% template_name)) {
      lag_times <- lag_time[template_name]
    } else {
      lag_times <- setNames(lag_time, template_name)
    }
  } else {
    stop(sprintf(
      "length of lag_time (%d) must be 1 or match length of template_name (%d)",
      length(lag_time), length(template_name)
    ))
  }

  all_motifs <- list()
  stats_per_template <- list()

  for (t_name in template_name) {
    t_pre <- pre_times[[t_name]]
    t_lag <- lag_times[[t_name]]
    t_matches <- x$templates$template_matches[[t_name]]

    if (is.null(t_matches) || nrow(t_matches) == 0) {
      if (verbose) {
        cat(sprintf("\nNo detections found for template '%s'. Skipping.\n", t_name))
      }
      stats_per_template[[t_name]] <- list(total_detections = 0, valid = 0, excluded = 0)
      next
    }

    if (verbose) {
      cat(sprintf("\nProcessing template '%s' (pre_time = %.2f, lag_time = %.2f):",
                  t_name, t_pre, t_lag))
    }

    # Determine unique days/subfolders to process
    days_to_process <- unique(t_matches$day_post_hatch)
    t_motifs <- list()

    for (current_day in days_to_process) {
      if (is.na(current_day)) {
        day_detections <- t_matches[is.na(t_matches$day_post_hatch), ]
        day_str <- "NA"
        sub_dir <- x$base_path
      } else {
        day_detections <- t_matches[!is.na(t_matches$day_post_hatch) & t_matches$day_post_hatch == current_day, ]
        day_str <- as.character(current_day)
        candidate_dir <- file.path(x$base_path, day_str)
        sub_dir <- if (dir.exists(candidate_dir)) candidate_dir else x$base_path
      }

      if (nrow(day_detections) == 0) next

      day_motifs <- find_motif.default(
        x = day_detections,
        pre_time = t_pre,
        lag_time = t_lag,
        wav_dir = sub_dir,
        add_path_attr = FALSE,
        verbose = FALSE,
        ...
      )

      if (!is.null(day_motifs) && nrow(day_motifs) > 0) {
        if (!"day_post_hatch" %in% names(day_motifs)) {
          day_motifs$day_post_hatch <- current_day
        }
        if (!"label" %in% names(day_motifs)) {
          day_motifs$label <- if ("label" %in% names(day_detections)) day_detections$label[1] else NA
        }
        day_motifs$subfolder <- day_str
        day_motifs$template_name <- t_name

        t_motifs[[day_str]] <- day_motifs

        excluded_motifs <- nrow(day_detections) - nrow(day_motifs)
        if (verbose) {
          cat(sprintf("\n  Day/subfolder '%s': %d valid motifs, %d excluded",
                      day_str, nrow(day_motifs), excluded_motifs))
        }
      } else {
        if (verbose) {
          cat(sprintf("\n  Day/subfolder '%s': 0 valid motifs, %d excluded",
                      day_str, nrow(day_detections)))
        }
      }
    }

    if (length(t_motifs) > 0) {
      t_all <- do.call(rbind, t_motifs)
      all_motifs[[t_name]] <- t_all
      stats_per_template[[t_name]] <- list(
        total_detections = nrow(t_matches),
        valid = nrow(t_all),
        excluded = nrow(t_matches) - nrow(t_all)
      )
    } else {
      stats_per_template[[t_name]] <- list(
        total_detections = nrow(t_matches),
        valid = 0,
        excluded = nrow(t_matches)
      )
    }
  }

  if (length(all_motifs) > 0) {
    final_motifs <- do.call(rbind, all_motifs)
    rownames(final_motifs) <- NULL

    # Fill NA day_post_hatch with user-supplied value(s)
    orig_day <- as.character(final_motifs$day_post_hatch)
    dph_numeric <- suppressWarnings(as.numeric(orig_day))
    if (any(is.na(dph_numeric))) {
      if (!is.null(day_post_hatch)) {
        na_idx <- which(is.na(dph_numeric))
        if (!is.null(names(day_post_hatch))) {
          for (i in na_idx) {
            folder <- orig_day[i]
            if (folder %in% names(day_post_hatch)) {
              dph_numeric[i] <- day_post_hatch[[folder]]
            } else {
              warning(sprintf("No day_post_hatch mapping for subfolder '%s'", folder))
            }
          }
        } else {
          dph_numeric[na_idx] <- day_post_hatch
        }
        n_filled <- sum(!is.na(dph_numeric[na_idx]))
        if (verbose) {
          message(sprintf("\nFilled %d day_post_hatch values from user-supplied argument", n_filled))
        }
      }
    }
    final_motifs$day_post_hatch <- dph_numeric

    # Convert to segment object
    final_motifs <- as_segment(final_motifs)

    # Store results in SAP object
    x$motifs <- final_motifs
    x$misc$last_modified <- Sys.time()

    if (verbose) {
      message(sprintf("\n\n=== Motif Extraction Summary ==="))
      for (t_name in names(stats_per_template)) {
        st <- stats_per_template[[t_name]]
        message(sprintf("Template '%s': %d detections -> %d valid motifs (%d excluded)",
                        t_name, st$total_detections, st$valid, st$excluded))
      }
      total_det <- sum(vapply(stats_per_template, function(s) s$total_detections, numeric(1)))
      total_val <- nrow(final_motifs)
      message(sprintf("Total: %d detections -> %d valid motifs (%d excluded)",
                      as.integer(total_det), total_val, as.integer(total_det - total_val)))
      message("Motif results stored in SAP object ($motifs).")
    }
  } else {
    warning("No valid motifs found across specified template(s)")
  }

  invisible(x)
}
