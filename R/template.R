
# Create Audio Clips------------------------------------------------------
# Update date : Feb. 7, 2025

#' Create Audio Clips from Sound Files
#'
#' @description
#' Creates audio clips from WAV files or SAP objects by extracting specified time segments.
#'
#' @param x An object to process, either a file path or SAP object
#' @param start_time Numeric start time(s) of the clip(s)
#' @param end_time Numeric end time(s) of the clip(s)
#' @param clip_name,clip_names Name(s) for the output clip(s)
#' @param unit Time unit ("second" or "millisecond")
#' @param indices For SAP objects: Numeric vector of indices to process
#' @param verbose For SAP objects: Whether to print progress messages
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For single WAV files:
#' \itemize{
#'   \item Validates input file and time parameters
#'   \item Creates templates directory if needed
#'   \item Extracts specified segment from audio file
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Creates clips for specified indices
#'   \item Updates template information in SAP object
#'   \item Maintains metadata about created clips
#' }
#'
#' @return
#' For default method: Character string containing path to created audio clip
#' For SAP objects: Updated SAP object with new template information
#'
#' @examples
#' \dontrun{
#' # Create clip from single WAV file
#' create_audio_clip("path/to/song.wav",
#'                   start_time = 10,
#'                   end_time = 20,
#'                   clip_name = "song_clip")
#'
#' # Create multiple clips from SAP object
#' create_audio_clip(sap_object,
#'                   indices = c(1, 2),
#'                   start_time = c(10, 20),
#'                   end_time = c(20, 30),
#'                   clip_names = c("clip1", "clip2"))
#'
#' # Create clip with millisecond units
#' create_audio_clip("song.wav",
#'                   start_time = 10000,
#'                   end_time = 20000,
#'                   unit = "millisecond")
#' }
#'
#' @rdname create_audio_clip
#' @export
create_audio_clip <- function(x, ...) {
  UseMethod("create_audio_clip")
}

#' @rdname create_audio_clip
#' @export
create_audio_clip.default <- function(x,
                                      ...,
                                      start_time,
                                      end_time,
                                      clip_name = NULL,
                                      unit = "second") {
  # Match unit argument
  unit <- match.arg(unit, choices = c("second", "millisecond"))

  # Normalize and validate file path
  x <- normalizePath(x, mustWork = TRUE)
  if (!file.exists(x)) {
    stop("File does not exist: ", x)
  }

  # Validate times
  if (is.null(start_time) || is.null(end_time)) {
    stop("start_time and end_time must be provided")
  }

  # Convert times to seconds if needed
  if (unit == "millisecond") {
    start_time <- start_time / 1000
    end_time <- end_time / 1000
  }

  # Set up output path and clip name
  templates_path <- file.path(dirname(x), "templates")
  if (!dir.exists(templates_path)) {
    dir.create(templates_path)
  }
  templates_path <- normalizePath(templates_path)

  # Generate clip name if not provided
  if (is.null(clip_name)) {
    clip_name <- paste0("clip_", basename(x))
  }

  # Construct full path for output clip
  clip_path <- file.path(templates_path, clip_name)

  # Create audio clip
  av::av_audio_convert(
    x,
    clip_path,
    start_time = start_time,
    total_time = end_time - start_time
  )

  # Normalize output path after creation
  clip_path <- normalizePath(clip_path)

  cat("Song clip", basename(clip_path), "is generated.\n")
  return(clip_path)
}

#' @rdname create_audio_clip
#' @return
#' Updated SAP object with new template information
#'
#' @export
create_audio_clip.Sap <- function(x,
                                  indices,
                                  start_time,
                                  end_time,
                                  clip_names,
                                  unit = "second",
                                  verbose = TRUE,
                                  ...) {
  if(verbose) message(sprintf("\n=== Starting Audio Clip Creation ===\n"))

  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  # Match unit argument
  unit <- match.arg(unit, choices = c("second", "millisecond"))

  # Validate inputs
  if (is.null(indices) || is.null(clip_names)) {
    stop("indices and clip_names must be provided")
  }

  if (length(start_time) != length(indices) ||
      length(end_time) != length(indices) ||
      length(clip_names) != length(indices)) {
    stop("start_time, end_time, and clip_names must match the length of indices")
  }

  # Remove .wav if users accidentally included it
  clip_names <- file_path_sans_ext(clip_names)

  # Set up templates directory in root directory
  templates_path <- file.path(x$base_path, "templates")
  if (!dir.exists(templates_path)) {
    dir.create(templates_path)
  }

  # Process each index
  for (i in seq_along(indices)) {
    index <- indices[i]

    # Get full path to source file
    song_file <- file.path(x$base_path,
                           x$metadata$day_post_hatch[index],
                           x$metadata$filename[index])

    # Create clip using default method with templates_path override
    clip_path <- file.path(templates_path, paste0(clip_names[i], ".wav"))

    # Create audio clip
    av::av_audio_convert(
      song_file,
      clip_path,
      start_time = if(unit == "millisecond") start_time[i]/1000 else start_time[i],
      total_time = if(unit == "millisecond") (end_time[i] - start_time[i])/1000
      else end_time[i] - start_time[i]
    )

    cat(sprintf("Audio clip '%s' is created \nOutput path: %s\n\n",
                clip_names[i],
                clip_path))

    # Create new row
    new_row <- data.frame(
      template_name = NA_character_,
      start_time = NA,
      end_time = NA,
      duration = NA,
      freq_min = NA,
      freq_max = NA,
      threshold = NA,
      clip_name = clip_names[i],
      clip_path = clip_path,
      source_file = x$metadata$filename[index],
      source_file_path = song_file,
      creation_date = Sys.time(),
      stringsAsFactors = FALSE
    )

    # Check if clip_name already exists
    existing_row <- which(x$templates$template_info$clip_name == clip_names[i])

    if (length(existing_row) > 0) {
      # Override existing row
      x$templates$template_info[existing_row, ] <- new_row
    } else {
      # Add new row
      x$templates$template_info <- rbind(x$templates$template_info, new_row)
    }
  }

  return(x)
}


# Create Templates ------------------------------------------------------
# Update date : Feb. 7, 2025

#' Create Correlation Templates for Song Analysis
#'
#' @description
#' Creates correlation templates from WAV files or SAP objects for song detection and analysis.
#'
#' @param x An object to process, either a file path or SAP object
#' @param template_name Character name for the template
#' @param start_time Numeric start time of template segment
#' @param end_time Numeric end time of template segment
#' @param freq_min Numeric minimum frequency in kHz (default: 0)
#' @param freq_max Numeric maximum frequency in kHz (default: 15)
#' @param threshold Numeric correlation threshold (default: 0.6)
#' @param write_template Logical whether to write template to disk
#' @param clip_name For SAP objects: Character name of the clip to use
#' @param verbose For SAP objects: Whether to print progress messages
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For WAV files:
#' \itemize{
#'   \item Validates input parameters
#'   \item Creates template using monitoR package
#'   \item Optionally writes template to disk
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Validates clip existence
#'   \item Creates template using specified parameters
#'   \item Updates SAP object with template information
#' }
#'
#' @return
#' For default method: A correlation template object from monitoR package
#' For SAP objects: Updated SAP object with new template information
#'
#' @examples
#' \dontrun{
#' # Create template from WAV file
#' template <- create_template("path/to/song.wav",
#'                            template_name = "template1",
#'                            start_time = 1.0,
#'                            end_time = 2.0,
#'                            freq_min = 2,
#'                            freq_max = 8)
#'
#' # Create and save template
#' template <- create_template("song.wav",
#'                            template_name = "template2",
#'                            start_time = 1.0,
#'                            end_time = 2.0,
#'                            write_template = TRUE)
#'
#' # Create template from SAP object
#' sap_obj <- create_template(sap_object,
#'                           template_name = "template1",
#'                           clip_name = "clip1",
#'                           freq_min = 2,
#'                           freq_max = 8)
#' }
#'
#' @seealso \code{\link{create_audio_clip}} for creating audio clips
#'
#' @rdname create_template
#' @export
create_template <- function(x, ...) {
  UseMethod("create_template")
}

#' @rdname create_template
#' @export
create_template.default <- function(x,              # x is wav file path
                                    template_name,    # mandatory
                                    start_time = NULL,
                                    end_time = NULL,
                                    freq_min = 0,
                                    freq_max = 15,
                                    threshold = 0.6,
                                    write_template = FALSE,
                                    ...) {
  # Validate file path
  if (!file.exists(x)) {
    stop("File does not exist: ", x)
  }

  # Validate template name
  if (missing(template_name) || is.null(template_name)) {
    stop("template_name must be provided")
  }

  # Validate time parameters if provided
  if ((!is.null(start_time) && is.null(end_time)) ||
      (is.null(start_time) && !is.null(end_time))) {
    stop("Both start_time and end_time must be provided together")
  }

  if (!is.null(start_time) && !is.null(end_time)) {
    if (start_time >= end_time) {
      stop("start_time must be less than end_time")
    }
  }

  # Validate frequency parameters if provided
  if ((!is.null(freq_min) && is.null(freq_max)) ||
      (is.null(freq_min) && !is.null(freq_max))) {
    stop("Both freq_min and freq_max must be provided together")
  }

  if (!is.null(freq_min) && !is.null(freq_max)) {
    if (freq_min >= freq_max) {
      stop("freq_min must be less than freq_max")
    }
  }

  # Validate threshold
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be a number between 0 and 1")
  }

  # Create template
  if (is.null(start_time) && is.null(end_time)) {
    template <- monitoR::makeCorTemplate(x,
                                         name = template_name,
                                         score.cutoff = threshold)
  } else {
    # If ANY parameters are provided, ALL must be provided
    if (any(is.null(c(start_time, end_time, freq_min, freq_max)))) {
      stop("All parameters (start_time, end_time, freq_min, freq_max) must be provided together")
    }

    # Create template with all parameters
    template <- monitoR::makeCorTemplate(x,
                                         t.lim = c(start_time, end_time),
                                         frq.lim = c(freq_min, freq_max),
                                         name = template_name,
                                         score.cutoff = threshold)
  }

  # Write template if requested
  if (write_template) {
    monitoR::writeCorTemplates(template, dir = dirname(x))
  }

  return(template)
}

#' @rdname create_template
#' @export
create_template.Sap <- function(x,              # x is Sap object
                                template_name,    # mandatory
                                clip_name,        # mandatory, to identify which clip to use
                                start_time = NULL,
                                end_time = NULL,
                                freq_min = 0,
                                freq_max = 15,
                                threshold = 0.6,
                                write_template = FALSE,
                                verbose = TRUE,
                                ...) {
  if(verbose) message(sprintf("\n=== Starting Template Creation ===\n"))

  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  # Validate template name
  if (missing(template_name) || is.null(template_name)) {
    stop("template_name must be provided")
  }

  # Validate clip name
  if (missing(clip_name) || is.null(clip_name)) {
    stop("clip_name must be provided")
  }

  # Check if clip exists in template_info
  clip_row <- match(clip_name, x$templates$template_info$clip_name)
  if (is.na(clip_row)) {
    stop("Clip not found: ", clip_name)
  }

  # Get clip path
  clip_path <- x$templates$template_info$clip_path[clip_row]

  # Create template using default method
  template <- create_template.default(
    x = clip_path,
    template_name = template_name,
    start_time = start_time,
    end_time = end_time,
    freq_min = freq_min,
    freq_max = freq_max,
    threshold = threshold,
    write_template = write_template
  )

  # Update template_info only if template is written
  if (write_template) {
    new_row <- data.frame(
      template_name = template_name,
      start_time = start_time,
      end_time = end_time,
      duration = end_time - start_time,
      freq_min = freq_min,
      freq_max = freq_max,
      threshold = threshold,
      clip_name = clip_name,
      clip_path = clip_path,
      source_file = x$templates$template_info$source_file[clip_row],
      source_file_path = x$templates$template_info$source_file_path[clip_row],
      creation_date = Sys.time(),
      stringsAsFactors = FALSE
    )

    # Check if template name already exists
    template_row <- which(x$templates$template_info$template_name == template_name)

    if (length(template_row) > 0) {
      # Override existing row
      x$templates$template_info[template_row, ] <- new_row
    } else {
      # Add new row
      x$templates$template_info <- rbind(x$templates$template_info, new_row)
    }

    # Store S4 template object
    x$templates$template_list[[template_name]] <- template

    # Initialize empty matches data frame for this template if it doesn't exist
    if (is.null(x$templates$template_matches[[template_name]])) {
      x$templates$template_matches[[template_name]] <- data.frame()
    }
  }

  return(x)
}


# Detect Templates ------------------------------------------------------
# Update date : Feb. 7, 2025

#' Detect Templates in Song Data
#'
#' @description
#' Performs template matching on audio files using correlation-based detection.
#'
#' @param x An object to process, either a file path or SAP object
#' @param template For default method: A template object created by create_template()
#' @param cor.method Correlation method ("pearson" or "spearman")
#' @param save_plot Whether to save detection plots
#' @param plot_dir For default method: Directory to save plots
#' @param proximity_window Time window in seconds to filter nearby detections (NULL to disable filtering).
#'        Only the detection with the highest score within each window is retained.
#' @param template_name For SAP objects: Name of template to use
#' @param day For SAP objects: Numeric vector of days to process
#' @param indices For SAP objects: Numeric vector of indices to process
#' @param threshold For SAP objects: New threshold value
#' @param cores For SAP objects: Number of cores for parallel processing
#' @param plot_percent For SAP objects: Percentage of files to plot (default: 10)
#' @param verbose For SAP objects: Whether to print progress messages
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For WAV files:
#' \itemize{
#'   \item Validates input file and template
#'   \item Performs correlation matching
#'   \item Finds peaks in correlation scores
#'   \item Optionally filters nearby detections (if proximity_window is set)
#'   \item Optionally saves detection plots
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Parallel processing support
#'   \item Day-specific processing
#'   \item Optional threshold adjustment
#'   \item Progress tracking and reporting
#'   \item Selective plot generation
#'   \item Filtering of nearby detections when proximity_window is specified
#' }
#'
#' @section Proximity Filtering:
#' When \code{proximity_window} is specified, the function will filter detections that occur
#' within the specified time window (in seconds). For each group of detections within the window,
#' only the one with the highest score is kept. This is useful for removing false positive detections.
#'
#' @return
#' For default method: A data frame containing detection results with columns:
#' \itemize{
#'   \item filename: Name of the processed file
#'   \item time: Time point of detection
#'   \item score: Correlation score
#' }
#'
#' For SAP objects: Updated SAP object with detection results stored in template_matches
#'
#' @examples
#' \dontrun{
#' # Detect template in single WAV file
#' detections <- detect_template("path/to/song.wav",
#'                              template = template_obj,
#'                              save_plot = TRUE)
#'
#' # Detect template in SAP object
#' sap_obj <- detect_template(sap_object,
#'                           template_name = "template1",
#'                           day = c(30, 40),
#'                           threshold = 0.7,
#'                           cores = 4)
#'
#' # Process specific indices with plots
#' sap_obj <- detect_template(sap_object,
#'                           template_name = "template1",
#'                           indices = 1:10,
#'                           save_plot = TRUE)
#'
#' # Filter nearby detections within 0.5 seconds
#' sap_obj <- detect_template(sap_object,
#'                           template_name = "template1",
#'                           proximity_window = 0.5)
#' }
#'
#' @seealso \code{\link{create_template}} for creating templates
#'
#' @rdname detect_template
#' @export
detect_template <- function(x, ...) {
  UseMethod("detect_template")
}

#' @rdname detect_template
#' @export
detect_template.default <- function(x,  # x is wav file path
                                    template,
                                    cor.method = "pearson",
                                    save_plot = FALSE,
                                    plot_dir = NULL,
                                    proximity_window = NULL,
                                    ...) {
  # Validate inputs
  if (!file.exists(x)) stop("File does not exist: ", x)
  if (is.null(template)) stop("template must be provided")

  # Set default plot directory
  if (save_plot && is.null(plot_dir)) {
    plot_dir <- file.path(dirname(x), "plots", "template_matches")
  }

  # Suppress monitoR output
  null_con <- file(tempfile(), open = "w")
  on.exit({
    sink(type = "output"); sink(type = "message"); close(null_con)
  })
  sink(null_con, type = "output"); sink(null_con, type = "message")

  # Perform correlation matching
  scores <- suppressWarnings(
    monitoR::corMatch(
      survey = x,
      templates = template,
      show.prog = FALSE,
      cor.method = cor.method,
      time.source = "fileinfo"
    )
  )

  # Find peaks with or without proximity filtering
  if (is.null(proximity_window)) {
    pks <- monitoR::findPeaks(score.obj = scores)
  } else {
    pks <- find_peaks_with_proximity(scores, proximity_window)
  }

  # Get detections
  detections <- monitoR::getDetections(pks, id = basename(x))

  # Early return if no detections
  if (is.null(detections) || nrow(detections) == 0) return(NULL)

  # Process final output
  detections <- detections %>%
    dplyr::rename(filename = .data$id) %>%
    dplyr::select(-.data$date.time)

  # Generate plot with filtered detections
  if (save_plot) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    plot_file <- file.path(plot_dir, paste0(tools::file_path_sans_ext(basename(x)), ".png"))

    tryCatch({
      # Open device
      png(plot_file, width = 1200, height = 800, res = 150)
      on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)

      # Get plot method once per function call
      if (!exists("plot_detections")) {
        plot_detections <- methods::getMethod("plot", "detectionList", where = asNamespace("monitoR"))
      }

      # Plot
      suppressMessages(plot_detections(pks, ask = FALSE))

      # Close device
      if (dev.cur() > 1) dev.off()
    }, error = function(e) {
      warning("Error creating plot for ", basename(x), ": ", e$message)
      if (dev.cur() > 1) dev.off()
    })
  }


  return(detections)
}

#' @rdname detect_template
#' @export
detect_template.Sap <- function(x,  # x is SAP object
                                day = NULL,
                                indices = NULL,
                                template_name,
                                threshold = NULL,
                                cores = NULL,
                                cor.method = "pearson",
                                save_plot = FALSE,
                                plot_percent = 10,
                                verbose = TRUE,
                                proximity_window = NULL,
                                ...) {
  if(verbose) message(sprintf("\n=== Starting Template Detection ==="))

  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  if (is.null(template_name) || !template_name %in% names(x$templates$template_list)) {
    stop("template_name must be provided and must exist in the SAP object")
  }

  # Get the template from SAP object
  template <- x$templates$template_list[[template_name]]

  # Update threshold if specified
  if (!is.null(threshold)) {
    original_threshold <- monitoR::templateCutoff(template)[[1]]

    # Only update if different
    if (original_threshold != threshold) {
      monitoR::templateCutoff(template) <- setNames(threshold, template_name)
      x$templates$template_list[[template_name]] <- template
      threshold_idx <- which(x$templates$template_info$template_name == template_name)
      if (length(threshold_idx) > 0) {
        x$templates$template_info$threshold[threshold_idx] <- threshold
      }
      cat(sprintf("\nTemplate threshold updated from %.2f to %.2f\n",
                  original_threshold, threshold))
    } else {
      cat(sprintf("\nTemplate threshold already at %.2f\n", threshold))
    }
  }

  # Filter metadata based on day
  if (!is.null(day)) {
    process_metadata <- x$metadata[x$metadata$day_post_hatch %in% day, ]
    days_to_process <- day
    if (nrow(process_metadata) == 0) {
      stop("No files found for specified day(s)")
    }
  } else {
    process_metadata <- x$metadata
    days_to_process <- unique(process_metadata$day_post_hatch)
  }

  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }

  # Create plot directories if save_plot is TRUE
  if (save_plot) {
    plots_dir <- file.path(x$base_path, "plots", "template_matches")
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    for (d in days_to_process) {
      day_dir <- file.path(plots_dir, paste0("day_", d))
      dir.create(day_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  # Process each day
  all_results <- list()

  for (current_day in days_to_process) {
    day_metadata <- process_metadata[process_metadata$day_post_hatch == current_day, ]
    if (!is.null(indices)) {
      valid_indices <- indices[indices <= nrow(day_metadata)]
      if (length(valid_indices) > 0) {
        day_metadata <- day_metadata[valid_indices, ]
      } else {
        cat(sprintf("\nNo valid indices for day %s.\n", current_day))
        next
      }
    }

    unique_files <- which(!duplicated(day_metadata$filename))
    cat(sprintf("\nProcessing %d files for day %s using %d cores.\n",
                length(unique_files), current_day, cores))

    # Determine which files to plot
    if (save_plot) {
      if (!is.null(indices)) {
        files_to_plot <- unique_files
      } else {
        n_plots <- ceiling(length(unique_files) * plot_percent / 100)
        files_to_plot <- sort(sample(unique_files, n_plots))
      }
    }

    process_file <- function(i) {
      tryCatch({
        should_plot <- save_plot && (i %in% files_to_plot)
        wavfile <- file.path(x$base_path,
                             day_metadata$day_post_hatch[i],
                             day_metadata$filename[i])
        plot_dir <- if (should_plot) {
          file.path(x$base_path, "plots", "template_matches",
                    paste0("day_", day_metadata$day_post_hatch[i]))
        } else NULL

        result <- detect_template.default(
          x = wavfile,
          template = template,
          cor.method = cor.method,
          save_plot = should_plot,
          plot_dir = plot_dir,
          proximity_window = proximity_window,
          ...
        )

        if (!is.null(result)) {
          result <- result |>
            dplyr::mutate(
              day_post_hatch = day_metadata$day_post_hatch[i],
              label = day_metadata$label[i],
              .after = filename
            )
        }
        return(result)
      }, error = function(e) {
        warning(sprintf("Error processing file %s: %s",
                        day_metadata$filename[i], e$message))
        return(NULL)
      })
    }

    # Parallel processing
    day_results <- parallel_apply(unique_files, process_file, cores)

    valid_detections <- day_results[!sapply(day_results, is.null)]
    if (length(valid_detections) > 0) {
      day_detections <- do.call(rbind, valid_detections)
      all_results[[as.character(current_day)]] <- day_detections
      cat(sprintf("\nProcessed files in day %s. Total detections: %d\n",
                  current_day, nrow(day_detections)))
    } else {
      cat(sprintf("\nNo detections found in day %s.\n", current_day))
    }
  }

  if (length(all_results) > 0) {
    final_results <- do.call(rbind, all_results)
    row.names(final_results) <- NULL
    x$templates$template_matches[[template_name]] <- final_results
    x$misc$last_modified <- Sys.time()
    message(sprintf("\nTotal detections across all days: %d", nrow(final_results)))
    message(sprintf("Access detection results via: Sap_object$templates$template_matches[[\"%s\"]]",
                    template_name))
  } else {
    warning(sprintf("No detections found for template '%s'", template_name))
  }

  invisible(x)
}

#' @keywords internal
find_peaks_with_proximity <- function(score.obj, proximity_window = NULL) {
  pks <- monitoR::findPeaks(score.obj = score.obj)

  # Process each template's detections
  for (tpl_name in names(pks@templates)) {
    dets <- pks@detections[[tpl_name]]

    # Skip if no detections
    if (is.null(dets) || nrow(dets) == 0) next

    # Sort by time and apply proximity filtering
    dets <- dets[order(dets$time), ]
    time_diff <- diff(dets$time)
    group_ids <- cumsum(c(TRUE, time_diff > proximity_window))

    filtered_dets <- dets %>%
      dplyr::mutate(group = group_ids) %>%
      dplyr::group_by(.data$group) %>%
      dplyr::slice_max(.data$score, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data$group)

    # Preserve original structure
    if (nrow(filtered_dets) > 0) {
      pks@detections[[tpl_name]] <- filtered_dets
    } else {
      pks@detections[[tpl_name]] <- data.frame()  # Empty but valid
    }
  }

  return(pks)
}
