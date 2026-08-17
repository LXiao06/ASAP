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
#'   start_time = 10,
#'   end_time = 20,
#'   clip_name = "song_clip"
#' )
#'
#' # Create multiple clips from SAP object
#' create_audio_clip(sap_object,
#'   indices = c(1, 2),
#'   start_time = c(10, 20),
#'   end_time = c(20, 30),
#'   clip_names = c("clip1", "clip2")
#' )
#'
#' # Create clip with millisecond units
#' create_audio_clip("song.wav",
#'   start_time = 10000,
#'   end_time = 20000,
#'   unit = "millisecond"
#' )
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
  if (verbose) message(sprintf("\n=== Starting Audio Clip Creation ===\n"))

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
    song_file <- file.path(
      x$base_path,
      x$metadata$day_post_hatch[index],
      x$metadata$filename[index]
    )

    # Create clip using default method with templates_path override
    clip_path <- file.path(templates_path, paste0(clip_names[i], ".wav"))

    # Create audio clip
    av::av_audio_convert(
      song_file,
      clip_path,
      start_time = if (unit == "millisecond") start_time[i] / 1000 else start_time[i],
      total_time = if (unit == "millisecond") {
        (end_time[i] - start_time[i]) / 1000
      } else {
        end_time[i] - start_time[i]
      }
    )

    cat(sprintf(
      "Audio clip '%s' is created \nOutput path: %s\n\n",
      clip_names[i],
      clip_path
    ))

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
#'   template_name = "template1",
#'   start_time = 1.0,
#'   end_time = 2.0,
#'   freq_min = 2,
#'   freq_max = 8
#' )
#'
#' # Create and save template
#' template <- create_template("song.wav",
#'   template_name = "template2",
#'   start_time = 1.0,
#'   end_time = 2.0,
#'   write_template = TRUE
#' )
#'
#' # Create template from SAP object
#' sap_obj <- create_template(sap_object,
#'   template_name = "template1",
#'   clip_name = "clip1",
#'   freq_min = 2,
#'   freq_max = 8
#' )
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
create_template.default <- function(x, # x is wav file path
                                    template_name, # mandatory
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
      score.cutoff = threshold
    )
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
      score.cutoff = threshold
    )
  }

  # Write template if requested
  if (write_template) {
    monitoR::writeCorTemplates(template, dir = dirname(x))
  }

  return(template)
}

#' @rdname create_template
#' @export
create_template.Sap <- function(x, # x is Sap object
                                template_name, # mandatory
                                clip_name, # mandatory, to identify which clip to use
                                start_time = NULL,
                                end_time = NULL,
                                freq_min = 0,
                                freq_max = 15,
                                threshold = 0.6,
                                write_template = FALSE,
                                verbose = TRUE,
                                ...) {
  if (verbose) message(sprintf("\n=== Starting Template Creation ===\n"))

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
#'        Can be a single numeric value (broadcast to all templates) or a vector matching \code{template_name} in length.
#'        Only the detection with the highest score within each window is retained.
#' @param template_name For SAP objects: Character vector of template name(s) to use
#' @param label For SAP objects: Vector of label(s) to process (character, numeric, or factor). NULL processes all labels.
#' @param indices For SAP objects: Numeric vector of indices to process
#' @param threshold For SAP objects: New threshold value(s) matching template_name
#' @param cores For SAP objects: Number of cores for parallel processing
#' @param plot_percent For SAP objects: Percentage of files to plot (default: 10)
#' @param verbose For SAP objects: Whether to print progress messages
#' @param use_preschedule For SAP objects: Whether to use pre-scheduling
#'        for parallel processing (default: FALSE). When TRUE, tasks are
#'        distributed to workers before processing starts, providing more
#'        consistent performance but may be slower if files have variable
#'        processing times. When FALSE (default), workers dynamically grab
#'        tasks from a queue, providing better load balancing but with
#'        slightly more overhead.
#' @param ... Additional arguments passed to specific methods (e.g. deprecated \code{day})
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
#'   \item Label-specific processing and non-numeric label support
#'   \item Multi-template matching: when multiple templates are supplied, they can be matched
#'         across all labels (when \code{label} is NULL or single) or paired 1-to-1 with \code{label}
#'   \item Optional threshold adjustment
#'   \item Progress tracking and reporting
#'   \item Selective plot generation
#'   \item Filtering of nearby detections with single or template-specific proximity windows
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
#'   template = template_obj,
#'   save_plot = TRUE
#' )
#'
#' # Detect template in SAP object for specific labels
#' sap_obj <- detect_template(sap_object,
#'   template_name = "template1",
#'   label = c("pre", "post"),
#'   threshold = 0.7,
#'   cores = 4
#' )
#'
#' # Multiple templates matched 1-to-1 with multiple labels
#' sap_obj <- detect_template(sap_object,
#'   template_name = c("tpl_pre", "tpl_post"),
#'   label = c("pre", "post")
#' )
#'
#' # Multiple templates with template-specific proximity windows
#' sap_obj <- detect_template(sap_object,
#'   template_name = c("tpl_pre", "tpl_post"),
#'   label = c("pre", "post"),
#'   proximity_window = c(0.5, 0.3)
#' )
#'
#' # Process specific indices with plots
#' sap_obj <- detect_template(sap_object,
#'   template_name = "template1",
#'   indices = 1:10,
#'   save_plot = TRUE
#' )
#'
#' # Filter nearby detections within 0.5 seconds
#' sap_obj <- detect_template(sap_object,
#'   template_name = "template1",
#'   proximity_window = 0.5
#' )
#'
#' # Use pre-scheduling for more consistent performance
#' sap_obj <- detect_template(sap_object,
#'   template_name = "template1",
#'   cores = 4,
#'   use_preschedule = TRUE
#' )
#' }
#'
#' @seealso \code{\link{create_template}} for creating templates
#'
#' @rdname detect_template
#' @export
detect_template <- function(x, ...) {
  UseMethod("detect_template")
}

# Render a monitoR detectionList without using monitoR's plot method.  Its
# current implementation formats the (double) axis tick positions with `%d`,
# which errors on recent versions of R.  This is intentionally private: it is
# only used when detect_template() is asked to save a plot.
.plot_template_detections <- function(x, flim = c(0, 12), t.each = 30) {
  survey <- x@survey
  survey_length <- length(survey@left) / survey@samp.rate
  n_plots <- ceiling(survey_length / t.each)
  t_start <- seq_len(n_plots) * t.each - t.each
  if (n_plots == 1L) t.each <- survey_length
  t_end <- t_start + t.each
  t_end[t_end > survey_length] <- survey_length
  t_start[n_plots] <- t_end[n_plots] - t.each

  amp <- x@survey.data[[1]]$amp
  t_bins <- x@survey.data[[1]]$t.bins
  frq_bins <- x@survey.data[[1]]$frq.bins
  template_names <- names(x@templates)
  colors <- rep(c("red", "blue", "green", "orange", "purple", "pink",
                  "darkgreen", "turquoise", "royalblue", "orchid4", "brown",
                  "salmon2"), length.out = length(template_names))
  names(colors) <- template_names
  score_limit <- c(0, max(vapply(x@scores, function(score) max(score$score), numeric(1))))

  old_par <- graphics::par(mar = c(1, 4, 1, 1), oma = c(6, 0, 0, 0), mfrow = c(2, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  for (i in seq_len(n_plots)) {
    times <- t_bins[t_bins >= t_start[i] & t_bins <= t_end[i]]
    amp_clip <- amp[, t_bins %in% times, drop = FALSE]
    graphics::image(
      x = times, y = frq_bins, z = t(amp_clip), ylim = flim,
      col = monitoR::gray.2(), xlab = "", ylab = "Frequency (kHz)",
      xaxt = "n", las = 1
    )

    for (j in template_names) {
      template <- x@templates[[j]]
      peaks <- x@detections[[j]]
      peaks <- peaks[peaks$time + template@duration >= t_start[i] &
                       peaks$time - template@duration <= t_end[i], , drop = FALSE]
      if (nrow(peaks) > 0L) {
        for (k in seq_len(nrow(peaks))) {
          x_left <- peaks$time[k] - template@duration / 2
          x_right <- peaks$time[k] + template@duration / 2
          graphics::polygon(
            x = c(x_left, x_left, x_right, x_right),
            y = c(template@frq.lim[1], template@frq.lim[2], template@frq.lim[2], template@frq.lim[1]),
            border = colors[j], lwd = 1
          )
        }
      }
    }

    graphics::plot(NULL, xlim = c(t_start[i], t_end[i]), ylim = score_limit,
                   xlab = "", ylab = "Score", type = "n", xaxs = "i", las = 1,
                   mgp = c(3, 1, 0))
    graphics::mtext("Time (s or min:sec)", 1, 2.5, outer = TRUE)
    x_axis <- graphics::par("xaxp")
    ticks <- seq(x_axis[1], x_axis[2], length.out = x_axis[3] + 1L)
    # `ticks` is double, even where every value is whole.  Convert its
    # components explicitly before applying the integer-only `%d` formatter.
    labels <- paste(
      sprintf("%02d", as.integer(floor(ticks / 60))),
      sprintf("%02d", as.integer(floor(ticks %% 60))), sep = ":"
    )
    graphics::axis(1, at = ticks, labels = labels, mgp = c(3, 1.9, 0))
    graphics::legend("topright", template_names, lty = 1, col = colors,
                     cex = 0.7)

    for (j in template_names) {
      template <- x@templates[[j]]
      score <- x@scores[[j]]
      peaks <- x@detections[[j]]
      score_clip <- score[score$time >= t_start[i] & score$time <= t_end[i], , drop = FALSE]
      peaks <- peaks[peaks$time + template@duration >= t_start[i] &
                       peaks$time - template@duration <= t_end[i], , drop = FALSE]
      graphics::lines(score_clip$time, score_clip$score, col = colors[j])
      graphics::abline(v = peaks$time, col = colors[j])
      if (is.vector(template@score.cutoff)) {
        graphics::abline(h = template@score.cutoff, lty = 2, col = colors[j])
      }
    }
  }
}

#' @rdname detect_template
#' @export
detect_template.default <- function(x, # x is wav file path
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
  # Safe on macOS: fork workers each get their own copy of the connection stack
  # Safe on Linux: PSOCK workers are isolated fresh R processes
  null_con <- file("/dev/null", open = "w")
  on.exit({
    sink(type = "output")
    sink(type = "message")
    close(null_con)
  })
  sink(null_con, type = "output")
  sink(null_con, type = "message")

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
    pks <- suppressMessages(monitoR::findPeaks(score.obj = scores))
  } else {
    pks <- suppressMessages(find_peaks_with_proximity(scores, proximity_window))
  }

  # Get detections
  detections <- monitoR::getDetections(pks, id = basename(x))

  # Drop detections with NA time values. These arise from flatline (all-zero)
  # regions in corrupted WAV files: the Pearson correlation denominator is 0
  # there, producing NaN scores whose peak positions are recorded as time = NA.
  # Real detections always have a well-defined time (frame_index / sample_rate)
  # and are never affected by this filter.
  n_na <- sum(is.na(detections$time))
  if (n_na > 0) {
    warning(sprintf(
      paste0(
        "%d detection(s) with NA time dropped from '%s'.\n",
        "  This usually indicates a flatline (silent/corrupted) region in the audio.\n",
        "  Check the file for stretches of zero-amplitude samples."
      ),
      n_na, basename(x)
    ))
    detections <- detections[!is.na(detections$time), ]
  }

  # Early return if no detections
  if (is.null(detections) || nrow(detections) == 0) {
    return(NULL)
  }

  # Process final output
  names(detections)[names(detections) == "id"] <- "filename"
  detections$date.time <- NULL

  # Generate plot with filtered detections
  if (save_plot) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    plot_file <- file.path(plot_dir, paste0(tools::file_path_sans_ext(basename(x)), ".png"))

    tryCatch(
      {
        saved_fd <- suppress_stderr()
        png(plot_file, width = 1200, height = 800, res = 150)
        restore_stderr(saved_fd)
        on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)

        suppressMessages(.plot_template_detections(pks))

        if (dev.cur() > 1) dev.off()
      },
      error = function(e) {
        warning("Error creating plot for ", basename(x), ": ", e$message)
        if (dev.cur() > 1) dev.off()
      }
    )
  }

  return(detections)
}

#' @rdname detect_template
#' @export
detect_template.Sap <- function(x, # x is SAP object
                                label = NULL,
                                indices = NULL,
                                template_name,
                                threshold = NULL,
                                cores = NULL,
                                cor.method = "pearson",
                                save_plot = FALSE,
                                plot_percent = 10,
                                verbose = TRUE,
                                proximity_window = NULL,
                                use_preschedule = FALSE,
                                ...) {
  if (verbose) message(sprintf("\n=== Starting Template Detection ==="))

  # Handle deprecation of `day` argument
  dots <- list(...)
  if (!is.null(dots$day)) {
    warning("Argument 'day' is deprecated; please use 'label' instead.")
    if (is.null(label)) {
      label <- dots$day
    }
  }

  # Validate inputs
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  if (missing(template_name) || is.null(template_name) || length(template_name) == 0) {
    stop("template_name must be provided")
  }

  missing_tpls <- setdiff(template_name, names(x$templates$template_list))
  if (length(missing_tpls) > 0) {
    stop("The following template(s) do not exist in the SAP object: ", paste(missing_tpls, collapse = ", "))
  }

  # Update thresholds if specified
  if (!is.null(threshold)) {
    if (length(threshold) == 1) {
      thresholds <- setNames(rep(threshold, length(unique(template_name))), unique(template_name))
    } else if (length(threshold) == length(template_name)) {
      thresholds <- setNames(threshold, template_name)
    } else {
      stop("length of threshold must be 1 or match length of template_name")
    }

    for (t_name in names(thresholds)) {
      t_val <- thresholds[[t_name]]
      tpl_obj <- x$templates$template_list[[t_name]]
      orig_thresh <- monitoR::templateCutoff(tpl_obj)[[1]]
      if (orig_thresh != t_val) {
        monitoR::templateCutoff(tpl_obj) <- setNames(t_val, t_name)
        x$templates$template_list[[t_name]] <- tpl_obj
        threshold_idx <- which(x$templates$template_info$template_name == t_name)
        if (length(threshold_idx) > 0) {
          x$templates$template_info$threshold[threshold_idx] <- t_val
        }
        if (verbose) {
          cat(sprintf("Template '%s' threshold updated from %.2f to %.2f\n", t_name, orig_thresh, t_val))
        }
      }
    }
  }

  # Parse proximity_window for single or template-specific values
  if (!is.null(proximity_window)) {
    if (!is.numeric(proximity_window)) {
      stop("proximity_window must be numeric")
    }
    if (length(proximity_window) == 1) {
      prox_windows <- setNames(rep(proximity_window, length(unique(template_name))), unique(template_name))
    } else if (length(proximity_window) == length(template_name)) {
      if (!is.null(names(proximity_window))) {
        prox_windows <- proximity_window
      } else {
        prox_windows <- setNames(proximity_window, template_name)
      }
    } else {
      stop(sprintf(
        "length of proximity_window (found %d) must be 1 or match length of template_name (%d)",
        length(proximity_window), length(template_name)
      ))
    }
  } else {
    prox_windows <- NULL
  }

  # Filter metadata based on label
  if (!is.null(label)) {
    process_metadata <- x$metadata[x$metadata$label %in% label, ]
    if (nrow(process_metadata) == 0) {
      stop("No files found for specified label(s)")
    }
    labels_to_process <- unique(process_metadata$label)
  } else {
    process_metadata <- x$metadata
    labels_to_process <- unique(process_metadata$label)
  }

  # Build mapping of label -> vector of template names
  label_templates <- list()
  if (is.null(label) || length(label) <= 1) {
    # Broad / Single-label dispatch: all templates applied to all targeted labels
    for (lbl in labels_to_process) {
      label_templates[[as.character(lbl)]] <- unique(template_name)
    }
  } else {
    # Multi-label specific dispatch: length(label) > 1 requires length(template_name) == length(label)
    if (length(template_name) != length(label)) {
      stop(sprintf(
        "When multiple labels are specified (length %d), template_name must have matching length (found length %d). You can repeat template names if needed.",
        length(label), length(template_name)
      ))
    }
    for (i in seq_along(label)) {
      lbl_str <- as.character(label[i])
      if (lbl_str %in% as.character(labels_to_process)) {
        label_templates[[lbl_str]] <- unique(c(label_templates[[lbl_str]], template_name[i]))
      }
    }
  }

  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }

  # On Linux, create one PSOCK cluster here and reuse it across all labels.
  psock_cl <- NULL
  if (Sys.info()["sysname"] != "Darwin" && cores > 1) {
    ensure_pkgs("parallel")
    psock_cl <- parallel::makeCluster(cores, type = "PSOCK")
    parallel::clusterEvalQ(psock_cl, {
      loadNamespace("ASAP")
      saved <- ASAP:::suppress_stderr()
      tmp <- tempfile(fileext = ".png")
      grDevices::png(tmp)
      plot(1, 1, main = "Init")
      grDevices::dev.off()
      unlink(tmp)
      ASAP:::restore_stderr(saved)
    })
    on.exit(parallel::stopCluster(psock_cl), add = TRUE)
  }

  # Create plot directories if save_plot is TRUE
  if (save_plot) {
    plots_dir <- file.path(x$base_path, "plots", "template_matches")
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    for (lbl in labels_to_process) {
      lbl_dir <- file.path(plots_dir, paste0("label_", sanitize_group_name(as.character(lbl))))
      dir.create(lbl_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  # Process each label group
  all_results <- list()

  for (current_label in labels_to_process) {
    label_str <- as.character(current_label)
    lbl_metadata <- process_metadata[process_metadata$label == current_label, ]
    if (!is.null(indices)) {
      valid_indices <- indices[indices <= nrow(lbl_metadata)]
      if (length(valid_indices) > 0) {
        lbl_metadata <- lbl_metadata[valid_indices, ]
      } else {
        if (verbose) cat(sprintf("\nNo valid indices for label %s.\n", label_str))
        next
      }
    }

    unique_files <- which(!duplicated(lbl_metadata$filename))
    tpls_for_label <- label_templates[[label_str]]
    if (is.null(tpls_for_label) || length(tpls_for_label) == 0) next

    if (verbose) {
      cat(sprintf(
        "\nProcessing %d files for label '%s' using template(s) [%s] on %d cores.\n",
        length(unique_files), label_str, paste(tpls_for_label, collapse = ", "), cores
      ))
    }

    # Prepare template object to pass to corMatch
    if (length(tpls_for_label) == 1) {
      template_to_use <- x$templates$template_list[[tpls_for_label]]
    } else {
      template_to_use <- monitoR::combineCorTemplates(
        unname(x$templates$template_list[tpls_for_label])
      )
    }

    # Subset proximity window to only the templates used for this label.
    # Passing a single scalar (for 1-template labels) or a named sub-vector
    # ensures find_peaks_with_proximity gets exactly the values it needs and
    # avoids confusion from extra entries for other labels' templates.
    label_prox_window <- if (!is.null(prox_windows)) {
      if (length(tpls_for_label) == 1) {
        unname(prox_windows[tpls_for_label])
      } else {
        prox_windows[tpls_for_label]
      }
    } else {
      NULL
    }

    # Determine which files to plot
    if (save_plot) {
      if (!is.null(indices)) {
        files_to_plot <- unique_files
      } else {
        n_plots <- ceiling(length(unique_files) * plot_percent / 100)
        files_to_plot <- sort(sample(unique_files, n_plots))
      }
    } else {
      files_to_plot <- integer(0)
    }

    # Explicitly capture all loop-varying values so the closure does not
    # accidentally reference a later iteration's bindings if R's lazy
    # evaluation or any future parallelism refactoring is involved.
    local_template    <- template_to_use
    local_prox        <- label_prox_window
    local_lbl_meta    <- lbl_metadata
    local_label_str   <- label_str
    local_files_plot  <- files_to_plot

    process_file <- function(i) {
      tryCatch(
        {
          should_plot <- save_plot && (i %in% local_files_plot)
          wavfile <- construct_wav_path(local_lbl_meta[i, ], wav_dir = x$base_path)
          plot_dir <- if (should_plot) {
            file.path(
              x$base_path, "plots", "template_matches",
              paste0("label_", sanitize_group_name(local_label_str))
            )
          } else {
            NULL
          }

          result <- detect_template.default(
            x = wavfile,
            template = local_template,
            cor.method = cor.method,
            save_plot = should_plot,
            plot_dir = plot_dir,
            proximity_window = local_prox,
            ...
          )

          if (!is.null(result)) {
            result <- result |>
              dplyr::mutate(
                label = local_lbl_meta$label[i],
                day_post_hatch = if ("day_post_hatch" %in% names(local_lbl_meta)) local_lbl_meta$day_post_hatch[i] else NA,
                .after = filename
              )
          }
          return(result)
        },
        error = function(e) {
          warning(sprintf(
            "Error processing file %s: %s",
            local_lbl_meta$filename[i], e$message
          ))
          return(NULL)
        }
      )
    }

    # Parallel processing
    lbl_results <- parallel_apply(unique_files, process_file, cores,
      use_preschedule = use_preschedule, cl = psock_cl
    )

    valid_detections <- lbl_results[!sapply(lbl_results, is.null)]
    if (length(valid_detections) > 0) {
      lbl_detections <- do.call(rbind, valid_detections)
      all_results[[label_str]] <- lbl_detections
      if (verbose) {
        cat(sprintf(
          "\nProcessed files for label '%s'. Total detections: %d\n",
          label_str, nrow(lbl_detections)
        ))
      }
    } else {
      if (verbose) cat(sprintf("\nNo detections found for label '%s'.\n", label_str))
    }
  }

  if (length(all_results) > 0) {
    final_results <- do.call(rbind, all_results)
    row.names(final_results) <- NULL

    for (t_name in unique(template_name)) {
      t_matches <- final_results[final_results$template == t_name, , drop = FALSE]
      row.names(t_matches) <- NULL
      x$templates$template_matches[[t_name]] <- t_matches
      if (verbose) {
        cat(sprintf("\nTotal detections for template '%s': %d\n", t_name, nrow(t_matches)))
        cat(sprintf(
          "Access detection results via: Sap_object$templates$template_matches[[\"%s\"]]\n",
          t_name
        ))
      }
    }
    x$misc$last_modified <- Sys.time()
  } else {
    for (t_name in unique(template_name)) {
      x$templates$template_matches[[t_name]] <- data.frame()
    }
    warning("No detections found for specified template(s)")
  }

  invisible(x)
}

#' @keywords internal
find_peaks_with_proximity <- function(score.obj, proximity_window = NULL) {
  pks <- monitoR::findPeaks(score.obj = score.obj)

  if (is.null(proximity_window)) {
    return(pks)
  }

  # Process each template's detections
  for (tpl_name in names(pks@templates)) {
    dets <- pks@detections[[tpl_name]]

    # Skip if no detections
    if (is.null(dets) || nrow(dets) == 0) next

    # Determine proximity window for this specific template
    win <- if (!is.null(names(proximity_window)) && tpl_name %in% names(proximity_window)) {
      proximity_window[[tpl_name]]
    } else if (length(proximity_window) == 1) {
      proximity_window[[1]]
    } else {
      NULL
    }

    if (is.null(win) || is.na(win) || win <= 0) next

    # Sort by time and apply proximity filtering using anchor-based grouping.
    # Each group is anchored at its first detection: a new group starts only
    # when the current peak is > proximity_window seconds from the *anchor*
    # (first peak of the group), not just the immediately preceding peak.
    # This prevents a long chain of closely-spaced peaks from being incorrectly
    # collapsed into one group even when they collectively span much longer than
    # proximity_window.
    dets <- dets[order(dets$time), ]
    n <- nrow(dets)
    group_ids <- integer(n)
    group_ids[1] <- 1L
    anchor_time <- dets$time[1]
    current_group <- 1L
    for (j in seq_len(n - 1) + 1L) {
      if (dets$time[j] - anchor_time > win) {
        current_group <- current_group + 1L
        anchor_time <- dets$time[j]
      }
      group_ids[j] <- current_group
    }

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
      pks@detections[[tpl_name]] <- data.frame() # Empty but valid
    }
  }

  return(pks)
}
