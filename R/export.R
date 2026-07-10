# Create motif/bout clips --------------------------------------------------------
# Update date : Mar. 26, 2026

#' Create Motif Audio Clips
#'
#' @description
#' Exports motif-level clips detected by \code{find_motif()} either as WAV files
#' in a directory tree or as datasets in a single HDF5 file.
#'
#' @param x A SAP object or a motif data frame.
#' @param indices Optional numeric vector of motif row indices to export.
#'   If NULL, all motifs are exported.
#' @param n_motifs Optional integer. Number of motifs to randomly sample when
#'   \code{indices} is NULL. Sampling is applied per day/subdirectory
#'   (\code{day_post_hatch}) when available. If NULL, all motifs are processed.
#' @param seed Integer random seed used when \code{n_motifs} sampling
#'   is applied. In this repository, sampling is fixed to \code{222}.
#' @param output_format Output type: \code{"wav"} or \code{"hdf5"}.
#' @param output_dir Directory where output files are written.
#' @param wav_dir For data-frame method: base directory containing source WAV files.
#' @param hdf5_filename File name used when \code{output_format = "hdf5"}.
#' @param metadata_filename Name of metadata CSV written to \code{output_dir}.
#' @param name_prefix Prefix used for generated clip file names. Can be any
#'   string (e.g. \code{"motif"} or \code{"bout"}). If \code{NULL}, defaults to
#'   \code{"motif"}.
#' @param keep_source_file_name Logical. If \code{TRUE}, uses the stem of the
#'   originating WAV file combined with the \code{selec} column (or sequential
#'   index) as the clip identifier (e.g. \code{S237_42685_001.wav}). This option
#'   is especially useful for Scenario A exports where tracing a clip back to its
#'   original recording is important. Overrides \code{name_prefix}.
#' @param amp_normalize Waveform amplitude normalization applied to exported clips:
#'   one of "none", "peak", or "rms" (default: "none")
#' @param noise_reduction Logical. If TRUE, denoise source WAV files before
#'   extracting motif clips. Denoised sources are written under
#'   \code{output_dir/motifs/denoised_sources} and used for clip extraction.
#' @param cores Number of CPU cores used for clip processing.
#' @param overwrite Logical. Overwrite existing output file(s) when TRUE.
#'   Default is TRUE.
#' @param write_metadata Logical. Write metadata CSV when TRUE.
#' @param verbose Logical. Print one export summary per day and overall totals.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Output layout for \code{output_format = "wav"}:
#' \preformatted{
#' output_dir/motifs/{bird_id}/{day_post_hatch}/{name_prefix}_001.wav
#' }
#'
#' Output layout for \code{output_format = "hdf5"}:
#' \preformatted{
#' output_dir/{hdf5_filename}
#'   /{bird_id}/{day_post_hatch}/{name_prefix}_001
#' }
#'
#' If \code{bird_id} is not available in motif rows, the function attempts to
#' infer it from metadata. Missing values are stored as \code{"unknown_bird"}.
#'
#' @return
#' For SAP input: updated SAP object with export summary in
#' \code{x$misc$motif_clip_exports}. For data-frame input: metadata data frame.
#'
#' @examples
#' \dontrun{
#' sap <- sap |>
#'   find_motif(template_name = "syllable_d", pre_time = 0.7, lag_time = 0.5) |>
#'   create_motif_clips(indices = 1:50, output_format = "wav")
#'
#' sap <- create_motif_clips(
#'   sap,
#'   output_format = "hdf5",
#'   hdf5_filename = "motifs.h5"
#' )
#'
#' # Export RMS-normalized motif clips
#' sap <- create_motif_clips(
#'   sap,
#'   output_format = "wav",
#'   amp_normalize = "rms"
#' )
#' }
#'
#' @export
create_motif_clips <- function(x, ...) {
  UseMethod("create_motif_clips")
}

#' @rdname create_motif_clips
#' @export
create_motif_clips.default <- function(x,
                                       wav_dir,
                                       indices = NULL,
                                       n_motifs = NULL,
                                       seed = 222,
                                       output_format = c("wav", "hdf5"),
                                       output_dir = NULL,
                                       hdf5_filename = "motifs.h5",
                                       metadata_filename = "metadata.csv",
                                       name_prefix = NULL,
                                       keep_source_file_name = FALSE,
                                       amp_normalize = c("none", "peak", "rms"),
                                       noise_reduction = FALSE,
                                       cores = NULL,
                                       overwrite = TRUE,
                                       write_metadata = TRUE,
                                       verbose = TRUE,
                                       ...) {
  if (!is.data.frame(x)) {
    stop("Input must be a data frame for default method")
  }
  if (is.null(name_prefix)) {
    name_prefix <- "motif"
  }

  required_cols <- c("filename", "start_time", "end_time")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (nrow(x) == 0) {
    stop("Input data frame is empty")
  }

  wav_dir <- normalizePath(wav_dir, mustWork = TRUE)
  output_format <- match.arg(output_format)
  amp_normalize <- .parse_amp_normalize(amp_normalize)
  if (!is.logical(noise_reduction) || length(noise_reduction) != 1) {
    stop("noise_reduction must be TRUE or FALSE")
  }

  if (is.null(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be provided as a single directory path")
  }

  output_dir <- file.path(output_dir, "motifs")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  clips <- subset_clip_rows(
    x,
    indices = indices,
    n_clips = n_motifs,
    seed = seed
  )
  clips <- prepare_clip_export_rows(clips)

  export_meta <- export_clip_rows(
    clips = clips,
    wav_dir = wav_dir,
    output_format = output_format,
    output_dir = output_dir,
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    keep_source_file_name = keep_source_file_name,
    amp_normalize = amp_normalize,
    noise_reduction = noise_reduction,
    cores = cores,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose,
    ...
  )

  return(export_meta)
}

#' @rdname create_motif_clips
#' @export
create_motif_clips.Sap <- function(x,
                                   indices = NULL,
                                   n_motifs = NULL,
                                   seed = 222,
                                   output_format = c("wav", "hdf5"),
                                   output_dir = NULL,
                                   hdf5_filename = "motifs.h5",
                                   metadata_filename = "metadata.csv",
                                   name_prefix = NULL,
                                   keep_source_file_name = FALSE,
                                   amp_normalize = c("none", "peak", "rms"),
                                   noise_reduction = FALSE,
                                   cores = NULL,
                                   overwrite = TRUE,
                                   write_metadata = TRUE,
                                   verbose = TRUE,
                                   ...) {
  if (verbose) message("\n=== Starting Motif Clip Extraction ===")

  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  if (!inherits(x$motifs, "segment") || nrow(x$motifs) == 0) {
    stop("No motifs found in SAP object. Run find_motif() first.")
  }

  output_format <- match.arg(output_format)
  amp_normalize <- .parse_amp_normalize(amp_normalize)
  if (!is.logical(noise_reduction) || length(noise_reduction) != 1) {
    stop("noise_reduction must be TRUE or FALSE")
  }
  if (is.null(name_prefix)) {
    name_prefix <- "motif"
  }

  if (is.null(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be provided as a single directory path")
  }

  output_dir <- file.path(output_dir, "motifs")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  clips <- subset_clip_rows(
    x$motifs,
    indices = indices,
    n_clips = n_motifs,
    seed = seed
  )
  clips <- add_bird_id_from_metadata(clips = clips, metadata = x$metadata)
  clips <- prepare_clip_export_rows(clips)

  export_meta <- export_clip_rows(
    clips = clips,
    wav_dir = x$base_path,
    output_format = output_format,
    output_dir = output_dir,
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    keep_source_file_name = keep_source_file_name,
    amp_normalize = amp_normalize,
    noise_reduction = noise_reduction,
    cores = cores,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose,
    ...
  )

  current_export <- list(
    creation_date = Sys.time(),
    output_format = output_format,
    amp_normalize = amp_normalize,
    noise_reduction = isTRUE(noise_reduction),
    output_dir = normalizePath(output_dir, mustWork = TRUE),
    n_requested = if (!is.null(indices)) {
      length(indices)
    } else if (!is.null(n_motifs)) {
      if ("day_post_hatch" %in% names(x$motifs)) {
        day_vals <- as.character(x$motifs$day_post_hatch)
        day_vals[is.na(day_vals) | day_vals == ""] <- "unknown_day"
        sum(pmin(as.integer(n_motifs), as.integer(table(day_vals))))
      } else {
        min(as.integer(n_motifs), nrow(x$motifs))
      }
    } else {
      nrow(x$motifs)
    },
    n_written = nrow(export_meta),
    selected_indices = as.integer(clips$.source_index),
    exported_indices = if (nrow(export_meta) > 0 && "source_index" %in% names(export_meta)) {
      as.integer(export_meta$source_index)
    } else {
      integer(0)
    },
    skipped_indices = if (nrow(export_meta) > 0 && "source_index" %in% names(export_meta)) {
      setdiff(as.integer(clips$.source_index), as.integer(export_meta$source_index))
    } else {
      as.integer(clips$.source_index)
    },
    metadata_path = if (write_metadata) {
      file.path(normalizePath(output_dir, mustWork = TRUE), metadata_filename)
    } else {
      NA_character_
    },
    metadata_df = export_meta
  )

  if (isTRUE(overwrite)) {
    x$misc$motif_clip_exports <- list(current_export)
  } else {
    prev <- x$misc$motif_clip_exports
    if (is.null(prev) || !is.list(prev)) prev <- list()
    x$misc$motif_clip_exports <- c(prev, list(current_export))
  }
  x$misc$last_modified <- Sys.time()

  invisible(x)
}

# Create bout clips -----------------------------------------------------------

#' Create Bout Audio Clips
#'
#' @description
#' Exports bout-level clips detected by \code{find_bout()} either as WAV files
#' in a directory tree or as datasets in a single HDF5 file. Bouts typically
#' span multiple motifs and are longer than individual motif clips.
#'
#' @param x A SAP object or a bouts data frame.
#' @param indices Optional numeric vector of bout row indices to export.
#'   If NULL, all bouts are exported.
#' @param n_bouts Optional integer. Number of bouts to randomly sample when
#'   \code{indices} is NULL. Sampling is applied per day/subdirectory
#'   (\code{day_post_hatch}) when available. If NULL, all bouts are processed.
#' @param seed Integer random seed used when \code{n_bouts} sampling is applied.
#'   Default is \code{222}.
#' @param output_format Output type: \code{"wav"} or \code{"hdf5"}.
#' @param output_dir Directory where output files are written.
#' @param wav_dir For data-frame method: base directory containing source WAV files.
#' @param hdf5_filename File name used when \code{output_format = "hdf5"}.
#' @param metadata_filename Name of metadata CSV written to \code{output_dir}.
#' @param name_prefix Prefix used for generated clip file names. Can be any
#'   string (e.g. \code{"bout"}). If \code{NULL}, defaults to \code{"bout"}. If
#'   \code{keep_source_file_name} is \code{TRUE}, this parameter is ignored.
#' @param keep_source_file_name Logical. If \code{TRUE}, uses the stem of the
#'   originating WAV file combined with the \code{selec} column (or sequential
#'   index) as the clip identifier (e.g. \code{S237_42685_001.wav}). This option
#'   is especially useful for Scenario A exports where tracing a clip back to its
#'   original recording is important. Overrides \code{name_prefix}.
#' @param amp_normalize Waveform amplitude normalization applied to exported clips:
#'   one of "none", "peak", or "rms" (default: "none").
#' @param margin Numeric vector of length 2 giving the time margin in seconds to
#'   prepend before each bout start and append after each bout end,
#'   e.g. \code{c(1, 2)} adds 1 s before and 2 s after.
#'   Both values must be non-negative. Default \code{c(0, 0)} preserves the
#'   original boundaries. Bouts whose adjusted \code{start_time} would fall
#'   below 0, or whose adjusted \code{end_time} would exceed the source WAV
#'   duration, are dropped and reported in the console summary.
#' @param cores Number of CPU cores used for clip processing.
#' @param overwrite Logical. Overwrite existing output file(s) when TRUE.
#'   Default is TRUE.
#' @param write_metadata Logical. Write metadata CSV when TRUE.
#' @param verbose Logical. Print one export summary per day and overall totals.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Output layout for \code{output_format = "wav"}:
#' \preformatted{
#' output_dir/bouts/{bird_id}/{day_post_hatch}/{name_prefix}_001.wav
#' }
#'
#' Output layout for \code{output_format = "hdf5"}:
#' \preformatted{
#' output_dir/{hdf5_filename}
#'   /{bird_id}/{day_post_hatch}/{name_prefix}_001
#' }
#'
#' If \code{bird_id} is not available in bout rows, the function attempts to
#' infer it from metadata. Missing values are stored as \code{"unknown_bird"}.
#'
#' @return
#' For SAP input: updated SAP object with export summary in
#' \code{x$misc$bout_clip_exports}. For data-frame input: metadata data frame.
#'
#' @examples
#' \dontrun{
#' # Export up to 50 bouts per day as WAV files
#' sap <- sap |>
#'   find_bout(min_duration = 0.4, summary = TRUE) |>
#'   create_bout_clips(
#'     n_bouts = 50,
#'     output_format = "wav",
#'     output_dir = "exported_bouts"
#'   )
#'
#' # Export all bouts to a single HDF5 file
#' sap <- create_bout_clips(
#'   sap,
#'   output_format = "hdf5",
#'   output_dir    = "exported_bouts",
#'   hdf5_filename = "bouts.h5"
#' )
#' }
#'
#' @export
create_bout_clips <- function(x, ...) {
  UseMethod("create_bout_clips")
}

#' @rdname create_bout_clips
#' @export
create_bout_clips.default <- function(x,
                                      wav_dir,
                                      indices = NULL,
                                      n_bouts = NULL,
                                      seed = 222,
                                      output_format = c("wav", "hdf5"),
                                      output_dir = NULL,
                                      hdf5_filename = "bouts.h5",
                                      metadata_filename = "metadata.csv",
                                      name_prefix = NULL,
                                      keep_source_file_name = FALSE,
                                      amp_normalize = c("none", "peak", "rms"),
                                      margin = c(0, 0),
                                      cores = NULL,
                                      overwrite = TRUE,
                                      write_metadata = TRUE,
                                      verbose = TRUE,
                                      ...) {
  if (!is.data.frame(x)) {
    stop("Input must be a data frame for default method")
  }
  if (is.null(name_prefix)) {
    name_prefix <- "bout"
  }

  required_cols <- c("filename", "start_time", "end_time")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (nrow(x) == 0) {
    stop("Input data frame is empty")
  }

  wav_dir <- normalizePath(wav_dir, mustWork = TRUE)
  output_format <- match.arg(output_format)
  amp_normalize <- .parse_amp_normalize(amp_normalize)

  if (is.null(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be provided as a single directory path")
  }

  output_dir <- file.path(output_dir, "bouts")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  bouts <- subset_clip_rows(x, indices = indices, n_clips = n_bouts, seed = seed)
  bouts <- prepare_clip_export_rows(bouts)

  margin_result <- apply_clip_margin(bouts, wav_dir = wav_dir, margin = margin)
  bouts <- margin_result$valid
  n_violated <- margin_result$n_violated

  if (n_violated > 0 && verbose) {
    message(sprintf(
      "Note: %d bout(s) dropped - adjusted time window falls outside source WAV bounds.",
      n_violated
    ))
  }

  if (nrow(bouts) == 0) {
    warning("No valid bouts remain after applying margin. Returning empty metadata.")
    return(data.frame())
  }

  export_meta <- export_clip_rows(
    clips = bouts,
    wav_dir = wav_dir,
    output_format = output_format,
    output_dir = output_dir,
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    keep_source_file_name = keep_source_file_name,
    amp_normalize = amp_normalize,
    cores = cores,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose,
    ...
  )

  return(export_meta)
}

#' @rdname create_bout_clips
#' @export
create_bout_clips.Sap <- function(x,
                                  indices = NULL,
                                  n_bouts = NULL,
                                  seed = 222,
                                  output_format = c("wav", "hdf5"),
                                  output_dir = NULL,
                                  hdf5_filename = "bouts.h5",
                                  metadata_filename = "metadata.csv",
                                  name_prefix = NULL,
                                  keep_source_file_name = FALSE,
                                  amp_normalize = c("none", "peak", "rms"),
                                  margin = c(0, 0),
                                  cores = NULL,
                                  overwrite = TRUE,
                                  write_metadata = TRUE,
                                  verbose = TRUE,
                                  ...) {
  if (verbose) message("\n=== Starting Bout Clip Extraction ===")

  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  if (is.null(x$bouts) || nrow(x$bouts) == 0) {
    stop("No bouts found in SAP object. Run find_bout() first.")
  }

  output_format <- match.arg(output_format)
  amp_normalize <- .parse_amp_normalize(amp_normalize)
  if (is.null(name_prefix)) {
    name_prefix <- "bout"
  }

  if (is.null(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be provided as a single directory path")
  }

  output_dir <- file.path(output_dir, "bouts")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  bouts <- subset_clip_rows(
    x$bouts,
    indices = indices,
    n_clips = n_bouts,
    seed = seed
  )
  bouts <- add_bird_id_from_metadata(clips = bouts, metadata = x$metadata)
  bouts <- prepare_clip_export_rows(bouts)

  margin_result <- apply_clip_margin(bouts, wav_dir = x$base_path, margin = margin)
  bouts <- margin_result$valid
  n_violated <- margin_result$n_violated

  if (n_violated > 0 && verbose) {
    message(sprintf(
      "Note: %d bout(s) dropped - adjusted time window falls outside source WAV bounds.",
      n_violated
    ))
  }

  if (nrow(bouts) == 0) {
    warning("No valid bouts remain after applying margin. Returning SAP object unchanged.")
    invisible(x)
  }

  export_meta <- export_clip_rows(
    clips = bouts,
    wav_dir = x$base_path,
    output_format = output_format,
    output_dir = output_dir,
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    keep_source_file_name = keep_source_file_name,
    amp_normalize = amp_normalize,
    cores = cores,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose,
    ...
  )

  current_export <- list(
    creation_date = Sys.time(),
    output_format = output_format,
    amp_normalize = amp_normalize,
    output_dir = normalizePath(output_dir, mustWork = TRUE),
    n_requested = if (!is.null(indices)) {
      length(indices)
    } else if (!is.null(n_bouts)) {
      if ("day_post_hatch" %in% names(x$bouts)) {
        day_vals <- as.character(x$bouts$day_post_hatch)
        day_vals[is.na(day_vals) | day_vals == ""] <- "unknown_day"
        sum(pmin(as.integer(n_bouts), as.integer(table(day_vals))))
      } else {
        min(as.integer(n_bouts), nrow(x$bouts))
      }
    } else {
      nrow(x$bouts)
    },
    n_written = nrow(export_meta),
    selected_indices = as.integer(bouts$.source_index),
    exported_indices = if (nrow(export_meta) > 0 && "source_index" %in% names(export_meta)) {
      as.integer(export_meta$source_index)
    } else {
      integer(0)
    },
    skipped_indices = if (nrow(export_meta) > 0 && "source_index" %in% names(export_meta)) {
      setdiff(as.integer(bouts$.source_index), as.integer(export_meta$source_index))
    } else {
      as.integer(bouts$.source_index)
    },
    metadata_path = if (write_metadata) {
      file.path(normalizePath(output_dir, mustWork = TRUE), metadata_filename)
    } else {
      NA_character_
    },
    metadata_df = export_meta
  )

  if (isTRUE(overwrite)) {
    x$misc$bout_clip_exports <- list(current_export)
  } else {
    prev <- x$misc$bout_clip_exports
    if (is.null(prev) || !is.list(prev)) prev <- list()
    x$misc$bout_clip_exports <- c(prev, list(current_export))
  }
  x$misc$last_modified <- Sys.time()

  invisible(x)
}

# Internal helpers ------------------------------------------------------------

#' Apply Time Margin to Clip Rows and Validate Against WAV Duration
#'
#' @description
#' Expands \code{start_time} and \code{end_time} in a clip data frame by the
#' requested margin, then validates each clip against its source WAV file
#' duration (read via header only). Clips whose adjusted
#' \code{start_time < 0} or adjusted \code{end_time > wav_duration} are
#' silently dropped and counted.
#'
#' @param clips Prepared clip data frame (must have \code{filename},
#'   \code{start_time}, \code{end_time}, and optionally \code{day_post_hatch}).
#' @param wav_dir Root directory containing source WAV files.
#' @param margin Numeric vector of length 2: seconds to prepend before start
#'   and append after end. Both values must be non-negative.
#'
#' @return A list with:
#'   \describe{
#'     \item{valid}{Clip data frame with validated, adjusted times.}
#'     \item{n_violated}{Integer count of dropped clips.}
#'   }
#'
#' @keywords internal
apply_clip_margin <- function(clips, wav_dir, margin) {
  # --- Validate margin -------------------------------------------------------
  if (!is.numeric(margin) || length(margin) != 2L) {
    stop("margin must be a numeric vector of length 2.")
  }
  if (any(is.na(margin)) || any(margin < 0)) {
    stop("Both margin values must be non-negative numbers.")
  }

  # Apply margin
  clips$start_time <- clips$start_time - margin[1]
  clips$end_time <- clips$end_time + margin[2]

  # If zero margin, skip expensive WAV-header reads
  if (margin[1] == 0 && margin[2] == 0) {
    return(list(valid = clips, n_violated = 0L))
  }

  # Check start_time < 0 (no file read needed)
  neg_start <- clips$start_time < 0

  # Check end_time > wav_duration (WAV header read, per unique file) 
  subdir <- if ("subfolder" %in% names(clips)) {
    as.character(clips$subfolder)
  } else if ("day_post_hatch" %in% names(clips)) {
    as.character(clips$day_post_hatch)
  } else {
    rep(NA_character_, nrow(clips))
  }
  
  has_subdir <- !is.na(subdir) & subdir != "unknown_day" & subdir != ""

  source_paths <- ifelse(
    has_subdir,
    file.path(wav_dir, subdir, as.character(clips$filename)),
    file.path(wav_dir, as.character(clips$filename))
  )

  unique_paths <- unique(source_paths)
  dur_map <- vapply(unique_paths, function(p) {
    if (!file.exists(p)) {
      return(NA_real_)
    }
    hdr <- tryCatch(
      tuneR::readWave(p, header = TRUE),
      error = function(e) NULL
    )
    if (is.null(hdr)) {
      return(NA_real_)
    }
    hdr$samples / hdr$sample.rate
  }, FUN.VALUE = numeric(1L))
  names(dur_map) <- unique_paths

  wav_dur <- dur_map[source_paths]
  exceed_end <- !is.na(wav_dur) & clips$end_time > wav_dur

  violated <- neg_start | exceed_end
  list(
    valid      = clips[!violated, , drop = FALSE],
    n_violated = sum(violated)
  )
}

#' Subset Clip Rows by Index or Random Sample
#'
#' @description
#' Subsets a clip data frame using explicit row indices or by randomly sampling
#' up to \code{n_clips} rows per day. Attaches a \code{.source_index} column
#' recording each row's original position in the input.
#'
#' @param rows Data frame of clip rows (e.g. motifs or bouts).
#' @param indices Optional integer vector of explicit row indices.
#' @param n_clips Optional integer. Maximum clips to sample per day when
#'   \code{indices} is NULL.
#' @param seed Integer random seed (fixed to 222 in this repository).
#'
#' @return Subset of \code{rows} with an added \code{.source_index} column.
#'
#' @keywords internal
subset_clip_rows <- function(rows, indices = NULL, n_clips = NULL, seed = 222) {
  if (is.null(indices)) {
    if (!is.null(n_clips)) {
      if (!is.numeric(n_clips) || length(n_clips) != 1 || is.na(n_clips)) {
        stop("n_clips must be a single numeric value")
      }

      n_clips <- as.integer(n_clips)
      if (n_clips < 1) {
        stop("n_clips must be >= 1")
      }

      if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
        stop("seed must be a single numeric value")
      }
      if (as.integer(seed) != 222L) {
        warning("seed is fixed to 222 in this repository; using 222")
      }
      set.seed(222L)

      if (n_clips >= nrow(rows)) {
        warning(sprintf(
          "n_clips (%d) >= available rows (%d); using all rows",
          n_clips,
          nrow(rows)
        ))
        rows$.source_index <- seq_len(nrow(rows))
        return(rows)
      }

      if ("day_post_hatch" %in% names(rows)) {
        day_vals <- as.character(rows$day_post_hatch)
        day_vals[is.na(day_vals) | day_vals == ""] <- "unknown_day"

        idx_by_day <- split(seq_len(nrow(rows)), day_vals)
        sampled_idx <- unlist(lapply(idx_by_day, function(idx) {
          sample(idx, size = min(n_clips, length(idx)), replace = FALSE)
        }), use.names = FALSE)
      } else {
        sampled_idx <- sample(seq_len(nrow(rows)), size = n_clips, replace = FALSE)
      }

      result <- rows[sampled_idx, , drop = FALSE]
      result$.source_index <- sampled_idx
      return(result)
    }

    rows$.source_index <- seq_len(nrow(rows))
    return(rows)
  }

  if (!is.numeric(indices) || anyNA(indices)) {
    stop("indices must be a numeric vector without NA values")
  }

  if (any(indices < 1) || any(indices > nrow(rows))) {
    stop(sprintf("indices must be between 1 and %d", nrow(rows)))
  }

  if (length(unique(indices)) != length(indices)) {
    warning("Duplicate indices detected; duplicates will be exported multiple times")
  }

  result <- rows[indices, , drop = FALSE]
  result$.source_index <- indices
  return(result)
}

#' Add Bird ID from SAP Metadata
#'
#' @description
#' Joins the \code{bird_id} column from a SAP metadata table into a clip data
#' frame using shared key columns (\code{filename}, \code{day_post_hatch},
#' \code{label}). Rows that cannot be matched are assigned
#' \code{"unknown_bird"}.
#'
#' @param clips Data frame of clip rows.
#' @param metadata SAP metadata data frame containing a \code{bird_id} column.
#'
#' @return \code{clips} with a populated \code{bird_id} column.
#'
#' @keywords internal
add_bird_id_from_metadata <- function(clips, metadata) {
  if ("bird_id" %in% names(clips) && all(!is.na(clips$bird_id))) {
    return(clips)
  }

  if (!is.data.frame(metadata) || !"bird_id" %in% names(metadata)) {
    clips$bird_id <- "unknown_bird"
    return(clips)
  }

  keys <- intersect(c("filename", "day_post_hatch", "label"), names(clips))
  md_keys <- intersect(c("filename", "day_post_hatch", "label"), names(metadata))
  keys <- intersect(keys, md_keys)

  if (length(keys) == 0) {
    clips$bird_id <- "unknown_bird"
    return(clips)
  }

  md_unique <- unique(metadata[, c(keys, "bird_id"), drop = FALSE])
  merged <- merge(clips, md_unique, by = keys, all.x = TRUE, sort = FALSE)

  if (!".source_index" %in% names(merged)) {
    merged$.source_index <- seq_len(nrow(merged))
  }
  merged <- merged[order(merged$.source_index), , drop = FALSE]
  merged$bird_id[is.na(merged$bird_id) | merged$bird_id == ""] <- "unknown_bird"
  merged
}

#' Prepare Clip Data Frame for Export
#'
#' @description
#' Ensures required columns (\code{day_post_hatch}, \code{label},
#' \code{bird_id}) are present, substituting sensible defaults for missing or
#' NA values.
#'
#' @param clips Data frame of clip rows.
#'
#' @return Standardised clip data frame ready for \code{export_clip_rows}.
#'
#' @keywords internal
prepare_clip_export_rows <- function(clips) {
  if (!"day_post_hatch" %in% names(clips)) {
    clips$day_post_hatch <- "unknown_day"
  }

  if (!"label" %in% names(clips)) {
    clips$label <- NA_character_
  }

  if (!"bird_id" %in% names(clips)) {
    clips$bird_id <- "unknown_bird"
  }

  clips$bird_id[is.na(clips$bird_id) | clips$bird_id == ""] <- "unknown_bird"
  clips$day_post_hatch[is.na(clips$day_post_hatch)] <- "unknown_day"
  clips
}

#' Core Clip Extraction and Writing Engine
#'
#' @description
#' Iterates over a prepared clip data frame, extracts each clip from its source
#' WAV file, applies optional amplitude normalisation or denoising, and writes
#' the result to WAV files or a single HDF5 archive. Returns a metadata data
#' frame describing every successfully written clip.
#'
#' @param clips Prepared clip data frame (output of \code{prepare_clip_export_rows}).
#' @param wav_dir Root directory containing source WAV files.
#' @param output_format One of \code{"wav"} or \code{"hdf5"}.
#' @param output_dir Directory where clips are written.
#' @param hdf5_filename File name for the HDF5 archive.
#' @param metadata_filename File name for the metadata CSV.
#' @param name_prefix Prefix for generated clip file names.
#' @param keep_source_file_name Logical; use WAV stem + index as clip name.
#' @param amp_normalize Amplitude normalisation method (\code{"none"},
#'   \code{"peak"}, or \code{"rms"}).
#' @param noise_reduction Logical; denoise source files before extraction.
#' @param cores Number of parallel cores.
#' @param overwrite Logical; overwrite existing files.
#' @param write_metadata Logical; write a metadata CSV.
#' @param verbose Logical; print per-day progress messages.
#'
#' @return Data frame of export metadata (one row per successfully written clip).
#'
#' @keywords internal
export_clip_rows <- function(clips,
                             wav_dir,
                             output_format,
                             output_dir,
                             hdf5_filename,
                             metadata_filename,
                             name_prefix,
                             keep_source_file_name,
                             amp_normalize,
                             noise_reduction = FALSE,
                             cores,
                             overwrite,
                             write_metadata,
                             verbose,
                             ...) {
  jobs <- clips
  jobs$day_summary <- as.character(jobs$day_post_hatch)
  jobs$day_summary[is.na(jobs$day_summary) | jobs$day_summary == ""] <- "unknown_day"
  jobs$bird_id_clean <- vapply(as.character(jobs$bird_id), sanitize_group_name, FUN.VALUE = character(1))
  jobs$day_clean <- vapply(as.character(jobs$day_post_hatch), sanitize_group_name, FUN.VALUE = character(1))
  group_key <- paste(jobs$bird_id_clean, jobs$day_clean, sep = "::")
  jobs$clip_seq <- stats::ave(seq_len(nrow(jobs)), group_key, FUN = seq_along)
  if (keep_source_file_name) {
    # Use WAV stem + a per-file sequential index as the clip identifier.
    # Bouts already carry a 'selec' column (their index within the source file).
    # Motifs do not have 'selec', so we compute an equivalent: a sequential
    # counter that restarts from 1 for every unique source WAV filename.
    wav_stems <- tools::file_path_sans_ext(as.character(jobs$filename))
    if ("selec" %in% names(jobs)) {
      file_idx <- as.integer(jobs$selec)
    } else {
      file_idx <- as.integer(
        stats::ave(seq_len(nrow(jobs)), jobs$filename, FUN = seq_along)
      )
    }
    jobs$clip_id <- sprintf("%s_%03d", wav_stems, file_idx)
  } else {
    jobs$clip_id <- sprintf("%s_%03d", name_prefix, jobs$clip_seq)
  }

  subdir <- if ("subfolder" %in% names(jobs)) {
    as.character(jobs$subfolder)
  } else if ("day_post_hatch" %in% names(jobs)) {
    as.character(jobs$day_post_hatch)
  } else {
    rep(NA_character_, nrow(jobs))
  }
  
  has_subdir <- !is.na(subdir) & subdir != "unknown_day" & subdir != ""
  
  jobs$source_path <- ifelse(
    has_subdir,
    file.path(wav_dir, subdir, as.character(jobs$filename)),
    file.path(wav_dir, as.character(jobs$filename))
  )

  if (isTRUE(noise_reduction)) {
    denoise_root <- file.path(output_dir, "denoised_sources")
    if (!dir.exists(denoise_root)) dir.create(denoise_root, recursive = TRUE, showWarnings = FALSE)
    unique_sources <- unique(jobs[, c("source_path", "day_clean"), drop = FALSE])
    denoise_map <- character(0)

    if (verbose) {
      message(sprintf("Denoising %d source file(s) before clip export...", nrow(unique_sources)))
    }

    for (i in seq_len(nrow(unique_sources))) {
      source_path <- as.character(unique_sources$source_path[[i]])
      day_dir <- as.character(unique_sources$day_clean[[i]])
      out_dir <- file.path(denoise_root, day_dir)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      if (!file.exists(source_path)) {
        warning("Skipping missing source file for denoise: ", source_path)
        next
      }

      clean_path <- tryCatch(
        {
          denoise_args <- list(...)
          default_args <- list(
            x = source_path,
            output_dir = out_dir,
            overwrite = overwrite,
            wl = 512L,
            ovlp = 75L,
            plot = FALSE,
            verbose = FALSE
          )
          
          call_args <- utils::modifyList(default_args, denoise_args)
          # Ensure core parameters cannot be overwritten by ...
          call_args$x <- source_path
          call_args$output_dir <- out_dir
          call_args$overwrite <- overwrite

          do.call(denoise.default, call_args)
        },
        error = function(e) {
          stop("Denoise failed for: ", source_path, "\n", conditionMessage(e))
        }
      )
      denoise_map[source_path] <- clean_path
    }

    mapped <- denoise_map[jobs$source_path]
    mapped[is.na(mapped)] <- jobs$source_path[is.na(mapped)]
    jobs$source_path <- unname(mapped)
  }

  process_one <- function(i) {
    row <- jobs[i, , drop = FALSE]
    source_path <- as.character(row$source_path)
    duration <- as.numeric(row$end_time) - as.numeric(row$start_time)

    if (!file.exists(source_path)) {
      warning("Skipping missing source file: ", source_path)
      return(list(status = "missing", day_summary = as.character(row$day_summary), metadata = NULL))
    }
    if (!is.finite(duration) || duration <= 0) {
      warning(sprintf("Skipping clip index %s with non-positive duration", row$.source_index))
      return(list(status = "invalid", day_summary = as.character(row$day_summary), metadata = NULL))
    }

    if (output_format == "wav") {
      out_dir <- file.path(output_dir, as.character(row$bird_id_clean), as.character(row$day_clean))
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      out_path <- file.path(out_dir, paste0(as.character(row$clip_id), ".wav"))

      if (file.exists(out_path) && !overwrite) {
        warning("Skipping existing file: ", out_path)
        return(list(status = "exists", day_summary = as.character(row$day_summary), metadata = NULL))
      }

      if (amp_normalize == "none") {
        av::av_audio_convert(
          source_path,
          out_path,
          start_time = as.numeric(row$start_time),
          total_time = as.numeric(duration),
          verbose = FALSE
        )
      } else {
        wave <- tuneR::readWave(
          filename = source_path,
          from = as.numeric(row$start_time),
          to = as.numeric(row$end_time),
          units = "seconds"
        )
        wave <- .normalize_wave_amplitude(
          wv = wave,
          method = amp_normalize,
          target_rms = 0.1,
          eps = 1e-8
        )
        tuneR::writeWave(wave, filename = out_path)
      }

      metadata <- data.frame(
        source_index = as.integer(row$.source_index),
        bird_id = as.character(row$bird_id_clean),
        day_post_hatch = as.character(row$day_post_hatch),
        label = as.character(row$label),
        filename = as.character(row$filename),
        start_time = as.numeric(row$start_time),
        end_time = as.numeric(row$end_time),
        duration = as.numeric(duration),
        amp_normalize = amp_normalize,
        noise_reduction = isTRUE(noise_reduction),
        clip_id = as.character(row$clip_id),
        output_path = normalizePath(out_path, mustWork = TRUE),
        stringsAsFactors = FALSE
      )
      return(list(status = "written", day_summary = as.character(row$day_summary), metadata = metadata))
    }

    wave <- tuneR::readWave(
      filename = source_path,
      from = as.numeric(row$start_time),
      to = as.numeric(row$end_time),
      units = "seconds"
    )
    if (amp_normalize != "none") {
      wave <- .normalize_wave_amplitude(
        wv = wave,
        method = amp_normalize,
        target_rms = 0.1,
        eps = 1e-8
      )
    }

    sample_vec <- if (isTRUE(wave@stereo)) {
      (as.numeric(wave@left) + as.numeric(wave@right)) / 2
    } else {
      as.numeric(wave@left)
    }
    norm_factor <- max(1, 2^(as.numeric(wave@bit) - 1))
    sample_vec <- as.single(sample_vec / norm_factor)

    metadata <- data.frame(
      source_index = as.integer(row$.source_index),
      bird_id = as.character(row$bird_id_clean),
      day_post_hatch = as.character(row$day_post_hatch),
      label = as.character(row$label),
      filename = as.character(row$filename),
      start_time = as.numeric(row$start_time),
      end_time = as.numeric(row$end_time),
      duration = as.numeric(duration),
      amp_normalize = amp_normalize,
      noise_reduction = isTRUE(noise_reduction),
      clip_id = as.character(row$clip_id),
      output_path = paste0("/", as.character(row$bird_id_clean), "/", as.character(row$day_clean), "/", as.character(row$clip_id)),
      stringsAsFactors = FALSE
    )

    list(
      status = "written",
      day_summary = as.character(row$day_summary),
      metadata = metadata,
      sample_vec = sample_vec,
      sample_rate = as.integer(wave@samp.rate),
      bird_id = as.character(row$bird_id_clean),
      day = as.character(row$day_clean),
      clip_id = as.character(row$clip_id)
    )
  }

  results <- parallel_apply(seq_len(nrow(jobs)), process_one, cores = cores)

  sample_rates <- integer(0)
  if (output_format == "hdf5") {
    if (verbose) {
      message("Parallel clip processing completed. Starting serialized HDF5 writing (this may take a while for large exports)...")
    }

    ensure_pkgs("hdf5r")
    h5_path <- file.path(output_dir, hdf5_filename)
    if (file.exists(h5_path)) {
      if (!overwrite) {
        stop("HDF5 file already exists: ", h5_path, ". Set overwrite = TRUE to replace.")
      }
      file.remove(h5_path)
    }
    h5 <- hdf5r::H5File$new(filename = h5_path, mode = "w")
    on.exit(h5$close_all(), add = TRUE)

    for (res in results) {
      if (!identical(res$status, "written")) next
      bird_group <- h5_get_or_create_group(h5, res$bird_id)
      day_group <- h5_get_or_create_group(bird_group, res$day)
      ds <- day_group$create_dataset(name = res$clip_id, robj = res$sample_vec)
      hdf5r::h5attr(ds, "sample_rate") <- as.integer(res$sample_rate)
      hdf5r::h5attr(ds, "source_file") <- as.character(res$metadata$filename[[1]])
      hdf5r::h5attr(ds, "start_time") <- as.numeric(res$metadata$start_time[[1]])
      hdf5r::h5attr(ds, "end_time") <- as.numeric(res$metadata$end_time[[1]])
      hdf5r::h5attr(ds, "label") <- as.character(res$metadata$label[[1]])
      hdf5r::h5attr(ds, "amp_normalize") <- as.character(amp_normalize)
      hdf5r::h5attr(ds, "noise_reduction") <- as.logical(isTRUE(noise_reduction))
      ds$close()
      sample_rates <- c(sample_rates, as.integer(res$sample_rate))
    }

    hdf5r::h5attr(h5, "creation_date") <- as.character(Sys.time())
    hdf5r::h5attr(h5, "n_clips") <- as.integer(sum(vapply(results, function(x) identical(x$status, "written"), logical(1))))
    hdf5r::h5attr(h5, "amp_normalize") <- as.character(amp_normalize)
    hdf5r::h5attr(h5, "noise_reduction") <- as.character(isTRUE(noise_reduction))
    if (length(sample_rates) > 0 && length(unique(sample_rates)) == 1) {
      hdf5r::h5attr(h5, "sample_rate") <- as.integer(sample_rates[[1]])
    } else if (length(sample_rates) > 0) {
      hdf5r::h5attr(h5, "sample_rate") <- "mixed"
    }
  }

  metadata_rows <- lapply(results, function(x) x$metadata)
  export_meta <- do.call(rbind, metadata_rows[!vapply(metadata_rows, is.null, logical(1))])
  if (is.null(export_meta)) {
    warning("No clips were exported")
    export_meta <- data.frame()
  }

  if (write_metadata && nrow(export_meta) > 0) {
    metadata_path <- file.path(output_dir, metadata_filename)
    utils::write.csv(export_meta, metadata_path, row.names = FALSE)
  }

  if (verbose) {
    day_stats <- data.frame(
      day = jobs$day_summary,
      status = vapply(results, function(x) x$status, character(1)),
      stringsAsFactors = FALSE
    )
    day_keys <- sort(unique(day_stats$day))
    for (day_key in day_keys) {
      d <- day_stats[day_stats$day == day_key, , drop = FALSE]
      message(sprintf(
        "Day %s: %d written / %d total (skipped: %d missing, %d invalid, %d existing)",
        day_key,
        sum(d$status == "written"),
        nrow(d),
        sum(d$status == "missing"),
        sum(d$status == "invalid"),
        sum(d$status == "exists")
      ))
    }
    message(sprintf("Exported %d clips to %s format (%s normalization)", nrow(export_meta), output_format, amp_normalize))
    message(sprintf("Output directory: %s", output_dir))
  }

  export_meta
}

#' Sanitize a String for Use as a Group Name
#'
#' @description
#' Replaces any character that is not alphanumeric, a dot, underscore, or
#' hyphen with an underscore. Returns \code{"unknown"} for empty or NA input.
#'
#' @param x Character scalar to sanitize.
#'
#' @return Sanitized character string.
#'
#' @keywords internal
sanitize_group_name <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "") {
    return("unknown")
  }
  gsub("[^A-Za-z0-9._-]", "_", x)
}

#' Get or Create an HDF5 Group
#'
#' @description
#' Returns an existing group from an HDF5 parent object, or creates it if it
#' does not already exist.
#'
#' @param parent An \code{H5File} or \code{H5Group} object.
#' @param group_name Character name of the group to retrieve or create.
#'
#' @return The requested HDF5 group object.
#'
#' @keywords internal
h5_get_or_create_group <- function(parent, group_name) {
  if (group_name %in% names(parent)) {
    return(parent[[group_name]])
  }
  parent$create_group(group_name)
}


# Export Feature Data Frame to CSV --------------------------------------------
# Update date: July 9, 2026

#' Export a Stored Feature Data Frame to CSV
#'
#' @description
#' Extracts a pre-computed feature data frame stored in
#' \code{x$features[[feature_type]][[slot_name]]} and writes it to a CSV file.
#' Works universally with any feature slot produced by ASAP analysis functions,
#' including spectral features, trajectory similarity/variability/maturation
#' scores, and any other data frame or list-of-data-frames stored in the
#' \code{features} slot of a SAP object.
#'
#' When the stored object is a plain data frame
#' (e.g., \code{spectral_feature}, \code{maturation_scores}), it is exported
#' directly.
#'
#' When the stored object is a list (e.g., \code{trajectory_similarity},
#' \code{trajectory_dispersion}, \code{trajectory_path_deviation}), use the
#' \code{sub_table} argument to select which element of the list to export
#' (e.g., \code{sub_table = "similarity"}, \code{"dispersion"}, \code{"deviation"}).
#'
#' @param x A SAP object containing the features to export
#' @param feature_type Character. The top-level feature key inside
#'   \code{x$features}, e.g. \code{"motif"}, \code{"syllable"}, \code{"bout"},
#'   or \code{"segment"}
#' @param slot_name Character. The sub-slot key inside
#'   \code{x$features[[feature_type]]}, e.g. \code{"spectral_feature"},
#'   \code{"maturation_scores"}, \code{"trajectory_similarity"},
#'   \code{"trajectory_dispersion"}, or \code{"trajectory_path_deviation"}
#' @param sub_table Character or \code{NULL}. When the stored object is a
#'   list (trajectory results), the name of the list element to extract
#'   as a data frame (default: \code{NULL})
#' @param merge_segments Logical. If \code{TRUE} (default) and the extracted
#'   data frame contains \code{filename}, \code{start_time}, and
#'   \code{end_time} columns, it is left-joined against the canonical segment
#'   metadata stored in \code{x[[segment_type]]} to enrich the export with
#'   per-segment columns. When those columns are absent, the merge is skipped
#' @param segment_type Character or \code{NULL}. The segment slot to use for
#'   metadata merging (one of \code{"motifs"}, \code{"syllables"},
#'   \code{"bouts"}, \code{"segments"}). If \code{NULL} (default), derived
#'   automatically from \code{feature_type} by appending \code{"s"}
#' @param time_match_digits Integer. Number of decimal places used when
#'   building the merge key (default: \code{3})
#' @param balance_labels Logical. If \code{TRUE}, labels whose row count falls
#'   below threshold are dropped before writing the CSV (default: \code{FALSE})
#' @param balance_threshold Numeric in \code{(0, 1]}. Fraction of the median
#'   label count used as the minimum retention threshold when
#'   \code{balance_labels = TRUE} (default: \code{0.5})
#' @param exclude_anomalies Logical. If \code{TRUE}, excludes rows belonging to
#'   anomalous labels stored under \code{x$features[[feature_type]][["anomalous_labels"]]}
#'   if available (default: \code{TRUE})
#' @param csv_filename Character. Output file name (required)
#' @param output_dir Character. Directory where the CSV is written (default: \code{"."})
#' @param overwrite Logical. Whether to overwrite an existing file (default:
#'   \code{FALSE})
#' @param verbose Logical. Print progress and summary messages (default:
#'   \code{TRUE})
#'
#' @details
#' When the stored object is a list, the following sub-tables are typical choices
#' for \code{sub_table}:
#' \itemize{
#'   \item \code{trajectory_similarity}: \code{"similarity"}, \code{"summary"}
#'   \item \code{trajectory_dispersion}: \code{"dispersion"}, \code{"pairwise"}, \code{"path_length"}, \code{"summary"}
#'   \item \code{trajectory_path_deviation}: \code{"deviation"}, \code{"deviation_stats"}, \code{"summary"}
#' }
#'
#' @return The exported data frame, invisibly.
#'
#' @examples
#' \dontrun{
#' # Export spectral features stored from a previous analyze_spectral() run
#' export_feature_csv(sap,
#'   feature_type = "motif",
#'   slot_name    = "spectral_feature",
#'   csv_filename = "spectral_features.csv"
#' )
#'
#' # Export trajectory similarity scores
#' export_feature_csv(sap,
#'   feature_type   = "motif",
#'   slot_name      = "trajectory_similarity",
#'   sub_table      = "similarity",
#'   merge_segments = FALSE,
#'   csv_filename   = "traj_similarity.csv"
#' )
#' }
#'
#' @export
export_feature_csv <- function(x,
                               feature_type,
                               slot_name,
                               sub_table         = NULL,
                               merge_segments    = TRUE,
                               segment_type      = NULL,
                               time_match_digits = 3,
                               balance_labels    = FALSE,
                               balance_threshold = 0.5,
                               exclude_anomalies = TRUE,
                               csv_filename      = NULL,
                               output_dir        = ".",
                               overwrite         = FALSE,
                               verbose           = TRUE) {
  # Validate inputs
  if (!inherits(x, "Sap")) stop("'x' must be a Sap object")

  if (is.null(csv_filename) || !nzchar(csv_filename)) {
    stop("'csv_filename' is required")
  }

  if (!is.character(feature_type) || length(feature_type) != 1 || !nzchar(feature_type)) {
    stop("'feature_type' must be a non-empty character scalar (e.g. \"motif\")")
  }

  if (!is.character(slot_name) || length(slot_name) != 1 || !nzchar(slot_name)) {
    stop("'slot_name' must be a non-empty character scalar (e.g. \"spectral_feature\")")
  }

  if (!is.numeric(balance_threshold) || balance_threshold <= 0 || balance_threshold > 1) {
    stop("'balance_threshold' must be a number in (0, 1]")
  }

  # Extract stored object
  feature_list <- x$features[[feature_type]]
  if (is.null(feature_list)) {
    avail <- names(x$features)
    stop(sprintf(
      "No features found for feature_type = \"%s\".\nAvailable feature types: %s",
      feature_type,
      if (length(avail)) paste(avail, collapse = ", ") else "<none>"
    ))
  }

  stored <- feature_list[[slot_name]]
  if (is.null(stored)) {
    avail <- names(feature_list)
    stop(sprintf(
      paste0(
        "Slot \"%s\" not found in x$features$%s.\n",
        "Available slots: %s\n",
        "Run the corresponding analysis function first."
      ),
      slot_name,
      feature_type,
      if (length(avail)) paste(avail, collapse = ", ") else "<none>"
    ))
  }

  # Resolve list results via sub_table
  if (is.list(stored) && !is.data.frame(stored)) {
    if (is.null(sub_table)) {
      df_elements <- names(stored)[vapply(stored, is.data.frame, logical(1))]
      stop(sprintf(
        paste0(
          "x$features$%s$%s is a list, not a plain data frame.\n",
          "Specify 'sub_table' to select which element to export.\n",
          "Data frame elements available: %s"
        ),
        feature_type, slot_name,
        if (length(df_elements)) paste(df_elements, collapse = ", ") else "<none>"
      ))
    }

    if (!sub_table %in% names(stored)) {
      df_elements <- names(stored)[vapply(stored, is.data.frame, logical(1))]
      stop(sprintf(
        paste0(
          "sub_table = \"%s\" not found in x$features$%s$%s.\n",
          "Data frame elements available: %s"
        ),
        sub_table, feature_type, slot_name,
        if (length(df_elements)) paste(df_elements, collapse = ", ") else "<none>"
      ))
    }

    df <- stored[[sub_table]]

    if (!is.data.frame(df)) {
      stop(sprintf(
        "x$features$%s$%s$%s is not a data frame (class: %s)",
        feature_type, slot_name, sub_table,
        paste(class(df), collapse = "/")
      ))
    }

    if (verbose) {
      message(sprintf(
        "Extracted sub_table = \"%s\" from x$features$%s$%s  [%d rows x %d cols]",
        sub_table, feature_type, slot_name, nrow(df), ncol(df)
      ))
    }

  } else {
    # Plain data frame (e.g., spectral_feature, maturation_scores)
    df <- as.data.frame(stored)

    if (verbose) {
      message(sprintf(
        "Extracted x$features$%s$%s  [%d rows x %d cols]",
        feature_type, slot_name, nrow(df), ncol(df)
      ))
    }
  }

  # Merge with segment metadata
  merge_key_cols <- c("filename", "start_time", "end_time")
  has_merge_keys <- all(merge_key_cols %in% names(df))

  if (merge_segments) {
    if (has_merge_keys) {
      # Derive segment type if not provided
      seg_type <- if (!is.null(segment_type)) {
        segment_type
      } else {
        paste0(feature_type, "s")
      }

      segs <- x[[seg_type]]

      if (is.null(segs) || !is.data.frame(segs) || nrow(segs) == 0) {
        if (verbose) {
          message(sprintf(
            "Segment slot x$%s is empty or missing; skipping metadata merge",
            seg_type
          ))
        }
      } else {
        if (verbose) {
          message(sprintf(
            "Merging with segment metadata from x$%s  [%d rows]...",
            seg_type, nrow(segs)
          ))
        }

        df <- merge_feature_with_segments(
          segs_df       = segs,
          feature_df    = df,
          time_digits   = time_match_digits
        )

        if (verbose) {
          message(sprintf("  -> %d rows after merge", nrow(df)))
        }
      }
    } else {
      if (verbose) {
        message(paste(
          "merge_segments = TRUE but the data frame lacks",
          "filename/start_time/end_time columns; skipping metadata merge",
          "(This is normal for trajectory summary tables keyed by label + rendition.)"
        ))
      }
    }
  }

  # Drop under-represented labels
  if (balance_labels) {
    if (!"label" %in% names(df)) {
      warning(paste(
        "'balance_labels = TRUE' requested but no 'label' column found",
        "Skipping label balancing"
      ))
    } else {
      df <- balance_feature_labels(
        df        = df,
        threshold = balance_threshold,
        verbose   = verbose
      )
    }
  }

  # Exclude anomalous labels
  if (exclude_anomalies) {
    anomalous_labels <- x$features[[feature_type]][["anomalous_labels"]]
    if (!is.null(anomalous_labels) && length(anomalous_labels) > 0) {
      if ("label" %in% names(df)) {
        before_count <- nrow(df)
        df <- df[!df$label %in% anomalous_labels, , drop = FALSE]
        if (verbose) {
          message(sprintf(
            "Excluded %d rows belonging to anomalous labels: %s",
            before_count - nrow(df),
            paste(anomalous_labels, collapse = ", ")
          ))
        }
      } else {
        warning("'exclude_anomalies = TRUE' requested but no 'label' column found in feature data")
      }
    }
  }

  # Write CSV
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  out_path <- file.path(output_dir, csv_filename)

  if (file.exists(out_path) && !overwrite) {
    stop(sprintf(
      "File already exists: %s\nUse overwrite = TRUE to replace it",
      normalizePath(out_path, mustWork = FALSE)
    ))
  }

  utils::write.csv(df, out_path, row.names = FALSE)

  if (verbose) {
    message(sprintf(
      "\nExported %d rows x %d columns -> %s",
      nrow(df), ncol(df),
      normalizePath(out_path, mustWork = TRUE)
    ))
  }

  invisible(df)
}

#' Merge Feature Data Frame with Canonical Segment Metadata
#'
#' @description
#' Internal helper used by \code{export_feature_csv.Sap()}.
#' Left-joins a feature data frame onto the canonical segment data frame.
#'
#' @param segs_df    Data frame of canonical segments
#' @param feature_df Data frame of computed features
#' @param time_digits Integer decimal places for time-key formatting
#'
#' @keywords internal
#' @noRd
merge_feature_with_segments <- function(segs_df, feature_df, time_digits = 3) {
  # Build composite string key
  make_key <- function(df) {
    fn  <- as.character(df$filename)
    st  <- sprintf(paste0("%.", time_digits, "f"), round(as.numeric(df$start_time), time_digits))
    et  <- sprintf(paste0("%.", time_digits, "f"), round(as.numeric(df$end_time),   time_digits))
    paste(fn, st, et, sep = "||")
  }

  segs_df$.merge_key    <- make_key(segs_df)
  feature_df$.merge_key <- make_key(feature_df)

  # Keep only new feature columns
  extra_cols  <- setdiff(names(feature_df), c(names(segs_df), ".merge_key"))
  feature_slim <- feature_df[, c(".merge_key", extra_cols), drop = FALSE]

  merged <- merge(
    segs_df,
    feature_slim,
    by      = ".merge_key",
    all.x   = TRUE,
    sort    = FALSE
  )
  merged$.merge_key <- NULL

  merged
}


#' Drop Under-Represented Labels for ML Balance
#'
#' @description
#' Internal helper used by \code{export_feature_csv.Sap()}.
#' Identifies labels whose row count is below threshold and removes them.
#'
#' @param df        Data frame with a label column
#' @param threshold Numeric fraction of median count used as minimum retention threshold
#' @param verbose   Logical print summary table
#'
#' @keywords internal
#' @noRd
balance_feature_labels <- function(df, threshold = 0.5, verbose = TRUE) {
  counts     <- table(df$label)
  med_count  <- stats::median(as.integer(counts))
  min_count  <- ceiling(threshold * med_count)

  keep_labels <- names(counts)[counts >= min_count]
  drop_labels <- names(counts)[counts < min_count]

  if (verbose) {
    cat("\n== Label balance summary ==========================================\n")
    cat(sprintf(
      "  Median label count : %d\n  Retention threshold: >= %d rows (%.0f%% of median)\n\n",
      as.integer(med_count), as.integer(min_count), threshold * 100
    ))

    # Print table header
    cat(sprintf("  %-20s  %8s  %10s\n", "Label", "N rows", "Status"))
    cat(sprintf("  %-20s  %8s  %10s\n", "-----", "------", "------"))

    all_labels <- sort(names(counts))
    for (lbl in all_labels) {
      status <- if (lbl %in% keep_labels) "KEEP" else "DROP (too few)"
      cat(sprintf("  %-20s  %8d  %10s\n", lbl, counts[[lbl]], status))
    }
    cat("\n")
  }

  if (length(drop_labels) == 0) {
    if (verbose) message("All labels meet the threshold; no labels dropped")
    return(df)
  }

  if (verbose) {
    message(sprintf(
      "Dropped %d label(s): %s",
      length(drop_labels),
      paste(drop_labels, collapse = ", ")
    ))
  }

  df_balanced <- df[df$label %in% keep_labels, , drop = FALSE]

  if (verbose) {
    message(sprintf(
      "Retained %d labels: %s  [%d rows total]",
      length(keep_labels),
      paste(sort(keep_labels), collapse = ", "),
      nrow(df_balanced)
    ))
  }

  df_balanced
}

