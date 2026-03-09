# Create motif/bout clips --------------------------------------------------------
# Update date : Mar. 6, 2026

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
#' @param name_prefix Prefix used for generated motif clip names.
#' @param amp_normalize Waveform amplitude normalization applied to exported clips:
#'   one of "none", "peak", or "rms" (default: "none")
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
                                       name_prefix = "motif",
                                       amp_normalize = c("none", "peak", "rms"),
                                       cores = NULL,
                                       overwrite = TRUE,
                                       write_metadata = TRUE,
                                       verbose = TRUE,
                                       ...) {
  if (!is.data.frame(x)) {
    stop("Input must be a data frame for default method")
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

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  motifs <- .subset_clip_rows(
    x,
    indices = indices,
    n_clips = n_motifs,
    seed = seed
  )
  motifs <- .prepare_clip_export_rows(motifs)

  export_meta <- .export_clip_rows(
    motifs = motifs,
    wav_dir = wav_dir,
    output_format = output_format,
    output_dir = output_dir,
    output_subdir = "motifs",
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    amp_normalize = amp_normalize,
    cores = cores,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose
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
                                   name_prefix = "motif",
                                   amp_normalize = c("none", "peak", "rms"),
                                   cores = NULL,
                                   overwrite = TRUE,
                                   write_metadata = TRUE,
                                   verbose = TRUE,
                                   ...) {
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  if (!inherits(x$motifs, "segment") || nrow(x$motifs) == 0) {
    stop("No motifs found in SAP object. Run find_motif() first.")
  }

  amp_normalize <- .parse_amp_normalize(amp_normalize)

  if (is.null(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be provided as a single directory path")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  motifs <- .subset_clip_rows(
    x$motifs,
    indices = indices,
    n_clips = n_motifs,
    seed = seed
  )
  motifs <- .add_bird_id_from_metadata(motifs = motifs, metadata = x$metadata)
  motifs <- .prepare_clip_export_rows(motifs)

  export_meta <- .export_clip_rows(
    motifs = motifs,
    wav_dir = x$base_path,
    output_format = match.arg(output_format),
    output_dir = output_dir,
    output_subdir = "motifs",
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    amp_normalize = amp_normalize,
    cores = cores,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose
  )

  current_export <- list(
    creation_date = Sys.time(),
    output_format = match.arg(output_format),
    amp_normalize = amp_normalize,
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
    selected_indices = as.integer(motifs$.source_index),
    exported_indices = if (nrow(export_meta) > 0 && "source_index" %in% names(export_meta)) {
      as.integer(export_meta$source_index)
    } else {
      integer(0)
    },
    skipped_indices = if (nrow(export_meta) > 0 && "source_index" %in% names(export_meta)) {
      setdiff(as.integer(motifs$.source_index), as.integer(export_meta$source_index))
    } else {
      as.integer(motifs$.source_index)
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
#' @param name_prefix Prefix used for generated bout clip names.
#' @param amp_normalize Waveform amplitude normalization applied to exported clips:
#'   one of "none", "peak", or "rms" (default: "none").
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
                                      name_prefix = "bout",
                                      amp_normalize = c("none", "peak", "rms"),
                                      cores = NULL,
                                      overwrite = TRUE,
                                      write_metadata = TRUE,
                                      verbose = TRUE,
                                      ...) {
  if (!is.data.frame(x)) {
    stop("Input must be a data frame for default method")
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

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  bouts <- .subset_clip_rows(x, indices = indices, n_clips = n_bouts, seed = seed)
  bouts <- .prepare_clip_export_rows(bouts)

  export_meta <- .export_clip_rows(
    motifs = bouts,
    wav_dir = wav_dir,
    output_format = output_format,
    output_dir = output_dir,
    output_subdir = "bouts",
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    amp_normalize = amp_normalize,
    cores = cores,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose
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
                                  name_prefix = "bout",
                                  amp_normalize = c("none", "peak", "rms"),
                                  cores = NULL,
                                  overwrite = TRUE,
                                  write_metadata = TRUE,
                                  verbose = TRUE,
                                  ...) {
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  if (is.null(x$bouts) || nrow(x$bouts) == 0) {
    stop("No bouts found in SAP object. Run find_bout() first.")
  }

  amp_normalize <- .parse_amp_normalize(amp_normalize)

  if (is.null(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be provided as a single directory path")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  bouts <- .subset_clip_rows(
    x$bouts,
    indices = indices,
    n_clips = n_bouts,
    seed = seed
  )
  bouts <- .add_bird_id_from_metadata(motifs = bouts, metadata = x$metadata)
  bouts <- .prepare_clip_export_rows(bouts)

  export_meta <- .export_clip_rows(
    motifs = bouts,
    wav_dir = x$base_path,
    output_format = match.arg(output_format),
    output_dir = output_dir,
    output_subdir = "bouts",
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    amp_normalize = amp_normalize,
    cores = cores,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose
  )

  current_export <- list(
    creation_date = Sys.time(),
    output_format = match.arg(output_format),
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

.subset_clip_rows <- function(rows, indices = NULL, n_clips = NULL, seed = 222) {
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

.add_bird_id_from_metadata <- function(motifs, metadata) {
  if ("bird_id" %in% names(motifs) && all(!is.na(motifs$bird_id))) {
    return(motifs)
  }

  if (!is.data.frame(metadata) || !"bird_id" %in% names(metadata)) {
    motifs$bird_id <- "unknown_bird"
    return(motifs)
  }

  keys <- intersect(c("filename", "day_post_hatch", "label"), names(motifs))
  md_keys <- intersect(c("filename", "day_post_hatch", "label"), names(metadata))
  keys <- intersect(keys, md_keys)

  if (length(keys) == 0) {
    motifs$bird_id <- "unknown_bird"
    return(motifs)
  }

  md_unique <- unique(metadata[, c(keys, "bird_id"), drop = FALSE])
  merged <- merge(motifs, md_unique, by = keys, all.x = TRUE, sort = FALSE)

  if (!".source_index" %in% names(merged)) {
    merged$.source_index <- seq_len(nrow(merged))
  }
  merged <- merged[order(merged$.source_index), , drop = FALSE]
  merged$bird_id[is.na(merged$bird_id) | merged$bird_id == ""] <- "unknown_bird"
  merged
}

.prepare_clip_export_rows <- function(motifs) {
  if (!"day_post_hatch" %in% names(motifs)) {
    motifs$day_post_hatch <- "unknown_day"
  }

  if (!"label" %in% names(motifs)) {
    motifs$label <- NA_character_
  }

  if (!"bird_id" %in% names(motifs)) {
    motifs$bird_id <- "unknown_bird"
  }

  motifs$bird_id[is.na(motifs$bird_id) | motifs$bird_id == ""] <- "unknown_bird"
  motifs$day_post_hatch[is.na(motifs$day_post_hatch)] <- "unknown_day"
  motifs
}

.export_clip_rows <- function(motifs,
                              wav_dir,
                              output_format,
                              output_dir,
                              output_subdir,
                              hdf5_filename,
                              metadata_filename,
                              name_prefix,
                              amp_normalize,
                              cores,
                              overwrite,
                              write_metadata,
                              verbose) {
  jobs <- motifs
  jobs$day_summary <- as.character(jobs$day_post_hatch)
  jobs$day_summary[is.na(jobs$day_summary) | jobs$day_summary == ""] <- "unknown_day"
  jobs$bird_id_clean <- vapply(as.character(jobs$bird_id), .sanitize_group_name, FUN.VALUE = character(1))
  jobs$day_clean <- vapply(as.character(jobs$day_post_hatch), .sanitize_group_name, FUN.VALUE = character(1))
  group_key <- paste(jobs$bird_id_clean, jobs$day_clean, sep = "::")
  jobs$clip_seq <- stats::ave(seq_len(nrow(jobs)), group_key, FUN = seq_along)
  jobs$clip_id <- sprintf("%s_%03d", name_prefix, jobs$clip_seq)

  has_day <- !is.na(jobs$day_post_hatch) & jobs$day_post_hatch != "unknown_day"
  jobs$source_path <- ifelse(
    has_day,
    file.path(wav_dir, as.character(jobs$day_post_hatch), as.character(jobs$filename)),
    file.path(wav_dir, as.character(jobs$filename))
  )

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
      out_dir <- file.path(output_dir, output_subdir, as.character(row$bird_id_clean), as.character(row$day_clean))
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
      bird_group <- .h5_get_or_create_group(h5, res$bird_id)
      day_group <- .h5_get_or_create_group(bird_group, res$day)
      ds <- day_group$create_dataset(name = res$clip_id, robj = res$sample_vec)
      hdf5r::h5attr(ds, "sample_rate") <- as.integer(res$sample_rate)
      hdf5r::h5attr(ds, "source_file") <- as.character(res$metadata$filename[[1]])
      hdf5r::h5attr(ds, "start_time") <- as.numeric(res$metadata$start_time[[1]])
      hdf5r::h5attr(ds, "end_time") <- as.numeric(res$metadata$end_time[[1]])
      hdf5r::h5attr(ds, "label") <- as.character(res$metadata$label[[1]])
      hdf5r::h5attr(ds, "amp_normalize") <- as.character(amp_normalize)
      ds$close()
      sample_rates <- c(sample_rates, as.integer(res$sample_rate))
    }

    hdf5r::h5attr(h5, "creation_date") <- as.character(Sys.time())
    hdf5r::h5attr(h5, "n_clips") <- as.integer(sum(vapply(results, function(x) identical(x$status, "written"), logical(1))))
    hdf5r::h5attr(h5, "amp_normalize") <- as.character(amp_normalize)
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

.sanitize_group_name <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "") {
    return("unknown")
  }
  gsub("[^A-Za-z0-9._-]", "_", x)
}

.h5_get_or_create_group <- function(parent, group_name) {
  if (group_name %in% names(parent)) {
    return(parent[[group_name]])
  }
  parent$create_group(group_name)
}


.merge_motif_metadata_with_spectral <- function(metadata_df,
                                                spectral_df,
                                                output_file,
                                                time_digits = 6,
                                                verbose = TRUE) {
  # metadata_df must have all canonical columns (it is the left / authoritative side).
  required_meta <- c("filename", "start_time", "end_time", "duration", "day_post_hatch", "label")
  if (!all(required_meta %in% names(metadata_df))) {
    stop("metadata_df must contain: filename, start_time, end_time, duration, day_post_hatch, label")
  }

  # Coerce spectral_df to a plain data.frame to guard against tibble / matrix
  # edge cases where names() might behave unexpectedly after parallel rbind.
  spectral_df <- as.data.frame(spectral_df)

  # spectral_df needs only filename + time columns for key-building.
  # duration is intentionally excluded: spectral_analysis() rounds it to
  # 1 decimal while metadata_df stores it at full precision, causing key
  # mismatches. filename + start_time + end_time alone uniquely identify a row.
  required_spec <- c("filename", "start_time", "end_time")
  missing_spec <- setdiff(required_spec, names(spectral_df))
  if (length(missing_spec) > 0) {
    stop("spectral_df must contain: ", paste(missing_spec, collapse = ", "))
  }

  # Build the merge key using only filename, start_time, end_time.
  # day_post_hatch / label are omitted (they come from metadata_df's left-join)
  # and duration is omitted (precision mismatch between the two data frames).
  make_key <- function(df) {
    filename <- as.character(df$filename)
    start_key <- sprintf(paste0("%.", time_digits, "f"), round(as.numeric(df$start_time), time_digits))
    end_key <- sprintf(paste0("%.", time_digits, "f"), round(as.numeric(df$end_time), time_digits))
    paste(filename, start_key, end_key, sep = "||")
  }

  metadata_df$.merge_key <- make_key(metadata_df)
  spectral_df$.merge_key <- make_key(spectral_df)

  # Keep metadata columns as canonical output columns; add only spectral-specific columns.
  keep_spec <- setdiff(names(spectral_df), intersect(names(spectral_df), names(metadata_df)))
  keep_spec <- union(".merge_key", keep_spec)
  spectral_reduced <- spectral_df[, keep_spec, drop = FALSE]

  merged <- merge(
    metadata_df,
    spectral_reduced,
    by = ".merge_key",
    all.x = TRUE,
    sort = FALSE
  )
  merged$.merge_key <- NULL

  utils::write.csv(merged, output_file, row.names = FALSE)

  if (verbose) {
    message(sprintf("Merged CSV path: %s", normalizePath(output_file, mustWork = TRUE)))
  }

  merged
}
