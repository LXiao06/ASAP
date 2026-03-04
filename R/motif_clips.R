# Create motif clips --------------------------------------------------------
# Update date : Mar. 2, 2026

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

  if (is.null(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be provided as a single directory path")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  motifs <- .subset_motif_rows(
    x,
    indices = indices,
    n_motifs = n_motifs,
    seed = seed
  )
  motifs <- .prepare_motif_export_rows(motifs)

  export_meta <- .export_motif_rows(
    motifs = motifs,
    wav_dir = wav_dir,
    output_format = output_format,
    output_dir = output_dir,
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
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

  if (is.null(output_dir) || !is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be provided as a single directory path")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  motifs <- .subset_motif_rows(
    x$motifs,
    indices = indices,
    n_motifs = n_motifs,
    seed = seed
  )
  motifs <- .add_bird_id_from_metadata(motifs = motifs, metadata = x$metadata)
  motifs <- .prepare_motif_export_rows(motifs)

  export_meta <- .export_motif_rows(
    motifs = motifs,
    wav_dir = x$base_path,
    output_format = match.arg(output_format),
    output_dir = output_dir,
    hdf5_filename = hdf5_filename,
    metadata_filename = metadata_filename,
    name_prefix = name_prefix,
    overwrite = overwrite,
    write_metadata = write_metadata,
    verbose = verbose
  )

  current_export <- list(
    creation_date = Sys.time(),
    output_format = match.arg(output_format),
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

.subset_motif_rows <- function(motifs, indices = NULL, n_motifs = NULL, seed = 222) {
  if (is.null(indices)) {
    if (!is.null(n_motifs)) {
      if (!is.numeric(n_motifs) || length(n_motifs) != 1 || is.na(n_motifs)) {
        stop("n_motifs must be a single numeric value")
      }

      n_motifs <- as.integer(n_motifs)
      if (n_motifs < 1) {
        stop("n_motifs must be >= 1")
      }

      if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
        stop("seed must be a single numeric value")
      }
      if (as.integer(seed) != 222L) {
        warning("seed is fixed to 222 in this repository; using 222")
      }
      set.seed(222L)

      if (n_motifs >= nrow(motifs)) {
        warning(sprintf(
          "n_motifs (%d) >= available motifs (%d); using all motifs",
          n_motifs,
          nrow(motifs)
        ))
        motifs$.source_index <- seq_len(nrow(motifs))
        return(motifs)
      }

      if ("day_post_hatch" %in% names(motifs)) {
        day_vals <- as.character(motifs$day_post_hatch)
        day_vals[is.na(day_vals) | day_vals == ""] <- "unknown_day"

        idx_by_day <- split(seq_len(nrow(motifs)), day_vals)
        sampled_idx <- unlist(lapply(idx_by_day, function(idx) {
          sample(idx, size = min(n_motifs, length(idx)), replace = FALSE)
        }), use.names = FALSE)
      } else {
        sampled_idx <- sample(seq_len(nrow(motifs)), size = n_motifs, replace = FALSE)
      }

      result <- motifs[sampled_idx, , drop = FALSE]
      result$.source_index <- sampled_idx
      return(result)
    }

    motifs$.source_index <- seq_len(nrow(motifs))
    return(motifs)
  }

  if (!is.numeric(indices) || anyNA(indices)) {
    stop("indices must be a numeric vector without NA values")
  }

  if (any(indices < 1) || any(indices > nrow(motifs))) {
    stop(sprintf("indices must be between 1 and %d", nrow(motifs)))
  }

  if (length(unique(indices)) != length(indices)) {
    warning("Duplicate indices detected; duplicates will be exported multiple times")
  }

  result <- motifs[indices, , drop = FALSE]
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

.prepare_motif_export_rows <- function(motifs) {
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

.export_motif_rows <- function(motifs,
                               wav_dir,
                               output_format,
                               output_dir,
                               hdf5_filename,
                               metadata_filename,
                               name_prefix,
                               overwrite,
                               write_metadata,
                               verbose) {
  counters <- new.env(parent = emptyenv())
  day_stats <- new.env(parent = emptyenv())
  metadata_rows <- vector("list", nrow(motifs))

  h5 <- NULL
  h5_path <- NULL
  sample_rates <- numeric(0)

  if (output_format == "hdf5") {
    ensure_pkgs("hdf5r")
    h5_path <- file.path(output_dir, hdf5_filename)
    if (file.exists(h5_path)) {
      if (!overwrite) {
        stop("HDF5 file already exists: ", h5_path, ". Set overwrite = TRUE to replace.")
      }
      file.remove(h5_path)
    }
    h5 <- hdf5r::H5File$new(filename = h5_path, mode = "w")
    on.exit(
      {
        if (!is.null(h5)) {
          h5$close_all()
        }
      },
      add = TRUE
    )
  }

  for (i in seq_len(nrow(motifs))) {
    row <- motifs[i, , drop = FALSE]
    day_summary <- as.character(row$day_post_hatch)
    if (is.na(day_summary) || day_summary == "") {
      day_summary <- "unknown_day"
    }
    if (!exists(day_summary, envir = day_stats, inherits = FALSE)) {
      assign(
        day_summary,
        c(total = 0L, written = 0L, missing = 0L, invalid = 0L, exists = 0L),
        envir = day_stats
      )
    }
    stats <- get(day_summary, envir = day_stats, inherits = FALSE)
    stats[["total"]] <- stats[["total"]] + 1L
    assign(day_summary, stats, envir = day_stats)

    has_day <- "day_post_hatch" %in% names(row) &&
      !is.na(row$day_post_hatch) &&
      row$day_post_hatch != "unknown_day"
    source_path <- if (has_day) {
      file.path(wav_dir, as.character(row$day_post_hatch), as.character(row$filename))
    } else {
      file.path(wav_dir, as.character(row$filename))
    }

    if (!file.exists(source_path)) {
      warning("Skipping missing source file: ", source_path)
      stats <- get(day_summary, envir = day_stats, inherits = FALSE)
      stats[["missing"]] <- stats[["missing"]] + 1L
      assign(day_summary, stats, envir = day_stats)
      next
    }
    duration <- row$end_time - row$start_time
    if (duration <= 0) {
      warning(sprintf(
        "Skipping motif index %s with non-positive duration",
        row$.source_index
      ))
      stats <- get(day_summary, envir = day_stats, inherits = FALSE)
      stats[["invalid"]] <- stats[["invalid"]] + 1L
      assign(day_summary, stats, envir = day_stats)
      next
    }

    bird_id <- .sanitize_group_name(as.character(row$bird_id))
    day <- .sanitize_group_name(as.character(row$day_post_hatch))

    counter_key <- paste(bird_id, day, sep = "::")
    if (!exists(counter_key, envir = counters, inherits = FALSE)) {
      assign(counter_key, 0L, envir = counters)
    }
    current_n <- get(counter_key, envir = counters, inherits = FALSE) + 1L
    assign(counter_key, current_n, envir = counters)

    clip_id <- sprintf("%s_%03d", name_prefix, current_n)

    if (output_format == "wav") {
      out_dir <- file.path(output_dir, "motifs", bird_id, day)
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
      }

      out_path <- file.path(out_dir, paste0(clip_id, ".wav"))
      if (file.exists(out_path) && !overwrite) {
        warning("Skipping existing file: ", out_path)
        stats <- get(day_summary, envir = day_stats, inherits = FALSE)
        stats[["exists"]] <- stats[["exists"]] + 1L
        assign(day_summary, stats, envir = day_stats)
        next
      }

      av::av_audio_convert(
        source_path,
        out_path,
        start_time = as.numeric(row$start_time),
        total_time = as.numeric(duration),
        verbose = FALSE
      )

      metadata_rows[[i]] <- data.frame(
        source_index = as.integer(row$.source_index),
        bird_id = bird_id,
        day_post_hatch = as.character(row$day_post_hatch),
        label = as.character(row$label),
        filename = as.character(row$filename),
        start_time = as.numeric(row$start_time),
        end_time = as.numeric(row$end_time),
        duration = as.numeric(duration),
        clip_id = clip_id,
        output_path = normalizePath(out_path, mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    } else {
      wave <- tuneR::readWave(
        filename = source_path,
        from = as.numeric(row$start_time),
        to = as.numeric(row$end_time),
        units = "seconds"
      )
      sample_rates <- c(sample_rates, wave@samp.rate)

      sample_vec <- if (isTRUE(wave@stereo)) {
        (as.numeric(wave@left) + as.numeric(wave@right)) / 2
      } else {
        as.numeric(wave@left)
      }

      norm_factor <- max(1, 2^(as.numeric(wave@bit) - 1))
      sample_vec <- as.single(sample_vec / norm_factor)

      bird_group <- .h5_get_or_create_group(h5, bird_id)
      day_group <- .h5_get_or_create_group(bird_group, day)
      ds <- day_group$create_dataset(name = clip_id, robj = sample_vec)
      hdf5r::h5attr(ds, "sample_rate") <- as.integer(wave@samp.rate)
      hdf5r::h5attr(ds, "source_file") <- as.character(row$filename)
      hdf5r::h5attr(ds, "start_time") <- as.numeric(row$start_time)
      hdf5r::h5attr(ds, "end_time") <- as.numeric(row$end_time)
      hdf5r::h5attr(ds, "label") <- as.character(row$label)
      ds$close()

      metadata_rows[[i]] <- data.frame(
        source_index = as.integer(row$.source_index),
        bird_id = bird_id,
        day_post_hatch = as.character(row$day_post_hatch),
        label = as.character(row$label),
        filename = as.character(row$filename),
        start_time = as.numeric(row$start_time),
        end_time = as.numeric(row$end_time),
        duration = as.numeric(duration),
        clip_id = clip_id,
        output_path = paste0("/", bird_id, "/", day, "/", clip_id),
        stringsAsFactors = FALSE
      )
    }

    stats <- get(day_summary, envir = day_stats, inherits = FALSE)
    stats[["written"]] <- stats[["written"]] + 1L
    assign(day_summary, stats, envir = day_stats)
  }

  export_meta <- do.call(rbind, metadata_rows[!sapply(metadata_rows, is.null)])
  if (is.null(export_meta)) {
    warning("No motif clips were exported")
    export_meta <- data.frame()
  }

  if (output_format == "hdf5" && !is.null(h5)) {
    hdf5r::h5attr(h5, "creation_date") <- as.character(Sys.time())
    hdf5r::h5attr(h5, "n_clips") <- as.integer(nrow(export_meta))
    if (length(sample_rates) > 0 && length(unique(sample_rates)) == 1) {
      hdf5r::h5attr(h5, "sample_rate") <- as.integer(sample_rates[[1]])
    } else if (length(sample_rates) > 0) {
      hdf5r::h5attr(h5, "sample_rate") <- "mixed"
    }
  }

  if (write_metadata && nrow(export_meta) > 0) {
    metadata_path <- file.path(output_dir, metadata_filename)
    utils::write.csv(export_meta, metadata_path, row.names = FALSE)
  }

  if (verbose) {
    day_keys <- sort(ls(envir = day_stats))
    for (day_key in day_keys) {
      stats <- get(day_key, envir = day_stats, inherits = FALSE)
      message(sprintf(
        "Day %s: %d written / %d total (skipped: %d missing, %d invalid, %d existing)",
        day_key,
        stats[["written"]],
        stats[["total"]],
        stats[["missing"]],
        stats[["invalid"]],
        stats[["exists"]]
      ))
    }
    message(sprintf("Exported %d motif clips to %s format", nrow(export_meta), output_format))
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
