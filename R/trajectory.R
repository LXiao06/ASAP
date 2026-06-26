# Song Trajectory Analysis
# Update date : Feb. 7, 2026

# Suppress R CMD check notes about internal functions used by trajectory analysis functions
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("across", "ends_with", "cor"))
}

#' Create Spectrogram Matrices for Song Trajectory Analysis
#'
#' @description
#' Creates spectrogram matrices from sliding windows for song trajectory analysis.
#'
#' @param x An object to process, either a data frame or SAP object
#' @param wav_dir Directory containing WAV files (for default method)
#' @param window_size Size of sliding window in seconds (default: 0.1)
#' @param step_size Step size between windows (default: 0.005)
#' @param wl Window length for spectrogram (default: 128)
#' @param ovlp Overlap percentage (default: 50)
#' @param fftw Logical, use FFTW or not (default: TRUE)
#' @param flim Frequency limits (default: c(1, 12))
#' @param cores Number of processing cores
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables', 'bouts', 'segments')
#' @param data_type For SAP objects: Type of data to analyze
#' @param clusters For SAP objects: Specific clusters to include
#' @param sample_percent For SAP objects: Percentage to sample
#' @param balanced For SAP objects: Whether to balance across groups
#' @param max_samples_per_label For SAP objects: Optional integer. Maximum
#'   number of samples to randomly analyze from each label when
#'   \code{balanced = FALSE}. If a label has fewer samples, all available
#'   samples are used. Ignored when \code{balanced = TRUE}.
#' @param labels For SAP objects: Specific labels to include
#' @param seed For SAP objects: Random seed
#' @param verbose For SAP objects: Whether to print progress
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' Creates trajectory matrix with the following steps:
#' \itemize{
#'   \item Generates sliding windows for each segment
#'   \item Computes spectrograms for each window
#'   \item Combines results into matrix form
#' }
#'
#' For SAP objects, additional features include:
#' \itemize{
#'   \item Support for different segment types
#'   \item Optional cluster/label filtering
#'   \item Balanced sampling options
#'   \item Results storage in features slot
#' }
#'
#' @return
#' For default method: A list containing:
#' \itemize{
#'   \item spectrogram_matrix: Matrix of spectrogram vectors
#'   \item sliding_windows: Data frame of window information
#' }
#'
#' For SAP objects: Updated SAP object with trajectory matrix stored in features slot
#'
#' @examples
#' \dontrun{
#' # Create trajectory matrix from segments
#' matrix <- create_trajectory_matrix(segments,
#'   wav_dir = "path/to/wavs",
#'   window_size = 0.1,
#'   step_size = 0.005
#' )
#'
#' # Create trajectory matrix from SAP object
#' sap_obj <- create_trajectory_matrix(sap_object,
#'   segment_type = "motifs",
#'   balanced = TRUE,
#'   sample_percent = 80
#' )
#'
#' # Cap unbalanced analysis to at most 200 samples per label
#' sap_obj <- create_trajectory_matrix(sap_object,
#'   segment_type = "motifs",
#'   max_samples_per_label = 200
#' )
#'
#' # Create matrix with specific clusters
#' sap_obj <- create_trajectory_matrix(sap_object,
#'   segment_type = "syllables",
#'   clusters = c(1, 2),
#'   labels = c("a", "b")
#' )
#' }
#'
#' @rdname create_trajectory_matrix
#' @export
create_trajectory_matrix <- function(x, ...) {
  UseMethod("create_trajectory_matrix")
}

#' @rdname create_trajectory_matrix
#' @export
create_trajectory_matrix.default <- function(
    x,
    wav_dir,
    window_size = 0.1,
    step_size = 0.005,
    wl = 128,
    ovlp = 50,
    fftw = TRUE,
    flim = c(1, 12),
    cores = NULL,
    ...) {
  # Input validation
  if (missing(x) || missing(wav_dir)) {
    stop("Both data frame containing audio segment information and wave file directory must be provided")
  }
  if (!is.data.frame(x)) {
    stop("Input x must be a data frame")
  }
  required_cols <- c("filename", "day_post_hatch", "label", "start_time", "end_time")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  if (!dir.exists(wav_dir)) {
    stop("Wave file directory does not exist: ", wav_dir)
  }

  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }
  sys_type <- Sys.info()["sysname"]

  # Generate constrained sliding windows
  sliding_windows <- do.call(rbind, lapply(
    seq_len(nrow(x)),
    function(i) create_sliding_window(i, x, window_size, step_size)
  ))

  # Split windows by original segment
  segments <- split(sliding_windows, sliding_windows$rendition)

  # Define core processing function
  process_segment <- function(segment) {
    tryCatch(
      {
        original_row <- x[segment$rendition[1], ]
        sound_path <- construct_wav_path(original_row, wav_dir = wav_dir)

        # Read audio with small buffer
        buffer <- max(window_size * 0.005, step_size)
        wave <- tuneR::readWave(
          sound_path,
          from = original_row$start_time,
          to = original_row$end_time + buffer,
          units = "seconds"
        )

        # Compute spectrogram
        sr <- wave@samp.rate
        if (fftw) ensure_pkgs("fftw")

        spg <- seewave::spectro(
          wave,
          f = sr,
          wl = wl,
          ovlp = ovlp,
          fftw = fftw,
          flim = flim,
          plot = FALSE,
          osc = FALSE,
          cont = FALSE
        )

        # Calculate time bin parameters
        time_vector <- spg$time
        time_per_bin <- time_vector[2] - time_vector[1]
        expected_bins <- floor(window_size / time_per_bin)

        # Process all windows in segment
        vectors <- lapply(seq_len(nrow(segment)), function(j) {
          sw <- segment[j, ]
          start_rel <- sw$.time
          end_rel <- start_rel + window_size

          # Find matching time bins
          idx <- which(time_vector >= start_rel & time_vector <= end_rel)

          # Extract exactly expected_bins bins
          idx <- idx[1:expected_bins]

          # Convert to vector and return
          as.vector(spg$amp[, idx, drop = FALSE])
        })

        # Combine vectors column-wise
        do.call(cbind, vectors)
      },
      error = function(e) {
        warning(sprintf(
          "Error processing segment %d: %s",
          segment$rendition[1], e$message
        ))
        return(NULL)
      }
    )
  }

  # Parallel processing setup
  message(sprintf("Generating spectrograms using %d cores...", cores))
  if (fftw) ensure_pkgs("fftw")
  spc_list <- parallel_apply(segments, process_segment, cores)
  # if (cores > 1) {
  #   if (sys_type == "Linux") {
  #     requireNamespace("pbmcapply", quietly = TRUE)
  #     spc_list <- pbmcapply::pbmclapply(
  #       segments,
  #       process_segment,
  #       mc.cores = cores
  #     )
  #   } else {
  #     # cl <- parallel::makeCluster(cores)
  #     # on.exit(parallel::stopCluster(cl))
  #     # parallel::clusterExport(
  #     #   cl,
  #     #   c("construct_wav_path", "create_sliding_window"),
  #     #   envir = environment()
  #     # )
  #     spc_list <- pbapply::pblapply(
  #       segments,
  #       process_segment,
  #       cl = cores
  #     )
  #   }
  # } else {
  #   spc_list <- pbapply::pblapply(segments, process_segment)
  # }

  # Combine results without building a second giant transposed copy.
  spc_list <- spc_list[!sapply(spc_list, is.null)]
  if (length(spc_list) == 0) {
    stop("No valid spectrograms were generated")
  }

  feature_counts <- vapply(spc_list, nrow, integer(1))
  window_counts <- vapply(spc_list, ncol, integer(1))
  if (length(unique(feature_counts)) != 1) {
    stop("Generated spectrogram blocks have inconsistent feature counts")
  }

  n_windows <- sum(window_counts)
  n_features <- feature_counts[1]
  estimated_gb <- (as.numeric(n_windows) * as.numeric(n_features) * 8) / 1024^3
  message(sprintf(
    "Combining %d windows x %d features into trajectory matrix (~%.1f GB)...",
    n_windows, n_features, estimated_gb
  ))

  gc(verbose = FALSE)
  spectrogram_matrix <- tryCatch(
    matrix(NA_real_, nrow = n_windows, ncol = n_features),
    error = function(e) {
      stop(sprintf(
        paste(
          "Cannot allocate trajectory matrix (~%.1f GB).",
          "Reduce sample_percent/max_samples_per_label, increase step_size,",
          "reduce window_size/ovlp/flim, or use an on-disk/chunked PCA workflow."
        ),
        estimated_gb
      ), call. = FALSE)
    }
  )
  row_start <- 1L
  while (length(spc_list) > 0) {
    block <- spc_list[[1]]
    row_end <- row_start + ncol(block) - 1L
    spectrogram_matrix[row_start:row_end, ] <- t(block)
    row_start <- row_end + 1L
    spc_list[[1]] <- NULL
  }
  gc(verbose = FALSE)

  cat("Spectrogram generation complete!\n")
  list(
    spectrogram_matrix = spectrogram_matrix,
    sliding_windows = sliding_windows
  )
}


#' @rdname create_trajectory_matrix
#' @export
create_trajectory_matrix.Sap <- function(
    x,
    segment_type = c("motifs", "syllables", "bouts", "segments"),
    data_type = NULL,
    clusters = NULL,
    sample_percent = NULL,
    balanced = FALSE,
    max_samples_per_label = NULL,
    labels = NULL,
    seed = 222,
    window_size = 0.1,
    step_size = 0.005,
    wl = 128,
    ovlp = 50,
    fftw = TRUE,
    flim = c(1, 12),
    cores = NULL,
    verbose = TRUE,
    ...) {
  if (verbose) {
    message(sprintf(
      "\n=== Starting trajectory analysis for %s using %s%s ===",
      segment_type[1],
      data_type[1],
      ifelse(is.null(clusters), "", sprintf(" (cluster %s)", clusters))
    ))
  }

  # Validate segment_type
  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type) # Remove 's' from end

  # Get appropriate data frame based on data_type
  if (is.null(data_type)) {
    segments_df <- x[[segment_type]]
    if (!inherits(segments_df, "segment") || nrow(segments_df) == 0) {
      stop("No segments found in the specified segment type")
    }
  } else if (data_type == "feat.embeds") {
    # Get feature data
    feature_data <- x$features[[feature_type]][[data_type]]

    if (is.null(feature_data)) {
      stop(sprintf("No %s data found in %s features", data_type, feature_type))
    }

    # Select only key columns and embedding columns from feat.embeds
    # This prevents duplicate metadata columns while restoring 'subfolder'
    embed_cols <- c(
      "filename", "start_time", "end_time", "label", "day_post_hatch",
      "UMAP1", "UMAP2", "cluster"
    )
    embed_cols <- intersect(embed_cols, names(feature_data))
    key_cols <- intersect(embed_cols, names(x[[segment_type]]))

    segments_df <- x[[segment_type]] |>
      dplyr::inner_join(
        feature_data[, embed_cols],
        by = key_cols,
        relationship = "many-to-many"
      )
  } else {
    stop("Invalid data_type. Must be NULL or 'feat.embeds'")
  }

  # Validate clusters
  if (!is.null(clusters)) {
    if (!"cluster" %in% colnames(segments_df)) {
      stop("No cluster column found in the data")
    }

    # Check which clusters are missing
    available_clusters <- unique(segments_df$cluster)
    missing_clusters <- setdiff(clusters, available_clusters)
    if (length(missing_clusters) > 0) {
      stop(sprintf(
        "The following clusters were not found: %s",
        paste(missing_clusters, collapse = ", ")
      ))
    }
  }

  # Preserve original row index for downstream traceability
  # (done before select_segments so .source_row persists through filtering)
  segments_df$.source_row <- seq_len(nrow(segments_df))

  # Apply segment selection
  segments_df <- select_segments(segments_df,
    clusters = clusters,
    balanced = balanced,
    sample_percent = sample_percent,
    max_samples_per_label = max_samples_per_label,
    labels = labels,
    seed = seed
  )

  # Call default method to generate spectrograms
  if (fftw) ensure_pkgs("fftw")

  result <- create_trajectory_matrix.default(
    x = segments_df,
    wav_dir = x$base_path,
    window_size = window_size,
    step_size = step_size,
    wl = wl,
    ovlp = ovlp,
    fftw = fftw,
    flim = flim,
    cores = cores,
    ...
  )

  # Save results in the Sap object using feature_type
  x$features[[feature_type]][["traj_mat"]] <- result$spectrogram_matrix
  x$features[[feature_type]][["traj.embeds"]] <- result$sliding_windows

  # Add parameters as attributes
  attrs <- list(
    window_size = window_size,
    step_size = step_size,
    wl = wl,
    ovlp = ovlp,
    flim = flim,
    data_type = data_type,
    clusters = clusters,
    balanced = balanced,
    sample_percent = sample_percent,
    max_samples_per_label = max_samples_per_label,
    labels = labels,
    seed = seed
  )
  attributes(x$features[[feature_type]][["traj.embeds"]]) <- c(
    attributes(x$features[[feature_type]][["traj.embeds"]]),
    attrs
  )

  invisible(x)
}

#' Create Sliding Windows for Segment
#'
#' @description
#' Internal function to generate overlapping time windows for a segment.
#'
#' @param i Rendition number
#' @param x Data frame with segment information
#' @param window_size Window size in seconds
#' @param step_size Step size in seconds
#' @param ... Additional arguments
#'
#' @return
#' Data frame containing sliding window information
#'
#' @keywords internal
create_sliding_window <- function(
    i,
    x,
    window_size = 0.1,
    step_size = 0.005,
    ...) {
  # Constrain window generation to ensure sufficient bins
  start <- x$start_time[i]
  end <- x$end_time[i]

  # Calculate max windows that fit in the segment
  max_windows <- floor((end - start - window_size) / step_size) + 1
  relative_time <- round((0:(max_windows - 1)) * step_size, digits = 6)
  window_start_time <- start + relative_time

  # Create sliding window data frame
  sliding_windows <- data.frame(
    filename = x$filename[i],
    day_post_hatch = x$day_post_hatch[i],
    label = x$label[i],
    rendition = i,
    selec = seq_along(window_start_time),
    start_time = window_start_time,
    end_time = window_start_time + window_size,
    .time = relative_time,
    .source_row = if (".source_row" %in% names(x)) x$.source_row[i] else NA_integer_,
    stringsAsFactors = FALSE
  )

  return(sliding_windows)
}


# Filter Trajectory Outliers
# Update date : May 1, 2026

#' Filter Outlier Renditions from Trajectory Embeddings
#'
#' @description
#' Identifies and removes entire renditions (trials) from trajectory embeddings
#' when a substantial fraction of their time steps are statistical outliers
#' relative to the group distribution. Outlier detection is performed per label,
#' per time step, and per embedding dimension using the IQR method.
#'
#' @param x An object to filter: a trajectory embeddings data frame or SAP object
#' @param dims Character vector of dimension columns to assess (default: c("PC1", "PC2"))
#' @param iqr_multiplier Multiplier applied to IQR to define outlier bounds (default: 4)
#' @param min_outlier_fraction Minimum fraction of time steps that must be outliers
#'   before the entire rendition is removed (default: 0.1, i.e., 10%)
#' @param labels Optional character vector to restrict filtering to specific labels
#' @param plot Logical. If \code{TRUE}, plot a heatmap of the filtered trajectories (default: \code{FALSE}).
#' @param balanced Logical. Whether to balance segment numbers across groups for the heatmap (default: \code{FALSE}).
#' @param ordered Logical. Whether to order segments for the heatmap (default: \code{FALSE}).
#' @param clusters Integer vector of clusters to include for the heatmap (default: \code{NULL}).
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables', 'bouts', 'segments')
#' @param verbose Whether to print a removal summary (default: TRUE)
#' @param ... Additional arguments
#'
#' @details
#' The outlier detection algorithm works as follows:
#' \itemize{
#'   \item For each label separately, at each time step, compute Q1, Q3, and IQR
#'     across all renditions for each specified dimension.
#'   \item A time step is flagged as an outlier for a given rendition if its value
#'     in \strong{any} dimension falls outside
#'     \code{[Q1 - iqr_multiplier * IQR, Q3 + iqr_multiplier * IQR]}.
#'   \item A rendition is removed entirely if the fraction of flagged time steps
#'     exceeds \code{min_outlier_fraction}.
#' }
#'
#' Outlier bounds are computed within each label to avoid penalizing conditions
#' with naturally different mean trajectories.
#'
#' @return
#' For the default method: the filtered trajectory embeddings data frame (same
#' structure as input, with outlier renditions removed).
#'
#' For SAP objects: the updated SAP object with filtered \code{traj.embeds}
#' stored back in \code{x$features[[feature_type]]$traj.embeds}.
#'
#' @examples
#' \dontrun{
#' # Filter from data frame directly
#' traj_clean <- filter_trajectory_outliers(
#'   sap$features$motif$traj.embeds,
#'   dims = c("PC1", "PC2"),
#'   iqr_multiplier = 4,
#'   min_outlier_fraction = 0.1
#' )
#'
#' # Filter from SAP object (updates traj.embeds in place)
#' sap <- filter_trajectory_outliers(sap)
#'
#' # Stricter: remove renditions with > 5% outlier time points
#' sap <- filter_trajectory_outliers(sap, min_outlier_fraction = 0.05)
#'
#' # Use filtered data directly with trajectory_dispersion
#' result <- trajectory_dispersion(traj_clean, dims = c("PC1", "PC2"))
#' }
#'
#' @rdname filter_trajectory_outliers
#' @export
filter_trajectory_outliers <- function(x, ...) {
  UseMethod("filter_trajectory_outliers")
}


#' @rdname filter_trajectory_outliers
#' @export
filter_trajectory_outliers.default <- function(x,
                                               dims = c("PC1", "PC2"),
                                               iqr_multiplier = 4,
                                               min_outlier_fraction = 0.1,
                                               labels = NULL,
                                               verbose = TRUE,
                                               ...) {
  if (!is.data.frame(x)) stop("Input must be a data frame")

  required_cols <- c("label", "rendition", ".time")
  missing_cols <- setdiff(c(required_cols, dims), names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Optionally restrict to specific labels
  if (!is.null(labels)) {
    missing_labels <- setdiff(labels, unique(x$label))
    if (length(missing_labels) > 0) {
      stop(sprintf("Labels not found: %s", paste(missing_labels, collapse = ", ")))
    }
    x <- x[x$label %in% labels, ]
  }

  all_labels <- unique(x$label)
  total_before <- length(unique(x$rendition))
  summary_rows <- list()

  kept_renditions <- unlist(lapply(all_labels, function(lbl) {
    lbl_data <- x[x$label == lbl, ]
    renditions <- unique(lbl_data$rendition)
    time_steps <- sort(unique(lbl_data$.time))

    # For each time step, compute IQR bounds per dimension
    bounds <- do.call(rbind, lapply(time_steps, function(t) {
      t_rows <- lbl_data[lbl_data$.time == t, , drop = FALSE]
      if (nrow(t_rows) < 4) {
        return(NULL)
      } # need enough points for IQR

      row_data <- data.frame(.time = t)
      for (d in dims) {
        vals <- t_rows[[d]]
        q1 <- stats::quantile(vals, 0.25, na.rm = TRUE)
        q3 <- stats::quantile(vals, 0.75, na.rm = TRUE)
        iqr <- q3 - q1
        row_data[[paste0(d, "_lo")]] <- q1 - iqr_multiplier * iqr
        row_data[[paste0(d, "_hi")]] <- q3 + iqr_multiplier * iqr
      }
      row_data
    }))

    if (is.null(bounds) || nrow(bounds) == 0) {
      return(renditions)
    }

    # Merge bounds back onto individual rows
    lbl_merged <- merge(lbl_data, bounds, by = ".time", all.x = TRUE)

    # Flag each row: outlier if ANY dim is out of bounds
    is_outlier <- Reduce(`|`, lapply(dims, function(d) {
      lo_col <- paste0(d, "_lo")
      hi_col <- paste0(d, "_hi")
      val <- lbl_merged[[d]]
      lo <- lbl_merged[[lo_col]]
      hi <- lbl_merged[[hi_col]]
      (!is.na(val) & !is.na(lo) & (val < lo | val > hi))
    }))
    lbl_merged$is_outlier <- is_outlier

    # Compute outlier fraction per rendition
    rend_stats <- do.call(rbind, lapply(renditions, function(r) {
      r_rows <- lbl_merged[lbl_merged$rendition == r, ]
      n_total <- nrow(r_rows)
      n_outlier <- sum(r_rows$is_outlier, na.rm = TRUE)
      data.frame(
        rendition = r,
        n_total = n_total,
        n_outlier = n_outlier,
        frac = n_outlier / n_total
      )
    }))

    keep <- rend_stats$rendition[rend_stats$frac <= min_outlier_fraction]
    remove <- rend_stats$rendition[rend_stats$frac > min_outlier_fraction]

    summary_rows[[lbl]] <<- data.frame(
      label = lbl,
      n_before = length(renditions),
      n_removed = length(remove),
      n_after = length(keep),
      pct_removed = round(100 * length(remove) / length(renditions), 1),
      stringsAsFactors = FALSE
    )

    keep
  }))

  # Filter data frame
  x_filtered <- x[x$rendition %in% kept_renditions, ]

  # Print summary
  if (verbose) {
    summary_df <- do.call(rbind, summary_rows)
    message("\n=== Trajectory Outlier Filtering ===")
    message(sprintf(
      "IQR multiplier : %.1f  |  Min outlier fraction: %.0f%%",
      iqr_multiplier, min_outlier_fraction * 100
    ))
    message(sprintf("Dimensions     : %s\n", paste(dims, collapse = ", ")))
    print(summary_df, row.names = FALSE)
    total_after <- length(unique(x_filtered$rendition))
    message(sprintf(
      "\nTotal renditions: %d -> %d  (removed %d, %.1f%%)",
      total_before, total_after,
      total_before - total_after,
      100 * (total_before - total_after) / total_before
    ))
  }

  return(x_filtered)
}


#' @rdname filter_trajectory_outliers
#' @export
filter_trajectory_outliers.Sap <- function(x,
                                           segment_type = c(
                                             "motifs", "syllables",
                                             "bouts", "segments"
                                           ),
                                           dims = c("PC1", "PC2"),
                                           iqr_multiplier = 4,
                                           min_outlier_fraction = 0.1,
                                           labels = NULL,
                                           plot = FALSE,
                                           # plot_heatmap arguments (only used when plot = TRUE)
                                           balanced = FALSE,
                                           ordered = FALSE,
                                           clusters = NULL,
                                           verbose = TRUE,
                                           ...) {
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  traj_embeds <- x$features[[feature_type]][["traj.embeds"]]
  if (is.null(traj_embeds)) {
    stop(sprintf(
      "Trajectory embeddings not found for '%s'. Run create_trajectory_matrix() first.",
      segment_type
    ))
  }

  missing_dims <- setdiff(dims, names(traj_embeds))
  if (length(missing_dims) > 0) {
    stop(sprintf(
      "Dimensions not found in traj.embeds: %s",
      paste(missing_dims, collapse = ", ")
    ))
  }

  filtered <- filter_trajectory_outliers.default(
    x                    = traj_embeds,
    dims                 = dims,
    iqr_multiplier       = iqr_multiplier,
    min_outlier_fraction = min_outlier_fraction,
    labels               = labels,
    verbose              = verbose
  )

  x$features[[feature_type]][["traj.embeds"]] <- filtered

  # Optionally plot heatmap of retained motifs
  if (plot) {
    if (!"selec" %in% names(filtered)) {
      warning("Cannot plot heatmap: 'selec' column not found in traj.embeds")
    } else {
      # Each rendition's first sliding window (selec == 1) carries the
      # original motif start_time, allowing us to match back to sap$motifs
      surviving_motifs <- filtered |>
        dplyr::filter(.data$selec == 1) |>
        dplyr::select(.data$filename, .data$label,
          motif_start = .data$start_time
        ) |>
        dplyr::distinct()

      # Build a temporary SAP with only the retained motifs
      x_temp <- x
      x_temp[[segment_type]] <- x[[segment_type]] |>
        dplyr::inner_join(surviving_motifs,
          by = c("filename", "label",
            "start_time" = "motif_start"
          )
        )

      if (nrow(x_temp[[segment_type]]) == 0) {
        warning("No motifs matched after filtering; skipping heatmap")
      } else {
        if (verbose) {
          message(sprintf(
            "\nPlotting heatmap for %d retained %s...",
            nrow(x_temp[[segment_type]]), segment_type
          ))
        }
        # Pass balanced/ordered/clusters explicitly; all other plot_heatmap
        # args (color_palette, msmooth, contrast, etc.) flow through ...
        # so plot_heatmap handles its own type conversions and defaults
        plot_heatmap(x_temp,
          segment_type = segment_type,
          balanced     = balanced,
          ordered      = ordered,
          clusters     = clusters,
          labels       = labels,
          ...
        )
      }
    }
  }

  invisible(x)
}


# # Trajectory Dispersion
# Update date : Jun. 23, 2026

#' Trajectory Dispersion Analysis
#'
#' @description
#' Quantifies trajectory dispersion across experimental conditions using three
#' complementary metrics computed in PCA or UMAP embedding space.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or SAP object
#' @param dims Character vector of dimension columns to use (default: c("PC1", "PC2"))
#' @param labels Optional character vector of labels to include
#' @param max_pairs Maximum number of pairwise comparisons per label (default: 5000)
#' @param seed Random seed for reproducible pair sampling (default: 222)
#' @param stats Logical. If \code{TRUE} (default), run Kruskal-Wallis and pairwise
#'   Wilcoxon tests. Set to \code{FALSE} to skip statistical testing and return
#'   \code{NULL} for \code{tests}.
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables', 'bouts', 'segments')
#' @param ... Additional arguments
#'
#' @details
#' Three metrics are computed:
#' \describe{
#'   \item{Mean Pairwise Distance}{Average Euclidean distance between all pairs of
#'     trajectories within each condition, measured at matched time points and averaged
#'     across time. Higher values indicate greater spread among renditions.}
#'   \item{Centroid Dispersion}{For each rendition, the mean distance to the condition's
#'     centroid trajectory at each time point. Provides a per-rendition measure of how
#'     far each trajectory deviates from the "average" path.}
#'   \item{Path Length}{Total Euclidean distance traveled along each rendition's trajectory
#'     through the embedding space. Captures the overall complexity and extent of each
#'     trajectory, independent of the group mean.}
#' }
#'
#' Statistical testing is performed automatically:
#' \itemize{
#'   \item Kruskal-Wallis test for overall group differences
#'   \item Pairwise Wilcoxon rank-sum tests with Bonferroni correction
#' }
#'
#' @return
#' For default method: A list (returned invisibly) with the following elements:
#' \itemize{
#'   \item \code{pairwise}: Data frame of pairwise distances (label, pair_id, mean_dist)
#'   \item \code{dispersion}: Data frame of centroid dispersion (label, rendition, dispersion)
#'   \item \code{path_length}: Data frame of path lengths (label, rendition, path_length)
#'   \item \code{summary}: Summary table with mean and SD for each metric per label
#'   \item \code{tests}: List of statistical test results (\code{NULL} when \code{stats = FALSE})
#'   \item \code{type}: Character string \code{"dispersion"}, used by
#'     \code{plot_trajectory_variability()} for dispatch
#' }
#' For SAP objects: The updated SAP object with results stored in
#'   \code{x$features[[feature_type]][["trajectory_dispersion"]]} (returned invisibly).
#'
#' @examples
#' \dontrun{
#' # From SAP object using PC dimensions
#' result <- trajectory_dispersion(sap)
#'
#' # Using UMAP dimensions
#' result <- trajectory_dispersion(sap, dims = c("UMAP1", "UMAP2"))
#'
#' # From trajectory embeddings data frame directly
#' result <- trajectory_dispersion(sap$features$motif$traj.embeds,
#'   dims = c("PC1", "PC2")
#' )
#'
#' # Access results
#' result$summary # summary table
#' result$tests # statistical tests
#' result$dispersion # per-rendition dispersion values
#' }
#'
#' @rdname trajectory_dispersion
#' @export
