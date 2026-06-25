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
    .time = relative_time
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
#' @param iqr_multiplier Multiplier applied to IQR to define outlier bounds (default: 3)
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
#'   iqr_multiplier = 3,
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
                                               iqr_multiplier = 3,
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
                                           iqr_multiplier = 3,
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
trajectory_dispersion <- function(x, ...) {
  UseMethod("trajectory_dispersion")
}


#' @rdname trajectory_dispersion
#' @export
trajectory_dispersion.default <- function(x,
                                          dims = c("PC1", "PC2"),
                                          labels = NULL,
                                          max_pairs = 5000,
                                          seed = 222,
                                          stats = TRUE,
                                          verbose = TRUE,
                                          ...) {
  # Input validation
  if (!is.data.frame(x)) stop("Input must be a data frame")

  required_cols <- c("label", "rendition", ".time")
  missing_cols <- setdiff(c(required_cols, dims), names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  ensure_pkgs("ggplot2", "dplyr")

  # Filter labels if specified
  if (!is.null(labels)) {
    available <- unique(x$label)
    missing_labels <- setdiff(labels, available)
    if (length(missing_labels) > 0) {
      stop(sprintf(
        "Labels not found: %s\nAvailable: %s",
        paste(missing_labels, collapse = ", "),
        paste(available, collapse = ", ")
      ))
    }
    x <- x[x$label %in% labels, ]
  }

  all_labels <- unique(x$label)
  if (length(all_labels) < 2) {
    stop("At least two labels are required for variability comparison")
  }

  # Auto-disable stats when there are many groups
  if (stats && length(all_labels) > 6) {
    message(sprintf(
      "Note: %d labels detected. Setting stats = FALSE (statistical tests are not meaningful for this many groups).",
      length(all_labels)
    ))
    stats <- FALSE
  }

  if (verbose) {
    message("\n=== Trajectory Dispersion Analysis ===")
    message(sprintf("Dimensions: %s", paste(dims, collapse = ", ")))
    message(sprintf("Labels: %s", paste(all_labels, collapse = ", ")))
  }

  set.seed(seed)

  # ==== Metric 1: Mean Pairwise Distance ====
  if (verbose) message("\nComputing mean pairwise trajectory distances...")

  pairwise_results <- do.call(rbind, lapply(all_labels, function(lbl) {
    lbl_data <- x[x$label == lbl, ]
    renditions <- sort(unique(lbl_data$rendition))
    n_rend <- length(renditions)

    if (n_rend < 2) {
      warning(sprintf("Label '%s' has fewer than 2 renditions, skipping pairwise", lbl))
      return(NULL)
    }

    # Generate all pairs
    pairs <- utils::combn(renditions, 2)
    n_pairs <- ncol(pairs)

    # Sample if too many pairs
    if (n_pairs > max_pairs) {
      idx <- sample(n_pairs, max_pairs)
      pairs <- pairs[, idx, drop = FALSE]
    }

    # Reshape label data into a list of matrices keyed by rendition
    rend_list <- split(lbl_data, lbl_data$rendition)
    rend_matrices <- lapply(rend_list, function(rd) {
      rd <- rd[order(rd$.time), ]
      list(time = rd$.time, coords = as.matrix(rd[, dims, drop = FALSE]))
    })

    # Compute mean distance for each step
    pair_dists <- vapply(seq_len(ncol(pairs)), function(p) {
      r1 <- rend_matrices[[as.character(pairs[1, p])]]
      r2 <- rend_matrices[[as.character(pairs[2, p])]]

      # Match by time
      common_times <- intersect(r1$time, r2$time)
      if (length(common_times) == 0) {
        return(NA_real_)
      }

      idx1 <- match(common_times, r1$time)
      idx2 <- match(common_times, r2$time)

      # Euclidean distance at each time point, then average
      diffs <- r1$coords[idx1, , drop = FALSE] - r2$coords[idx2, , drop = FALSE]
      mean(sqrt(rowSums(diffs^2)), na.rm = TRUE)
    }, numeric(1))

    data.frame(
      label = lbl,
      pair_id = seq_along(pair_dists),
      mean_dist = pair_dists,
      stringsAsFactors = FALSE
    )
  }))

  if (verbose) message(sprintf("  %d total pairwise comparisons", nrow(pairwise_results)))

  # ==== Metric 2: Centroid Dispersion ====
  if (verbose) message("Computing centroid dispersion...")

  # Compute centroids per label per time point
  centroids <- x |>
    dplyr::group_by(.data$label, .data$.time) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(dims), mean),
      .groups = "drop"
    )

  centroid_dim_names <- paste0("centroid_", dims)
  names(centroids)[match(dims, names(centroids))] <- centroid_dim_names

  # Merge centroids with original data
  merged <- dplyr::left_join(x, centroids, by = c("label", ".time"))

  # Compute distance to centroid for each row
  merged$dist_to_centroid <- sqrt(rowSums(
    (as.matrix(merged[, dims]) - as.matrix(merged[, centroid_dim_names]))^2
  ))

  # Average across time for each rendition
  dispersion_results <- merged |>
    dplyr::group_by(.data$label, .data$rendition) |>
    dplyr::summarise(
      dispersion = mean(.data$dist_to_centroid, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()

  # ==== Metric 3: Path Length ====
  if (verbose) message("Computing trajectory path lengths...")

  path_length_results <- do.call(rbind, lapply(all_labels, function(lbl) {
    lbl_data <- x[x$label == lbl, ]
    renditions <- unique(lbl_data$rendition)

    do.call(rbind, lapply(renditions, function(r) {
      rend_data <- lbl_data[lbl_data$rendition == r, ]
      rend_data <- rend_data[order(rend_data$.time), ]
      coords <- as.matrix(rend_data[, dims, drop = FALSE])
      diffs <- diff(coords)
      path_len <- sum(sqrt(rowSums(diffs^2)))
      data.frame(
        label = lbl, rendition = r, path_length = path_len,
        stringsAsFactors = FALSE
      )
    }))
  }))

  # ==== Summary Table ====
  summary_pw <- pairwise_results |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(
      pairwise_mean = mean(.data$mean_dist, na.rm = TRUE),
      pairwise_sd = sd(.data$mean_dist, na.rm = TRUE),
      .groups = "drop"
    )

  summary_disp <- dispersion_results |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(
      dispersion_mean = mean(.data$dispersion, na.rm = TRUE),
      dispersion_sd = sd(.data$dispersion, na.rm = TRUE),
      .groups = "drop"
    )

  summary_pl <- path_length_results |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(
      path_length_mean = mean(.data$path_length, na.rm = TRUE),
      path_length_sd = sd(.data$path_length, na.rm = TRUE),
      .groups = "drop"
    )

  summary_table <- summary_pw |>
    dplyr::left_join(summary_disp, by = "label") |>
    dplyr::left_join(summary_pl, by = "label")

  # ==== Statistical Tests ====
  if (stats) {
    # Kruskal-Wallis omnibus tests
    test_pw <- stats::kruskal.test(mean_dist ~ label, data = pairwise_results)
    test_disp <- stats::kruskal.test(dispersion ~ label, data = dispersion_results)
    test_pl <- stats::kruskal.test(path_length ~ label, data = path_length_results)

    # Pairwise Wilcoxon post-hoc tests
    posthoc_pw <- stats::pairwise.wilcox.test(
      pairwise_results$mean_dist, pairwise_results$label,
      p.adjust.method = "bonferroni"
    )
    posthoc_disp <- stats::pairwise.wilcox.test(
      dispersion_results$dispersion, dispersion_results$label,
      p.adjust.method = "bonferroni"
    )
    posthoc_pl <- stats::pairwise.wilcox.test(
      path_length_results$path_length, path_length_results$label,
      p.adjust.method = "bonferroni"
    )

    tests <- list(
      kruskal = list(pairwise = test_pw, dispersion = test_disp, path_length = test_pl),
      posthoc = list(pairwise = posthoc_pw, dispersion = posthoc_disp, path_length = posthoc_pl)
    )

    if (verbose) {
      message("\n--- Summary ---")
      print(summary_table)
      message("\nKruskal-Wallis tests:")
      message(sprintf(
        "  Pairwise distance:   chi-sq = %.2f, p = %.2e",
        test_pw$statistic, test_pw$p.value
      ))
      message(sprintf(
        "  Centroid dispersion: chi-sq = %.2f, p = %.2e",
        test_disp$statistic, test_disp$p.value
      ))
      message(sprintf(
        "  Path length:         chi-sq = %.2f, p = %.2e",
        test_pl$statistic, test_pl$p.value
      ))
      message("\nPairwise Wilcoxon post-hoc (Bonferroni adjusted):")
      message("  Pairwise distance:")
      print(posthoc_pw$p.value)
      message("  Centroid dispersion:")
      print(posthoc_disp$p.value)
      message("  Path length:")
      print(posthoc_pl$p.value)
    }
  } else {
    tests <- NULL
    if (verbose) {
      message("\n--- Summary ---")
      print(summary_table)
    }
  }

  # Return results
  invisible(list(
    type = "dispersion",
    dims = dims,
    pairwise = pairwise_results,
    dispersion = dispersion_results,
    path_length = path_length_results,
    summary = summary_table,
    tests = tests
  ))
}


#' @rdname trajectory_dispersion
#' @export
trajectory_dispersion.Sap <- function(x,
                                      segment_type = c(
                                        "motifs", "syllables",
                                        "bouts", "segments"
                                      ),
                                      dims = c("PC1", "PC2"),
                                      labels = NULL,
                                      max_pairs = 5000,
                                      seed = 222,
                                      stats = TRUE,
                                      verbose = TRUE,
                                      ...) {
  # Validate
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  traj_embeds <- x$features[[feature_type]][["traj.embeds"]]
  if (is.null(traj_embeds)) {
    stop(sprintf(
      "Trajectory embeddings not found for %s. Run create_trajectory_matrix() first.",
      segment_type
    ))
  }

  # Check that requested dims exist
  missing_dims <- setdiff(dims, names(traj_embeds))
  if (length(missing_dims) > 0) {
    stop(sprintf(
      "Dimensions not found in traj.embeds: %s\nAvailable: %s",
      paste(missing_dims, collapse = ", "),
      paste(
        setdiff(
          names(traj_embeds),
          c(
            "filename", "day_post_hatch", "label", "rendition",
            "selec", "start_time", "end_time", ".time"
          )
        ),
        collapse = ", "
      )
    ))
  }

  result <- trajectory_dispersion.default(
    x = traj_embeds,
    dims = dims,
    labels = labels,
    max_pairs = max_pairs,
    seed = seed,
    stats = stats,
    verbose = verbose,
    ...
  )

  # Write back to SAP object
  x$features[[feature_type]][["trajectory_dispersion"]] <- result

  invisible(x)
}



# # Trajectory Path Deviation
# Update date : Jun. 23, 2026

#' Trajectory Path Deviation Analysis
#'
#' @description
#' Quantifies rendition-to-rendition trajectory deviation by measuring residual
#' spread around each label's mean trajectory and decomposing that spread into
#' orthogonal and parallel components relative to the local trajectory tangent.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or SAP object
#' @param dims Character vector of dimension columns to use (default: c("PC1", "PC2"))
#' @param trim_fraction Trim fraction for robust mean trajectory estimation
#'   (default: 0.1)
#' @param min_coverage Minimum fraction of renditions that must cover a time step
#'   for it to contribute to the reference trajectory (default: 0.5)
#' @param time_digits Number of decimal places used to bin \code{.time} before
#'   grouping and matching trajectories (default: \code{6}).
#' @param labels Optional character vector of labels to include
#' @param stats Logical. If \code{TRUE} (default), run Kruskal-Wallis and pairwise
#'   Wilcoxon tests. Set to \code{FALSE} to skip statistical testing and return
#'   \code{NULL} for \code{tests}.
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables',
#'   'bouts', 'segments')
#' @param ... Additional arguments
#'
#' @details
#' For each label, the function builds a robust mean trajectory in the requested
#' dimensions after binning \code{.time} to \code{time_digits}, estimates a
#' local tangent vector at each retained time step, and decomposes each
#' rendition's residual into:
#' \describe{
#'   \item{Total RMS Residual}{Overall deviation from the label-specific mean trajectory}
#'   \item{Orthogonal RMS Residual}{Deviation perpendicular to the local tangent;
#'     a direct measure of trajectory width / jitter around the backbone}
#'   \item{Parallel RMS Residual}{Deviation along the local tangent; often reflects
#'     phase or advance-lag variability}
#' }
#'
#' @return
#' For default method: A list (returned invisibly) with the following elements:
#' \itemize{
#'   \item \code{type}: Character string \code{"path_deviation"}, used by
#'     \code{plot_trajectory_variability()} for dispatch
#'   \item \code{width}: Per-rendition width metrics
#'   \item \code{summary}: Summary table with mean and SD for each metric per label
#'   \item \code{mean_trajectories}: Label-specific mean trajectories
#'   \item \code{tangent_vectors}: Label-specific unit tangent vectors
#'   \item \code{tests}: Kruskal-Wallis and pairwise Wilcoxon tests (\code{NULL}
#'     when \code{stats = FALSE} or only one label is present)
#' }
#' For SAP objects: The updated SAP object with results stored in
#'   \code{x$features[[feature_type]][["trajectory_path_deviation"]]} (returned invisibly).
#'
#' Use \code{\link{plot_trajectory_variability}(result)} to visualise the output.
#'
#' @examples
#' \dontrun{
#' result <- trajectory_path_deviation(sap, dims = c("PC1", "PC2"))
#' result$summary
#' result$width
#' }
#'
#' @export
trajectory_path_deviation <- function(x, ...) {
  UseMethod("trajectory_path_deviation")
}


#' @rdname trajectory_path_deviation
#' @export
trajectory_path_deviation.default <- function(x,
                                              dims = c("PC1", "PC2"),
                                              trim_fraction = 0.1,
                                              min_coverage = 0.5,
                                              time_digits = 6,
                                              labels = NULL,
                                              stats = TRUE,
                                              verbose = TRUE,
                                              ...) {
  if (!is.data.frame(x)) stop("Input must be a data frame")
  if (length(dims) < 2) {
    stop("Use at least two dimensions so orthogonal width is well-defined")
  }

  required_cols <- c("label", "rendition", ".time")
  missing_cols <- setdiff(c(required_cols, dims), names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  ensure_pkgs("ggplot2", "dplyr")

  if (!is.null(labels)) {
    missing_labels <- setdiff(labels, unique(x$label))
    if (length(missing_labels) > 0) {
      stop(sprintf("Labels not found: %s", paste(missing_labels, collapse = ", ")))
    }
    x <- x[x$label %in% labels, ]
  }

  x <- bin_trajectory_time_data(x, dims, time_digits)

  all_labels <- unique(x$label)
  if (length(all_labels) == 0) stop("No labels available after filtering")

  # Auto-disable stats when there are many groups
  if (stats && length(all_labels) > 6) {
    message(sprintf(
      "Note: %d labels detected. Setting stats = FALSE (statistical tests are not meaningful for this many groups).",
      length(all_labels)
    ))
    stats <- FALSE
  }

  if (verbose) {
    message("\n=== Trajectory Path Deviation Analysis ===")
    message(sprintf("Dimensions    : %s", paste(dims, collapse = ", ")))
    message(sprintf("Trim fraction : %.0f%% each tail", trim_fraction * 100))
    message(sprintf(
      "Min coverage  : %.0f%% of renditions per time step",
      min_coverage * 100
    ))
    message(sprintf("Time binning  : %d decimal places", time_digits))
    message(sprintf("Labels        : %s\n", paste(all_labels, collapse = ", ")))
  }

  compute_tangents <- function(coords) {
    n <- nrow(coords)
    p <- ncol(coords)
    tangents <- matrix(NA_real_, nrow = n, ncol = p)

    if (n < 2) {
      return(tangents)
    }

    for (i in seq_len(n)) {
      step <- if (i == 1) {
        coords[2, ] - coords[1, ]
      } else if (i == n) {
        coords[n, ] - coords[n - 1, ]
      } else {
        coords[i + 1, ] - coords[i - 1, ]
      }

      step_norm <- sqrt(sum(step^2))
      if (is.finite(step_norm) && step_norm > 0) {
        tangents[i, ] <- step / step_norm
      }
    }

    tangents
  }

  mean_trajectories <- list()
  tangent_vectors <- list()

  width_results <- do.call(rbind, lapply(all_labels, function(lbl) {
    lbl_data <- x[x$label == lbl, ]
    renditions <- unique(lbl_data$rendition)
    n_rend <- length(renditions)
    all_times <- sort(unique(lbl_data$.time))

    time_coverage <- vapply(all_times, function(t) {
      sum(lbl_data$.time == t) / n_rend
    }, numeric(1))
    common_times <- all_times[time_coverage >= min_coverage]

    if (length(common_times) < 2) {
      warning(sprintf("Label '%s': too few shared time steps; skipping", lbl))
      return(NULL)
    }

    mean_traj <- do.call(rbind, lapply(common_times, function(t) {
      t_vals <- lbl_data[lbl_data$.time == t, dims, drop = FALSE]
      vapply(dims, function(d) {
        mean(t_vals[[d]], trim = trim_fraction, na.rm = TRUE)
      }, numeric(1))
    }))
    colnames(mean_traj) <- dims

    tangents <- compute_tangents(mean_traj)
    colnames(tangents) <- dims

    mean_trajectories[[lbl]] <<- data.frame(.time = common_times, mean_traj)
    tangent_vectors[[lbl]] <<- data.frame(.time = common_times, tangents)

    reference_df <- data.frame(.time = common_times, mean_traj, tangents)
    mean_cols <- paste0(dims, "_mean")
    tangent_cols <- paste0(dims, "_tan")
    names(reference_df) <- c(".time", mean_cols, tangent_cols)

    if (verbose) {
      message(sprintf(
        "  %s: %d renditions | reference = %d time steps",
        lbl, n_rend, nrow(reference_df)
      ))
    }

    do.call(rbind, lapply(renditions, function(r) {
      rend_data <- lbl_data[lbl_data$rendition == r, ]
      rend_data <- rend_data[order(rend_data$.time), c("label", "rendition", ".time", dims)]
      aligned <- merge(rend_data, reference_df, by = ".time")

      if (nrow(aligned) < 2) {
        return(NULL)
      }

      residual_mat <- as.matrix(aligned[, dims, drop = FALSE]) -
        as.matrix(aligned[, mean_cols, drop = FALSE])
      tangent_mat <- as.matrix(aligned[, tangent_cols, drop = FALSE])

      total_sq <- rowSums(residual_mat^2)
      parallel <- rowSums(residual_mat * tangent_mat)
      parallel_sq <- parallel^2
      parallel_sq[rowSums(is.na(tangent_mat)) > 0] <- NA_real_
      orth_sq <- pmax(total_sq - parallel_sq, 0)

      data.frame(
        label = lbl,
        rendition = r,
        n_time = nrow(aligned),
        coverage = nrow(aligned) / nrow(reference_df),
        total_rms = sqrt(mean(total_sq, na.rm = TRUE)),
        orthogonal_rms = sqrt(mean(orth_sq, na.rm = TRUE)),
        parallel_rms = sqrt(mean(parallel_sq, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))
  }))

  if (is.null(width_results) || nrow(width_results) == 0) {
    stop("No valid renditions available for width analysis")
  }

  summary_df <- width_results |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(
      n = dplyr::n(),
      total_rms_mean = mean(.data$total_rms, na.rm = TRUE),
      total_rms_sd = sd(.data$total_rms, na.rm = TRUE),
      orthogonal_rms_mean = mean(.data$orthogonal_rms, na.rm = TRUE),
      orthogonal_rms_sd = sd(.data$orthogonal_rms, na.rm = TRUE),
      parallel_rms_mean = mean(.data$parallel_rms, na.rm = TRUE),
      parallel_rms_sd = sd(.data$parallel_rms, na.rm = TRUE),
      mean_coverage = mean(.data$coverage, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()

  tests <- NULL

  if (stats && length(unique(width_results$label)) > 1) {
    test_total <- stats::kruskal.test(total_rms ~ label, data = width_results)
    test_orth <- stats::kruskal.test(orthogonal_rms ~ label, data = width_results)
    test_parallel <- stats::kruskal.test(parallel_rms ~ label, data = width_results)

    posthoc_total <- stats::pairwise.wilcox.test(
      width_results$total_rms,
      width_results$label,
      p.adjust.method = "bonferroni"
    )
    posthoc_orth <- stats::pairwise.wilcox.test(
      width_results$orthogonal_rms,
      width_results$label,
      p.adjust.method = "bonferroni"
    )
    posthoc_parallel <- stats::pairwise.wilcox.test(
      width_results$parallel_rms,
      width_results$label,
      p.adjust.method = "bonferroni"
    )

    tests <- list(
      kruskal = list(total = test_total, orthogonal = test_orth, parallel = test_parallel),
      posthoc = list(total = posthoc_total, orthogonal = posthoc_orth, parallel = posthoc_parallel)
    )

    if (verbose) {
      message("\n--- Summary ---")
      print(summary_df)
      message("\nKruskal-Wallis tests:")
      message(sprintf(
        "  Total RMS:      chi-sq = %.2f, p = %.2e",
        test_total$statistic, test_total$p.value
      ))
      message(sprintf(
        "  Orthogonal RMS: chi-sq = %.2f, p = %.2e",
        test_orth$statistic, test_orth$p.value
      ))
      message(sprintf(
        "  Parallel RMS:   chi-sq = %.2f, p = %.2e",
        test_parallel$statistic, test_parallel$p.value
      ))
      message("\nPairwise Wilcoxon post-hoc (Bonferroni adjusted):")
      message("  Total RMS:")
      print(posthoc_total$p.value)
      message("  Orthogonal RMS:")
      print(posthoc_orth$p.value)
      message("  Parallel RMS:")
      print(posthoc_parallel$p.value)
    }
  } else if (verbose) {
    message("\n--- Summary ---")
    print(summary_df)
  }

  invisible(list(
    type              = "path_deviation",
    dims              = dims,
    min_coverage      = min_coverage,
    time_digits       = time_digits,
    width             = width_results,
    summary           = summary_df,
    mean_trajectories = mean_trajectories,
    tangent_vectors   = tangent_vectors,
    tests             = tests
  ))
}


#' @rdname trajectory_path_deviation
#' @export
trajectory_path_deviation.Sap <- function(x,
                                          segment_type = c(
                                            "motifs", "syllables",
                                            "bouts", "segments"
                                          ),
                                          dims = c("PC1", "PC2"),
                                          trim_fraction = 0.1,
                                          min_coverage = 0.5,
                                          time_digits = 6,
                                          labels = NULL,
                                          stats = TRUE,
                                          verbose = TRUE,
                                          ...) {
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  traj_embeds <- x$features[[feature_type]][["traj.embeds"]]
  if (is.null(traj_embeds)) {
    stop(sprintf(
      "Trajectory embeddings not found for %s. Run create_trajectory_matrix() first.",
      segment_type
    ))
  }

  missing_dims <- setdiff(dims, names(traj_embeds))
  if (length(missing_dims) > 0) {
    stop(sprintf(
      "Dimensions not found in traj.embeds: %s\nAvailable: %s",
      paste(missing_dims, collapse = ", "),
      paste(
        setdiff(
          names(traj_embeds),
          c(
            "filename", "day_post_hatch", "label", "rendition",
            "selec", "start_time", "end_time", ".time"
          )
        ),
        collapse = ", "
      )
    ))
  }

  result <- trajectory_path_deviation.default(
    x = traj_embeds,
    dims = dims,
    trim_fraction = trim_fraction,
    min_coverage = min_coverage,
    time_digits = time_digits,
    labels = labels,
    stats = stats,
    verbose = verbose,
    ...
  )

  # Write back to SAP object
  x$features[[feature_type]][["trajectory_path_deviation"]] <- result

  invisible(x)
}



# Trajectory UMAP Occupancy
# Update date : May 2, 2026

#' Trajectory UMAP Occupancy Analysis
#'
#' @description
#' Quantifies rendition-to-rendition diversity in local UMAP occupancy using a
#' shared 2D grid and same-label neighborhood structure. This is intended as a
#' complementary analysis to trajectory width in PCA space when the biological
#' question is about local latent-state exploration rather than calibrated
#' geometric variance.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or SAP object
#' @param dims Character vector of UMAP columns to use (default: c("UMAP1", "UMAP2"))
#' @param grid_n Number of bins per UMAP axis for occupancy calculations
#'   (default: 40)
#' @param k Number of same-label nearest neighbors used for local dispersion
#'   (default: 15)
#' @param peripheral_quantile Quantile of global occupied-bin density used to
#'   define peripheral bins (default: 0.2)
#' @param labels Optional character vector of labels to include
#' @param stats Logical. If \code{TRUE} (default), run Kruskal-Wallis and pairwise
#'   Wilcoxon tests. Set to \code{FALSE} to skip statistical testing and return
#'   \code{NULL} for \code{tests}.
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables',
#'   'bouts', 'segments')
#' @param ... Additional arguments
#'
#' @details
#' The function computes four per-rendition metrics in a shared UMAP grid:
#' \describe{
#'   \item{Occupied Fraction}{Fraction of grid bins visited by the rendition}
#'   \item{Occupancy Entropy}{Shannon entropy of the rendition's occupancy over
#'     grid bins, normalized to \code{[0, 1]} using the full shared grid}
#'   \item{Peripheral Fraction}{Fraction of the rendition's points falling in
#'     globally low-density bins of the shared UMAP manifold}
#'   \item{Same-Label kNN Dispersion}{Average distance from each point to its
#'     same-label nearest neighbors in UMAP space}
#' }
#'
#' @return A list (returned invisibly) with the following elements:
#' \itemize{
#'   \item \code{type}: Character string \code{"umap_occupancy"}, used by
#'     \code{plot_trajectory_variability()} for dispatch
#'   \item \code{occupancy}: Per-rendition occupancy metrics
#'   \item \code{summary}: Summary table with mean and SD for each metric per label
#'   \item \code{annotated_points}: Original data with occupancy annotations
#'   \item \code{bin_counts}: Shared-grid counts per label and bin
#'   \item \code{grid_info}: Grid settings and peripheral threshold metadata
#'   \item \code{tests}: Kruskal-Wallis and pairwise Wilcoxon tests (\code{NULL}
#'     when \code{stats = FALSE} or only one label is present)
#' }
#'
#' Use \code{\link{plot_trajectory_variability}(result)} to visualise the output.
#'
#' @examples
#' \dontrun{
#' result <- trajectory_umap_occupancy(sap)
#' result$summary
#' head(result$occupancy)
#' }
#'
#' @export
trajectory_umap_occupancy <- function(x, ...) {
  UseMethod("trajectory_umap_occupancy")
}


#' @rdname trajectory_umap_occupancy
#' @export
trajectory_umap_occupancy.default <- function(x,
                                              dims = c("UMAP1", "UMAP2"),
                                              grid_n = 40,
                                              k = 15,
                                              peripheral_quantile = 0.2,
                                              labels = NULL,
                                              stats = TRUE,
                                              verbose = TRUE,
                                              ...) {
  if (!is.data.frame(x)) stop("Input must be a data frame")
  if (length(dims) != 2) stop("Use exactly two dimensions for 2D UMAP occupancy")
  if (grid_n < 5) stop("grid_n must be at least 5")
  if (k < 1) stop("k must be at least 1")

  required_cols <- c("label", "rendition")
  missing_cols <- setdiff(c(required_cols, dims), names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  ensure_pkgs("ggplot2", "dplyr", "RANN")

  if (!is.null(labels)) {
    missing_labels <- setdiff(labels, unique(x$label))
    if (length(missing_labels) > 0) {
      stop(sprintf("Labels not found: %s", paste(missing_labels, collapse = ", ")))
    }
    x <- x[x$label %in% labels, ]
  }

  x <- x[stats::complete.cases(x[, dims, drop = FALSE]), , drop = FALSE]
  all_labels <- unique(x$label)
  if (length(all_labels) == 0) stop("No valid rows available after filtering")

  # Auto-disable stats when there are many groups
  if (stats && length(all_labels) > 6) {
    message(sprintf(
      "Note: %d labels detected. Setting stats = FALSE (statistical tests are not meaningful for this many groups).",
      length(all_labels)
    ))
    stats <- FALSE
  }

  if (verbose) {
    message("\n=== Trajectory UMAP Occupancy Analysis ===")
    message(sprintf("Dimensions          : %s", paste(dims, collapse = ", ")))
    message(sprintf("Grid size           : %d x %d", grid_n, grid_n))
    message(sprintf("Same-label kNN      : %d", k))
    message(sprintf("Peripheral quantile : %.0f%%", peripheral_quantile * 100))
    message(sprintf("Labels              : %s\n", paste(all_labels, collapse = ", ")))
  }

  x_min <- min(x[[dims[1]]], na.rm = TRUE)
  x_max <- max(x[[dims[1]]], na.rm = TRUE)
  y_min <- min(x[[dims[2]]], na.rm = TRUE)
  y_max <- max(x[[dims[2]]], na.rm = TRUE)

  if (x_min == x_max) {
    x_min <- x_min - 0.5
    x_max <- x_max + 0.5
  }
  if (y_min == y_max) {
    y_min <- y_min - 0.5
    y_max <- y_max + 0.5
  }

  x_breaks <- seq(x_min, x_max, length.out = grid_n + 1)
  y_breaks <- seq(y_min, y_max, length.out = grid_n + 1)
  bin_area <- diff(x_breaks)[1] * diff(y_breaks)[1]
  total_bins <- grid_n * grid_n

  x$bin_x <- cut(x[[dims[1]]],
    breaks = x_breaks, include.lowest = TRUE,
    labels = FALSE
  )
  x$bin_y <- cut(x[[dims[2]]],
    breaks = y_breaks, include.lowest = TRUE,
    labels = FALSE
  )
  x$bin_id <- paste(x$bin_x, x$bin_y, sep = "_")

  global_bin_counts <- sort(table(x$bin_id), decreasing = TRUE)
  peripheral_cut <- as.numeric(stats::quantile(
    as.numeric(global_bin_counts),
    peripheral_quantile,
    type = 1,
    names = FALSE
  ))
  peripheral_bins <- names(global_bin_counts)[global_bin_counts <= peripheral_cut]
  x$is_peripheral <- x$bin_id %in% peripheral_bins

  x$knn_same_label <- NA_real_
  for (lbl in all_labels) {
    idx <- which(x$label == lbl)
    coords <- as.matrix(x[idx, dims, drop = FALSE])
    n_lbl <- nrow(coords)
    if (n_lbl < 2) next

    k_eff <- min(k + 1, n_lbl)
    nn <- RANN::nn2(coords, k = k_eff)
    x$knn_same_label[idx] <- rowMeans(nn$nn.dists[, -1, drop = FALSE], na.rm = TRUE)
  }

  occupancy_results <- do.call(rbind, lapply(all_labels, function(lbl) {
    lbl_data <- x[x$label == lbl, , drop = FALSE]
    renditions <- unique(lbl_data$rendition)

    do.call(rbind, lapply(renditions, function(r) {
      rend_data <- lbl_data[lbl_data$rendition == r, , drop = FALSE]
      counts <- table(rend_data$bin_id)
      probs <- as.numeric(counts) / sum(counts)
      entropy <- -sum(probs * log(probs))

      data.frame(
        label = lbl,
        rendition = r,
        n_points = nrow(rend_data),
        occupied_bins = length(counts),
        occupied_fraction = length(counts) / total_bins,
        occupied_area = length(counts) * bin_area,
        occupancy_entropy = entropy,
        occupancy_entropy_norm = if (total_bins > 1) entropy / log(total_bins) else NA_real_,
        peripheral_fraction = mean(rend_data$is_peripheral, na.rm = TRUE),
        knn_dispersion = mean(rend_data$knn_same_label, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))

  summary_df <- occupancy_results |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(
      n = dplyr::n(),
      occupied_fraction_mean = mean(.data$occupied_fraction, na.rm = TRUE),
      occupied_fraction_sd = sd(.data$occupied_fraction, na.rm = TRUE),
      occupancy_entropy_norm_mean = mean(.data$occupancy_entropy_norm, na.rm = TRUE),
      occupancy_entropy_norm_sd = sd(.data$occupancy_entropy_norm, na.rm = TRUE),
      peripheral_fraction_mean = mean(.data$peripheral_fraction, na.rm = TRUE),
      peripheral_fraction_sd = sd(.data$peripheral_fraction, na.rm = TRUE),
      knn_dispersion_mean = mean(.data$knn_dispersion, na.rm = TRUE),
      knn_dispersion_sd = sd(.data$knn_dispersion, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()

  bin_counts <- x |>
    dplyr::group_by(.data$label, .data$bin_x, .data$bin_y, .data$bin_id) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    as.data.frame()

  tests <- NULL
  if (stats && length(unique(occupancy_results$label)) > 1) {
    test_occ <- stats::kruskal.test(occupied_fraction ~ label, data = occupancy_results)
    test_ent <- stats::kruskal.test(occupancy_entropy_norm ~ label, data = occupancy_results)
    test_per <- stats::kruskal.test(peripheral_fraction ~ label, data = occupancy_results)
    test_knn <- stats::kruskal.test(knn_dispersion ~ label, data = occupancy_results)

    posthoc_occ <- stats::pairwise.wilcox.test(
      occupancy_results$occupied_fraction,
      occupancy_results$label,
      p.adjust.method = "bonferroni",
      exact = FALSE
    )
    posthoc_ent <- stats::pairwise.wilcox.test(
      occupancy_results$occupancy_entropy_norm,
      occupancy_results$label,
      p.adjust.method = "bonferroni",
      exact = FALSE
    )
    posthoc_per <- stats::pairwise.wilcox.test(
      occupancy_results$peripheral_fraction,
      occupancy_results$label,
      p.adjust.method = "bonferroni",
      exact = FALSE
    )
    posthoc_knn <- stats::pairwise.wilcox.test(
      occupancy_results$knn_dispersion,
      occupancy_results$label,
      p.adjust.method = "bonferroni",
      exact = FALSE
    )

    tests <- list(
      kruskal = list(
        occupied_fraction = test_occ,
        entropy = test_ent,
        peripheral_fraction = test_per,
        knn_dispersion = test_knn
      ),
      posthoc = list(
        occupied_fraction = posthoc_occ,
        entropy = posthoc_ent,
        peripheral_fraction = posthoc_per,
        knn_dispersion = posthoc_knn
      )
    )

    if (verbose) {
      message("\n--- Summary ---")
      print(summary_df)
      message("\nKruskal-Wallis tests:")
      message(sprintf(
        "  Occupied fraction:  chi-sq = %.2f, p = %.2e",
        test_occ$statistic, test_occ$p.value
      ))
      message(sprintf(
        "  Occupancy entropy:  chi-sq = %.2f, p = %.2e",
        test_ent$statistic, test_ent$p.value
      ))
      message(sprintf(
        "  Peripheral usage:   chi-sq = %.2f, p = %.2e",
        test_per$statistic, test_per$p.value
      ))
      message(sprintf(
        "  kNN dispersion:     chi-sq = %.2f, p = %.2e",
        test_knn$statistic, test_knn$p.value
      ))
      message("\nPairwise Wilcoxon post-hoc (Bonferroni adjusted):")
      message("  Occupied fraction:")
      print(posthoc_occ$p.value)
      message("  Occupancy entropy:")
      print(posthoc_ent$p.value)
      message("  Peripheral usage:")
      print(posthoc_per$p.value)
      message("  kNN dispersion:")
      print(posthoc_knn$p.value)
    }
  } else if (verbose) {
    message("\n--- Summary ---")
    print(summary_df)
  }

  invisible(list(
    type = "umap_occupancy",
    dims = dims,
    occupancy = occupancy_results,
    summary = summary_df,
    annotated_points = x,
    bin_counts = bin_counts,
    grid_info = list(
      dims                = dims,
      grid_n              = grid_n,
      total_bins          = total_bins,
      x_breaks            = x_breaks,
      y_breaks            = y_breaks,
      bin_area            = bin_area,
      peripheral_quantile = peripheral_quantile,
      peripheral_cut      = peripheral_cut,
      peripheral_bins     = peripheral_bins
    ),
    tests = tests
  ))
}


#' @rdname trajectory_umap_occupancy
#' @export
trajectory_umap_occupancy.Sap <- function(x,
                                          segment_type = c(
                                            "motifs", "syllables",
                                            "bouts", "segments"
                                          ),
                                          dims = c("UMAP1", "UMAP2"),
                                          grid_n = 40,
                                          k = 15,
                                          peripheral_quantile = 0.2,
                                          labels = NULL,
                                          stats = TRUE,
                                          verbose = TRUE,
                                          ...) {
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  traj_embeds <- x$features[[feature_type]][["traj.embeds"]]
  if (is.null(traj_embeds)) {
    stop(sprintf(
      "Trajectory embeddings not found for %s. Run create_trajectory_matrix() first.",
      segment_type
    ))
  }

  missing_dims <- setdiff(dims, names(traj_embeds))
  if (length(missing_dims) > 0) {
    stop(sprintf(
      "Dimensions not found in traj.embeds: %s\nAvailable: %s",
      paste(missing_dims, collapse = ", "),
      paste(
        setdiff(
          names(traj_embeds),
          c(
            "filename", "day_post_hatch", "label", "rendition",
            "selec", "start_time", "end_time", ".time"
          )
        ),
        collapse = ", "
      )
    ))
  }

  result <- trajectory_umap_occupancy.default(
    x = traj_embeds,
    dims = dims,
    grid_n = grid_n,
    k = k,
    peripheral_quantile = peripheral_quantile,
    labels = labels,
    stats = stats,
    verbose = verbose,
    ...
  )

  # Write back to SAP object
  x$features[[feature_type]][["trajectory_umap_occupancy"]] <- result

  invisible(x)
}


# Unified Trajectory Variability Plot
# Update date : Jun. 22, 2026

#' Plot Trajectory Variability Results
#'
#' @description
#' A unified plotting function for results produced by
#' \code{\link{trajectory_dispersion}}, \code{\link{trajectory_path_deviation}},
#' or \code{\link{trajectory_umap_occupancy}}. Dispatches to the appropriate
#' panel layout based on the \code{type} field embedded in \code{result} or
#' retrieved from the SAP object and draws significance brackets when statistical
#' tests are present.
#'
#' @param x A list returned by \code{trajectory_dispersion()},
#'   \code{trajectory_path_deviation()}, or
#'   \code{trajectory_umap_occupancy()} (for the default method); or a SAP object
#'   (for the Sap method).
#' @param palette RColorBrewer palette name (default: \code{"Set1"}). When the
#'   number of labels exceeds the palette's maximum, colours are interpolated
#'   automatically via \code{colorRampPalette}.
#' @param max_annotations Maximum number of pairwise significance brackets to
#'   draw per panel (default: \code{10}). When more comparisons exist, the
#'   most significant pairs are retained and a message is issued.
#' @param show_cv Logical. If \code{TRUE}, add coefficient-of-variation panels
#'   for centroid dispersion and trajectory path length when plotting
#'   \code{trajectory_dispersion()} results (default: \code{FALSE}).
#' @param segment_type For SAP objects: Type of segments to visualize ('motifs', 'syllables', 'bouts', 'segments')
#' @param variability_type For SAP objects: Which computed variability type to plot ('dispersion', 'path_deviation', 'umap_occupancy')
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' Panel layouts by result type:
#' \describe{
#'   \item{\code{"dispersion"}}{3 panels: Mean Pairwise Distance \u00b7 Centroid
#'     Dispersion \u00b7 Path Length}
#'   \item{\code{"path_deviation"}}{3 panels: Total RMS \u00b7 Orthogonal RMS \u00b7
#'     Parallel RMS}
#'   \item{\code{"umap_occupancy"}}{4 panels: Occupied Fraction \u00b7 Occupancy
#'     Entropy \u00b7 Peripheral Fraction \u00b7 kNN Dispersion}
#' }
#' When \code{show_cv = TRUE} for \code{"dispersion"} results, two additional
#' panels show SD divided by mean for centroid dispersion and path length.
#'
#' Each panel displays a violin + box plot coloured by label. When
#' statistical tests are not \code{NULL}, Kruskal-Wallis p-values are shown
#' as subtitles and pairwise Wilcoxon p-values appear as significance brackets
#' above the data.
#'
#' @return The assembled \pkg{patchwork} object, printed as a side-effect and
#'   returned invisibly so the caller can save or further modify it.
#'
#' @examples
#' \dontrun{
#' # Plotting directly from a result list
#' result <- trajectory_dispersion(sap$features$motif$traj.embeds, dims = c("PC1", "PC2"))
#' plot_trajectory_variability(result)
#'
#' # Plotting from a SAP object with pre-computed results
#' sap <- trajectory_dispersion(sap)
#' plot_trajectory_variability(sap, segment_type = "motifs", variability_type = "dispersion")
#' }
#'
#' @export
plot_trajectory_variability <- function(x, ...) {
  UseMethod("plot_trajectory_variability")
}


#' @rdname plot_trajectory_variability
#' @export
plot_trajectory_variability.default <- function(x,
                                                palette = "Set1",
                                                max_annotations = 10,
                                                show_cv = FALSE,
                                                ...) {
  result <- x
  # ---- Validate input ----
  if (!is.list(result) || is.null(result$type)) {
    stop(paste(
      "'result' must be a list returned by trajectory_dispersion(),",
      "trajectory_path_deviation(), or trajectory_umap_occupancy().",
      "It must contain a 'type' element."
    ))
  }
  valid_types <- c("dispersion", "path_deviation", "umap_occupancy")
  if (!result$type %in% valid_types) {
    stop(sprintf(
      "Unknown result type '%s'. Expected one of: %s",
      result$type, paste(valid_types, collapse = ", ")
    ))
  }

  ensure_pkgs("ggplot2", "patchwork")

  dims <- result$dims
  many_labs <- FALSE

  # ---- Dispatch on type ----
  if (result$type == "dispersion") {
    # ---- trajectory_dispersion ----
    pw <- result$pairwise
    dis <- result$dispersion
    pl <- result$path_length
    tst <- result$tests

    labs_order <- sort_labels(unique(as.character(pw$label)))
    many_labs <- length(labs_order) > 6
    pal_map <- make_pal(labs_order, palette)

    kw_pw <- if (!is.null(tst)) fmt_p(tst$kruskal$pairwise$p.value) else NULL
    kw_dis <- if (!is.null(tst)) fmt_p(tst$kruskal$dispersion$p.value) else NULL
    kw_pl <- if (!is.null(tst)) fmt_p(tst$kruskal$path_length$p.value) else NULL

    if (many_labs) {
      p1 <- trend_panel(pw, "label", "mean_dist", "Mean Pairwise Distance", kw_pw, labs_order)
      p2 <- trend_panel(dis, "label", "dispersion", "Centroid Dispersion", kw_dis, labs_order)
      p3 <- trend_panel(pl, "label", "path_length", "Trajectory Path Length", kw_pl, labs_order)
    } else {
      p1 <- panel(pw, "label", "mean_dist", "Mean Pairwise Distance", kw_pw, pal_map, labs_order)
      p2 <- panel(dis, "label", "dispersion", "Centroid Dispersion", kw_dis, pal_map, labs_order)
      p3 <- panel(pl, "label", "path_length", "Trajectory Path Length", kw_pl, pal_map, labs_order)
      if (!is.null(tst)) {
        p1 <- add_brackets(p1, brackets(pw$mean_dist, tst$posthoc$pairwise, labs_order, max_annotations))
        p2 <- add_brackets(p2, brackets(dis$dispersion, tst$posthoc$dispersion, labs_order, max_annotations))
        p3 <- add_brackets(p3, brackets(pl$path_length, tst$posthoc$path_length, labs_order, max_annotations))
      }
    }

    if (show_cv) {
      cv <- trajectory_cv_summary(result$summary, labs_order)
      # Create CV panels with black line connections and extra y-axis padding
      p4 <- panel(cv, "label", "dispersion_cv", NULL, NULL, pal_map, labs_order, 
                  y_label = "Dispersion CV") +
        ggplot2::geom_line(ggplot2::aes(group = 1), color = "black", linewidth = 0.5) +
        ggplot2::geom_point(color = "black", size = 2) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1)))
      p5 <- panel(cv, "label", "path_length_cv", NULL, NULL, pal_map, labs_order,
                  y_label = "Path Length CV") +
        ggplot2::geom_line(ggplot2::aes(group = 1), color = "black", linewidth = 0.5) +
        ggplot2::geom_point(color = "black", size = 2) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1)))
      # Layout: Row 1 = all three main panels, Row 2 = empty + two CV panels
      panels <- (p1 + p2 + p3) / (patchwork::plot_spacer() + p4 + p5)
    } else {
      panels <- p1 + p2 + p3
    }

    combined <- panels +
      patchwork::plot_annotation(
        title = "Trajectory Dispersion Comparison",
        subtitle = paste("Dimensions:", paste(dims, collapse = " + ")),
        theme = ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = "white", color = NA)
        )
      )
  } else if (result$type == "path_deviation") {
    # ---- trajectory_path_deviation ----
    wd <- result$width
    tst <- result$tests

    labs_order <- sort_labels(unique(as.character(wd$label)))
    many_labs <- length(labs_order) > 6
    pal_map <- make_pal(labs_order, palette)

    kw_tot <- if (!is.null(tst)) fmt_p(tst$kruskal$total$p.value) else NULL
    kw_orth <- if (!is.null(tst)) fmt_p(tst$kruskal$orthogonal$p.value) else NULL
    kw_par <- if (!is.null(tst)) fmt_p(tst$kruskal$parallel$p.value) else NULL

    if (many_labs) {
      p1 <- trend_panel(wd, "label", "total_rms", "Total RMS Residual", kw_tot, labs_order)
      p2 <- trend_panel(wd, "label", "orthogonal_rms", "Orthogonal RMS Residual", kw_orth, labs_order)
      p3 <- trend_panel(wd, "label", "parallel_rms", "Parallel RMS Residual", kw_par, labs_order)
    } else {
      p1 <- panel(wd, "label", "total_rms", "Total RMS Residual", kw_tot, pal_map, labs_order)
      p2 <- panel(wd, "label", "orthogonal_rms", "Orthogonal RMS Residual", kw_orth, pal_map, labs_order)
      p3 <- panel(wd, "label", "parallel_rms", "Parallel RMS Residual", kw_par, pal_map, labs_order)
      if (!is.null(tst)) {
        p1 <- add_brackets(p1, brackets(wd$total_rms, tst$posthoc$total, labs_order, max_annotations))
        p2 <- add_brackets(p2, brackets(wd$orthogonal_rms, tst$posthoc$orthogonal, labs_order, max_annotations))
        p3 <- add_brackets(p3, brackets(wd$parallel_rms, tst$posthoc$parallel, labs_order, max_annotations))
      }
    }

    combined <- (p1 + p2 + p3) +
      patchwork::plot_annotation(
        title = "Trajectory Path Deviation Comparison",
        subtitle = paste("Dimensions:", paste(dims, collapse = " + ")),
        theme = ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = "white", color = NA)
        )
      )
  } else {
    # ---- trajectory_umap_occupancy ----
    occ <- result$occupancy
    tst <- result$tests

    labs_order <- sort_labels(unique(as.character(occ$label)))
    many_labs <- length(labs_order) > 6
    pal_map <- make_pal(labs_order, palette)

    kw_occ <- if (!is.null(tst)) fmt_p(tst$kruskal$occupied_fraction$p.value) else NULL
    kw_ent <- if (!is.null(tst)) fmt_p(tst$kruskal$entropy$p.value) else NULL
    kw_per <- if (!is.null(tst)) fmt_p(tst$kruskal$peripheral_fraction$p.value) else NULL
    kw_knn <- if (!is.null(tst)) fmt_p(tst$kruskal$knn_dispersion$p.value) else NULL

    if (many_labs) {
      p1 <- trend_panel(occ, "label", "occupied_fraction", "Occupied Fraction", kw_occ, labs_order)
      p2 <- trend_panel(occ, "label", "occupancy_entropy_norm", "Occupancy Entropy", kw_ent, labs_order)
      p3 <- trend_panel(occ, "label", "peripheral_fraction", "Peripheral Fraction", kw_per, labs_order)
      p4 <- trend_panel(occ, "label", "knn_dispersion", "Same-Label kNN Dispersion", kw_knn, labs_order)
    } else {
      p1 <- panel(occ, "label", "occupied_fraction", "Occupied Fraction", kw_occ, pal_map, labs_order, jitter = TRUE)
      p2 <- panel(occ, "label", "occupancy_entropy_norm", "Occupancy Entropy", kw_ent, pal_map, labs_order, jitter = TRUE)
      p3 <- panel(occ, "label", "peripheral_fraction", "Peripheral Fraction", kw_per, pal_map, labs_order, jitter = TRUE)
      p4 <- panel(occ, "label", "knn_dispersion", "Same-Label kNN Dispersion", kw_knn, pal_map, labs_order, jitter = TRUE)
      if (!is.null(tst)) {
        p1 <- add_brackets(p1, brackets(occ$occupied_fraction, tst$posthoc$occupied_fraction, labs_order, max_annotations))
        p2 <- add_brackets(p2, brackets(occ$occupancy_entropy_norm, tst$posthoc$entropy, labs_order, max_annotations))
        p3 <- add_brackets(p3, brackets(occ$peripheral_fraction, tst$posthoc$peripheral_fraction, labs_order, max_annotations))
        p4 <- add_brackets(p4, brackets(occ$knn_dispersion, tst$posthoc$knn_dispersion, labs_order, max_annotations))
      }
    }

    combined <- (p1 + p2 + p3 + p4) +
      patchwork::plot_annotation(
        title = "Trajectory UMAP Occupancy Comparison",
        subtitle = paste("Dimensions:", paste(dims, collapse = " + ")),
        theme = ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = "white", color = NA)
        )
      )
  }

  withCallingHandlers(
    print(combined),
    warning = function(w) {
      if (grepl("annotation$theme is not a valid theme", w$message, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  invisible(combined)
}


#' @rdname plot_trajectory_variability
#' @export
plot_trajectory_variability.Sap <- function(x,
                                            segment_type = c("motifs", "syllables", "bouts", "segments"),
                                            variability_type = c("dispersion", "path_deviation", "umap_occupancy"),
                                            palette = "Set1",
                                            max_annotations = 10,
                                            show_cv = FALSE,
                                            ...) {
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")
  segment_type <- match.arg(segment_type)
  variability_type <- match.arg(variability_type)

  feature_type <- sub("s$", "", segment_type)
  slot_name <- paste0("trajectory_", variability_type)

  result <- x$features[[feature_type]][[slot_name]]
  if (is.null(result)) {
    stop(sprintf(
      "No %s results found for %s. Run trajectory_%s() first.",
      variability_type,
      segment_type,
      variability_type
    ))
  }

  plot_trajectory_variability.default(
    x = result,
    palette = palette,
    max_annotations = max_annotations,
    show_cv = show_cv,
    ...
  )
}


#' Bin Trajectory Time Values
#'
#' @description
#' Internal helper used by trajectory analyses to make nominally identical
#' sliding-window time points match exactly across renditions.
#'
#' @keywords internal
#' @noRd
bin_trajectory_time_data <- function(x, dims, time_digits = 6) {
  if (!is.numeric(time_digits) || length(time_digits) != 1 || is.na(time_digits) ||
      time_digits < 0 || time_digits != floor(time_digits)) {
    stop("time_digits must be a non-negative whole number")
  }

  x$.time <- round(as.numeric(x$.time), digits = time_digits)

  x |>
    dplyr::group_by(.data$label, .data$rendition, .data$.time) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(dims), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$label, .data$rendition, .data$.time) |>
    as.data.frame()
}


#' Sort Trajectory Labels
#'
#' @description
#' Internal helper used by trajectory plotting and convergence functions to
#' sort labels numerically when possible and alphabetically otherwise.
#'
#' @keywords internal
#' @noRd
sort_labels <- function(labels) {
  nums <- suppressWarnings(as.numeric(labels))
  if (!anyNA(nums)) labels[order(nums)] else sort(labels)
}


#' Make Trajectory Plot Palette
#'
#' @description
#' Internal helper used by trajectory plotting functions to expand a
#' RColorBrewer palette to any number of labels.
#'
#' @keywords internal
#' @noRd
make_pal <- function(labels, palette) {
  n <- length(labels)
  pal_info <- RColorBrewer::brewer.pal.info
  max_col <- if (palette %in% rownames(pal_info)) {
    pal_info[palette, "maxcolors"]
  } else {
    8L
  }
  if (n <= max_col) {
    cols <- RColorBrewer::brewer.pal(max(n, 3L), palette)[seq_len(n)]
  } else {
    cols <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(max_col, palette)
    )(n)
  }
  stats::setNames(cols, labels)
}


#' Format Trajectory Plot P Values
#'
#' @description
#' Internal helper used by trajectory plotting functions to format p-values for
#' subtitles and pairwise brackets.
#'
#' @keywords internal
#' @noRd
fmt_p <- function(p) {
  if (is.null(p) || is.na(p)) {
    return("p = NA")
  }
  if (p < 0.001) sprintf("p = %.1e", p) else sprintf("p = %.3f", p)
}


#' Get Pairwise P Value
#'
#' @description
#' Internal helper used by trajectory plotting functions to look up a p-value
#' from a symmetric pairwise matrix.
#'
#' @keywords internal
#' @noRd
get_p <- function(pmat, g1, g2) {
  if (is.null(pmat)) {
    return(NA_real_)
  }
  rn <- rownames(pmat)
  cn <- colnames(pmat)
  if (g1 %in% rn && g2 %in% cn) {
    return(pmat[g1, g2])
  }
  if (g2 %in% rn && g1 %in% cn) {
    return(pmat[g2, g1])
  }
  NA_real_
}


#' Build Trajectory Significance Brackets
#'
#' @description
#' Internal helper used by trajectory plotting functions to build pairwise
#' significance bracket annotations.
#'
#' @keywords internal
#' @noRd
brackets <- function(values, posthoc_obj, labels_in_order, max_annotations) {
  comps <- utils::combn(labels_in_order, 2, simplify = FALSE)

  if (length(comps) > max_annotations) {
    pvals <- vapply(comps, function(comp) {
      get_p(posthoc_obj$p.value, comp[1], comp[2])
    }, numeric(1))
    keep_idx <- order(pvals)[seq_len(max_annotations)]
    message(sprintf(
      "  %d pairwise comparisons available; showing %d most significant.",
      length(comps), max_annotations
    ))
    comps <- comps[keep_idx]
  }

  y_range <- range(values, na.rm = TRUE)
  span <- diff(y_range)
  if (!is.finite(span) || span == 0) span <- max(abs(y_range), 1, na.rm = TRUE)
  base_y <- y_range[2] + 0.08 * span
  step_y <- 0.12 * span

  ann <- do.call(rbind, lapply(seq_along(comps), function(i) {
    comp <- comps[[i]]
    p_val <- get_p(posthoc_obj$p.value, comp[1], comp[2])
    data.frame(
      x1 = match(comp[1], labels_in_order),
      x2 = match(comp[2], labels_in_order),
      y = base_y + (i - 1) * step_y,
      lbl = fmt_p(p_val),
      stringsAsFactors = FALSE
    )
  }))
  ann$y_text <- ann$y + 0.025 * span
  ann$y_tip <- ann$y - 0.020 * span
  ann$y_max <- max(ann$y_text) + 0.06 * span
  ann
}


#' Add Trajectory Significance Brackets
#'
#' @description
#' Internal helper used by trajectory plotting functions to add bracket layers
#' to a ggplot object.
#'
#' @keywords internal
#' @noRd
add_brackets <- function(p, ann) {
  if (is.null(ann) || nrow(ann) == 0) {
    return(p)
  }
  p +
    ggplot2::geom_segment(
      data = ann,
      ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, yend = .data$y),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = ann,
      ggplot2::aes(x = .data$x1, xend = .data$x1, y = .data$y_tip, yend = .data$y),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = ann,
      ggplot2::aes(x = .data$x2, xend = .data$x2, y = .data$y_tip, yend = .data$y),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = ann,
      ggplot2::aes(
        x = (.data$x1 + .data$x2) / 2,
        y = .data$y_text,
        label = .data$lbl
      ),
      inherit.aes = FALSE, size = 3
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, ann$y_max[1]))
}


#' Build Trajectory Distribution Panel
#'
#' @description
#' Internal helper used by trajectory plotting functions to create a violin and
#' box plot panel.
#'
#' @keywords internal
#' @noRd
panel <- function(df, x_col, y_col, title, subtitle, pal_map,
                  labs_order, jitter = FALSE, y_label = y_col) {
  df[[x_col]] <- factor(df[[x_col]], levels = labs_order)
  n_labs <- length(labs_order)
  x_breaks <- if (n_labs >= 5) {
    labs_order[round(seq(1, n_labs, length.out = 5))]
  } else {
    labs_order
  }
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x    = .data[[x_col]],
      y    = .data[[y_col]],
      fill = .data[[x_col]]
    )
  ) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::labs(title = title, subtitle = subtitle, y = y_label, x = NULL) +
    ggplot2::scale_fill_manual(values = pal_map) +
    ggplot2::scale_x_discrete(breaks = x_breaks) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position  = "none",
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA)
    )
  if (jitter) {
    p <- p + ggplot2::geom_jitter(width = 0.12, alpha = 0.35, size = 0.9)
  }
  p
}


#' Build Trajectory Trend Panel
#'
#' @description
#' Internal helper used by trajectory plotting functions to create a mean and
#' standard-error trend panel.
#'
#' @keywords internal
#' @noRd
trend_panel <- function(df, x_col, y_col, title, subtitle, labs_order,
                        y_label = y_col) {
  agg <- do.call(rbind, lapply(labs_order, function(lbl) {
    vals <- df[[y_col]][as.character(df[[x_col]]) == lbl]
    n <- sum(!is.na(vals))
    data.frame(
      label = lbl,
      mean = mean(vals, na.rm = TRUE),
      se = if (n > 1) sd(vals, na.rm = TRUE) / sqrt(n) else 0,
      stringsAsFactors = FALSE
    )
  }))
  agg$label <- factor(agg$label, levels = labs_order)

  n_labs <- length(labs_order)
  x_breaks <- if (n_labs >= 5) {
    labs_order[round(seq(1, n_labs, length.out = 5))]
  } else {
    labs_order
  }

  ggplot2::ggplot(
    agg,
    ggplot2::aes(x = .data$label, y = .data$mean, group = 1)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$mean - .data$se,
        ymax = .data$mean + .data$se
      ),
      alpha = 0.2, fill = "steelblue"
    ) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.9) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      y        = y_label,
      x        = NULL
    ) +
    ggplot2::scale_x_discrete(breaks = x_breaks) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position  = "none",
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA)
    )
}


#' Build Trajectory CV Summary
#'
#' @description
#' Internal helper used by trajectory dispersion plotting to calculate
#' coefficient-of-variation values from the stored summary table.
#'
#' @keywords internal
#' @noRd
trajectory_cv_summary <- function(summary_df, labs_order) {
  required <- c(
    "label",
    "dispersion_mean",
    "dispersion_sd",
    "path_length_mean",
    "path_length_sd"
  )
  missing_cols <- setdiff(required, names(summary_df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Cannot plot CV panels; missing summary columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  summary_df <- summary_df[match(labs_order, as.character(summary_df$label)), ]
  data.frame(
    label = labs_order,
    dispersion_cv = ifelse(
      is.na(summary_df$dispersion_mean) | summary_df$dispersion_mean == 0,
      NA_real_,
      summary_df$dispersion_sd / summary_df$dispersion_mean
    ),
    path_length_cv = ifelse(
      is.na(summary_df$path_length_mean) | summary_df$path_length_mean == 0,
      NA_real_,
      summary_df$path_length_sd / summary_df$path_length_mean
    ),
    stringsAsFactors = FALSE
  )
}


#' Build Trajectory CV Panel
#'
#' @description
#' Internal helper used by trajectory dispersion plotting to create compact
#' coefficient-of-variation panels.
#'
#' @keywords internal
#' @noRd
cv_panel <- function(df, x_col, y_col, title, pal_map, labs_order) {
  df[[x_col]] <- factor(df[[x_col]], levels = labs_order)
  n_labs <- length(labs_order)
  x_breaks <- if (n_labs >= 5) {
    labs_order[round(seq(1, n_labs, length.out = 5))]
  } else {
    labs_order
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      fill = .data[[x_col]]
    )
  ) +
    ggplot2::geom_col(alpha = 0.75, width = 0.7, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::labs(title = title, y = "CV (SD / mean)", x = NULL) +
    ggplot2::scale_fill_manual(values = pal_map) +
    ggplot2::scale_x_discrete(breaks = x_breaks) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}


# Trajectory Convergence Analysis
# Update date : Jun. 23, 2026

#' Trajectory Convergence to Reference Path
#'
#' @description
#' Measures song maturity and developmental convergence by quantifying how
#' individual trial PC trajectories approach a reference path (e.g.,
#' adult/crystallized song). Uses complementary distance metrics (RMS,
#' Fr\'echet, DTW) with unified variability scaling and correlation analysis.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or SAP object
#' @param dims Character vector of dimension columns to use (default: c("PC1", "PC2"))
#' @param reference_label Character. Label to use as reference trajectory
#'   (default: NULL, uses last sorted label)
#' @param labels Optional character vector of labels to include
#' @param metrics Character vector of metrics to compute
#'   (default: c("rms", "frechet", "dtw", "correlation"))
#' @param interpolate_n Optional integer. If NULL (default), use matched time
#'   points. If integer, resample paths to interpolate_n equally spaced points
#' @param trim_fraction Numeric. Trim fraction for robust reference path
#'   building (default: 0.1, removes 10% from each tail)
#' @param min_coverage Minimum fraction of reference-label renditions that must
#'   cover a binned time step for it to contribute to the reference trajectory
#'   (default: 0.5).
#' @param time_digits Number of decimal places used to bin \code{.time} before
#'   grouping and matching trajectories (default: \code{6}).
#' @param similarity_baseline How distance similarities are transformed.
#'   \code{"reference"} (default) treats distances within the metric-specific
#'   reference-label median as converged before applying the exponential
#'   transform. \code{"zero"} uses the raw normalized distance.
#' @param similarity_scale_multiplier Multiplier applied to metric-specific
#'   reference scales before transforming distance similarities (default:
#'   \code{1}).
#' @param stats Logical. If \code{TRUE} (default), run Kruskal-Wallis and
#'   pairwise Wilcoxon tests on non-reference labels. Set to \code{FALSE} to
#'   skip statistical testing
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables',
#'   'bouts', 'segments')
#' @param ... Additional arguments
#'
#' @details
#' **Convergence Metrics:**
#' \itemize{
#'   \item RMS distance: Pointwise Euclidean distance, emphasizes persistent
#'     deviations across the full trajectory
#'   \item Fr\'echet distance: Curve shape distance; handles variable-length paths
#'   \item DTW distance: Dynamic time warping; captures shape similarity under
#'     timing shifts
#'   \item Correlation: Per-dimension Pearson r, averaged; intrinsically
#'     normalized (scale-free, range \eqn{[-1, 1]})
#' }
#'
#' **Reference Building:**
#' The reference trajectory is the robust trimmed-mean path of the reference
#' label. Trim fraction (default 0.1) removes outlier reference renditions at
#' each time point to avoid distortion. Time values are first rounded to
#' \code{time_digits}, and only binned time points covered by at least
#' \code{min_coverage} of reference renditions are retained.
#'
#' **Variability Scaling:**
#' Metric-specific reference scales are computed as the median reference-label
#' rendition distance to the reference path for RMS, Fr\'echet, and DTW. These
#' scales are used to normalize distance metrics (not correlation, which is
#' already scale-free).
#'
#' **Similarity Scores:**
#' Distance metrics are normalized by their metric-specific reference scale. By
#' default, the reference scale is multiplied by
#' \code{similarity_scale_multiplier = 1} before transformation. By
#' default, \code{similarity_baseline = "reference"} transforms excess distance
#' beyond the scaled reference threshold:
#' \code{similarity = exp(-pmax(normalized_distance - 1, 0))}.
#' With \code{similarity_baseline = "zero"}, the transform is:
#' \code{similarity = exp(-normalized_distance)}.
#'
#' Interpretation with the default reference baseline: similarity approximately
#' equals 1 means within the scaled reference threshold, similarity
#' approximately equals 0.37 means one additional scaled reference unit beyond
#' that threshold, and similarity less than 0.1 means far from the reference
#' path.
#'
#' **Statistical Testing:**
#' Reference label is included in results and plots but excluded from
#' Kruskal-Wallis and pairwise Wilcoxon tests (it defines the target, not a
#' convergence point).
#'
#' @return
#' For default method: A list (returned invisibly) with:
#' \itemize{
#'   \item \code{type}: Character string \code{"convergence"}
#'   \item \code{dims}: Requested dimensions
#'   \item \code{reference_label}: Reference label used
#'   \item \code{trim_fraction}: Trim fraction applied
#'   \item \code{min_coverage}: Minimum coverage threshold applied to the reference path
#'   \item \code{time_digits}: Decimal places used to bin \code{.time}
#'   \item \code{reference_scale}: Scalar variability scale
#'   \item \code{reference_scales}: Metric-specific variability scales
#'   \item \code{similarity_baseline}: Similarity transform baseline
#'   \item \code{similarity_scale_multiplier}: Similarity scale multiplier
#'   \item \code{reference_path}: Data frame of reference trajectory
#'   \item \code{metrics}: Metrics computed
#'   \item \code{interpolate_n}: Interpolation info
#'   \item \code{similarity}: Per-rendition results with all metric columns
#'   \item \code{summary}: Per-label summary statistics
#'   \item \code{tests}: Statistical test results (\code{NULL} if stats=FALSE
#'     or only one non-reference label)
#' }
#' For SAP objects: Updated object with results stored in
#'   \code{x$features[[feature_type]][["trajectory_convergence"]]}
#'
#' @examples
#' \dontrun{
#' # Compute convergence from trajectory embeddings
#' result <- trajectory_convergence(sap$features$motif$traj.embeds,
#'   dims = c("PC1", "PC2")
#' )
#'
#' # From SAP object
#' sap <- trajectory_convergence(sap)
#'
#' # Custom reference and metrics
#' result <- trajectory_convergence(sap,
#'   reference_label = "Adult",
#'   metrics = c("rms", "correlation")
#' )
#'
#' # Access results
#' result$summary # per-label statistics
#' result$tests # statistical tests
#' }
#'
#' @seealso \code{\link{plot_trajectory_convergence}} for visualization,
#'   \code{\link{trajectory_path_deviation}} for width-based trajectory analysis
#'
#' @rdname trajectory_convergence
#' @keywords internal
#' @export
trajectory_convergence <- function(x, ...) {
  UseMethod("trajectory_convergence")
}


#' @rdname trajectory_convergence
#' @export
trajectory_convergence.default <- function(x,
                                           dims = c("PC1", "PC2"),
                                           reference_label = NULL,
                                           labels = NULL,
                                           metrics = c("rms", "frechet", "dtw", "correlation"),
                                           interpolate_n = NULL,
                                           trim_fraction = 0.1,
                                           min_coverage = 0.5,
                                           time_digits = 6,
                                           similarity_baseline = c("reference", "zero"),
                                           similarity_scale_multiplier = 1,
                                           stats = TRUE,
                                           verbose = TRUE,
                                           ...) {
  # Input validation
  if (!is.data.frame(x)) stop("Input must be a data frame")

  required_cols <- c("label", "rendition", ".time")
  missing_cols <- setdiff(c(required_cols, dims), names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Validate metrics
  valid_metrics <- c("rms", "frechet", "dtw", "correlation")
  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop(sprintf(
      "Invalid metrics: %s. Choose from: %s",
      paste(invalid_metrics, collapse = ", "),
      paste(valid_metrics, collapse = ", ")
    ))
  }
  similarity_baseline <- match.arg(similarity_baseline)
  if (!is.numeric(similarity_scale_multiplier) ||
      length(similarity_scale_multiplier) != 1 ||
      is.na(similarity_scale_multiplier) ||
      similarity_scale_multiplier <= 0) {
    stop("similarity_scale_multiplier must be a positive number")
  }

  ensure_pkgs("ggplot2", "dplyr")
  if ("dtw" %in% metrics) ensure_pkgs("dtw")

  # Filter labels if specified
  if (!is.null(labels)) {
    available <- unique(x$label)
    missing_labels <- setdiff(labels, available)
    if (length(missing_labels) > 0) {
      stop(sprintf(
        "Labels not found: %s\nAvailable: %s",
        paste(missing_labels, collapse = ", "),
        paste(available, collapse = ", ")
      ))
    }
    x <- x[x$label %in% labels, ]
  }

  x <- bin_trajectory_time_data(x, dims, time_digits)

  all_labels <- unique(x$label)
  if (length(all_labels) == 0) stop("No labels available after filtering")

  # Resolve reference label
  if (is.null(reference_label)) {
    sorted_labels <- sort_labels(all_labels)
    reference_label <- sorted_labels[length(sorted_labels)]
  } else {
    if (!reference_label %in% all_labels) {
      stop(sprintf("Reference label '%s' not found in data", reference_label))
    }
  }

  # Auto-disable stats when there are many groups
  if (stats && length(all_labels) > 6) {
    message(sprintf(
      "Note: %d labels detected. Setting stats = FALSE (statistical tests not meaningful for this many groups).",
      length(all_labels)
    ))
    stats <- FALSE
  }

  if (verbose) {
    message("\n=== Trajectory Convergence Analysis ===")
    message(sprintf("Dimensions       : %s", paste(dims, collapse = ", ")))
    message(sprintf("Reference label  : %s", reference_label))
    message(sprintf("Trim fraction    : %.0f%% each tail", trim_fraction * 100))
    message(sprintf(
      "Min coverage     : %.0f%% of reference renditions per time step",
      min_coverage * 100
    ))
    message(sprintf("Time binning     : %d decimal places", time_digits))
    message(sprintf("Similarity base  : %s", similarity_baseline))
    message(sprintf("Similarity scale : %.3g", similarity_scale_multiplier))
    message(sprintf("Metrics          : %s", paste(metrics, collapse = ", ")))
    if (!is.null(interpolate_n)) {
      message(sprintf("Interpolation    : %d points", interpolate_n))
    } else {
      message("Interpolation    : matched time points")
    }
    message(sprintf("Labels           : %s\n", paste(all_labels, collapse = ", ")))
  }

  # ---- Build reference trajectory ----
  ref_data <- x[x$label == reference_label, ]
  ref_renditions <- unique(ref_data$rendition)
  n_ref_rend <- length(ref_renditions)
  ref_all_times <- sort(unique(ref_data$.time))

  ref_coverage <- vapply(ref_all_times, function(t) {
    sum(ref_data$.time == t) / n_ref_rend
  }, numeric(1))
  ref_times <- ref_all_times[ref_coverage >= min_coverage]

  if (length(ref_times) < 2) {
    stop(sprintf(
      "Reference label '%s' has fewer than 2 binned time steps with %.0f%% coverage",
      reference_label,
      min_coverage * 100
    ))
  }

  reference_path_list <- lapply(ref_times, function(t) {
    t_vals <- ref_data[ref_data$.time == t, dims, drop = FALSE]
    means <- vapply(dims, function(d) mean(t_vals[[d]], trim = trim_fraction, na.rm = TRUE), numeric(1))
    data.frame(.time = t, as.list(means), check.names = FALSE)
  })

  reference_path <- as.data.frame(do.call(rbind, reference_path_list))
  colnames(reference_path) <- c(".time", dims)

  if (verbose) {
    message(sprintf(
      "Reference trajectory: %d time steps (%d raw binned time steps)",
      nrow(reference_path),
      length(ref_all_times)
    ))
  }

  # ---- Interpolate paths if needed ----
  interpolate_path <- function(path_data, dims, n_interp) {
    if (is.null(n_interp)) {
      return(as.matrix(path_data[, dims, drop = FALSE]))
    }

    # Normalize time to [0, 1]
    times <- path_data$.time
    if (length(times) < 2) return(NULL)

    min_t <- min(times)
    max_t <- max(times)
    if (min_t == max_t) return(NULL)

    norm_times <- (times - min_t) / (max_t - min_t)
    grid <- seq(0, 1, length.out = n_interp)

    # Interpolate each dimension
    coords <- as.matrix(path_data[, dims, drop = FALSE])
    interp_coords <- matrix(NA, nrow = n_interp, ncol = length(dims))
    colnames(interp_coords) <- dims

    for (i in seq_along(dims)) {
      interp_coords[, i] <- stats::approx(norm_times, coords[, i],
        xout = grid, rule = 2
      )$y
    }

    interp_coords
  }

  # ---- Helper functions for distance metrics ----
  frechet_dist <- function(P, Q) {
    n <- nrow(P)
    m <- nrow(Q)
    ca <- matrix(NA_real_, n, m)

    d <- function(i, j) sqrt(sum((P[i, ] - Q[j, ])^2))

    ca[1, 1] <- d(1, 1)
    for (i in 2:n) ca[i, 1] <- max(ca[i - 1, 1], d(i, 1))
    for (j in 2:m) ca[1, j] <- max(ca[1, j - 1], d(1, j))
    for (i in 2:n) {
      for (j in 2:m) {
        ca[i, j] <- max(min(ca[i - 1, j], ca[i - 1, j - 1], ca[i, j - 1]), d(i, j))
      }
    }

    ca[n, m]
  }

  rms_distance <- function(trial, ref) {
    if (is.null(trial) || is.null(ref) || nrow(trial) < 1 || nrow(ref) < 1) {
      return(list(rms = NA_real_, mean_dist = NA_real_, max_dist = NA_real_))
    }

    if (!is.null(interpolate_n)) {
      # Already interpolated, same length
      pointwise <- sqrt(rowSums((trial - ref)^2))
    } else {
      # Inner join on time
      trial_times <- trial[, ".time"]
      ref_times <- ref[, ".time"]
      common_times <- intersect(trial_times, ref_times)

      if (length(common_times) == 0) {
        return(list(rms = NA_real_, mean_dist = NA_real_, max_dist = NA_real_))
      }

      trial_aligned <- trial[trial[, ".time"] %in% common_times, dims, drop = FALSE]
      ref_aligned <- ref[ref[, ".time"] %in% common_times, dims, drop = FALSE]

      if (nrow(trial_aligned) < 1) {
        return(list(rms = NA_real_, mean_dist = NA_real_, max_dist = NA_real_))
      }

      pointwise <- sqrt(rowSums((as.matrix(trial_aligned) - as.matrix(ref_aligned))^2))
    }

    list(
      rms = sqrt(mean(pointwise^2, na.rm = TRUE)),
      mean_dist = mean(pointwise, na.rm = TRUE),
      max_dist = max(pointwise, na.rm = TRUE)
    )
  }

  align_to_reference <- function(rend_data) {
    if (!is.null(interpolate_n)) {
      rend_matrix <- interpolate_path(rend_data, dims, interpolate_n)
      ref_matrix <- interpolate_path(reference_path, dims, interpolate_n)
      if (is.null(rend_matrix) || is.null(ref_matrix)) {
        return(NULL)
      }
      return(list(rend = rend_matrix, ref = ref_matrix, n_time = nrow(rend_matrix)))
    }

    common_times <- sort(intersect(rend_data$.time, reference_path$.time))
    if (length(common_times) == 0) {
      return(NULL)
    }

    rend_aligned <- as.matrix(rend_data[rend_data$.time %in% common_times, dims, drop = FALSE])
    ref_aligned <- as.matrix(reference_path[reference_path$.time %in% common_times, dims, drop = FALSE])
    if (nrow(rend_aligned) < 1 || nrow(ref_aligned) < 1) {
      return(NULL)
    }

    list(rend = rend_aligned, ref = ref_aligned, n_time = nrow(rend_aligned))
  }

  metric_distances <- function(rend_data) {
    aligned <- align_to_reference(rend_data)
    empty_result <- list(
      n_time = 0L,
      rms_dist = NA_real_,
      rms_mean_dist = NA_real_,
      rms_max_dist = NA_real_,
      frechet_dist = NA_real_,
      dtw_dist = NA_real_,
      correlation = NA_real_
    )
    if (is.null(aligned)) {
      return(empty_result)
    }

    pointwise <- sqrt(rowSums((aligned$rend - aligned$ref)^2))
    result <- empty_result
    result$n_time <- aligned$n_time
    result$rms_dist <- sqrt(mean(pointwise^2, na.rm = TRUE))
    result$rms_mean_dist <- mean(pointwise, na.rm = TRUE)
    result$rms_max_dist <- max(pointwise, na.rm = TRUE)

    if ("frechet" %in% metrics) {
      result$frechet_dist <- frechet_dist(aligned$rend, aligned$ref)
    }

    if ("dtw" %in% metrics) {
      result$dtw_dist <- tryCatch(
        dtw::dtw(aligned$rend, aligned$ref)$normalizedDistance,
        error = function(e) NA_real_
      )
    }

    if ("correlation" %in% metrics) {
      cors <- vapply(seq_along(dims), function(i) {
        stats::cor(aligned$rend[, i], aligned$ref[, i], use = "complete.obs")
      }, numeric(1))
      result$correlation <- mean(cors, na.rm = TRUE)
    }

    result
  }

  transform_distance <- function(distance, scale) {
    if (!is.finite(distance) || !is.finite(scale) || scale <= 0) {
      return(list(normalized = NA_real_, excess_normalized = NA_real_, similarity = NA_real_))
    }
    normalized <- distance / (scale * similarity_scale_multiplier)
    excess_normalized <- if (similarity_baseline == "reference") {
      max(normalized - 1, 0)
    } else {
      normalized
    }
    list(
      normalized = normalized,
      excess_normalized = excess_normalized,
      similarity = exp(-excess_normalized)
    )
  }

  # ---- Compute reference scale ----
  ref_scales <- vapply(ref_renditions, function(r) {
    rend_data <- ref_data[ref_data$rendition == r, ]
    metric_distances(rend_data)$rms_dist
  }, numeric(1))

  reference_metric_results <- do.call(rbind, lapply(ref_renditions, function(r) {
    rend_data <- ref_data[ref_data$rendition == r, ]
    distances <- metric_distances(rend_data)
    data.frame(
      rendition = r,
      rms = distances$rms_dist,
      frechet = distances$frechet_dist,
      dtw = distances$dtw_dist,
      stringsAsFactors = FALSE
    )
  }))

  reference_scales <- c(
    rms = stats::median(reference_metric_results$rms, na.rm = TRUE),
    frechet = stats::median(reference_metric_results$frechet, na.rm = TRUE),
    dtw = stats::median(reference_metric_results$dtw, na.rm = TRUE)
  )

  reference_scale <- reference_scales[["rms"]]
  if (!is.finite(reference_scale) || reference_scale == 0) {
    reference_scale <- stats::median(ref_scales[ref_scales > 0], na.rm = TRUE)
  }
  if (!is.finite(reference_scale) || reference_scale == 0) {
    reference_scale <- 1.0
    if (verbose) warning("Reference scale could not be computed; using default scale = 1.0")
  }
  reference_scales[["rms"]] <- reference_scale

  for (metric in names(reference_scales)) {
    if (!is.finite(reference_scales[[metric]]) || reference_scales[[metric]] <= 0) {
      reference_scales[[metric]] <- reference_scale
    }
  }

  if (verbose) {
    message("Reference scales  :")
    message(sprintf("  RMS     : %.4f", reference_scales[["rms"]]))
    message(sprintf("  Fr\u00e9chet : %.4f", reference_scales[["frechet"]]))
    message(sprintf("  DTW     : %.4f\n", reference_scales[["dtw"]]))
  }

  # ---- Compute convergence metrics for all trials ----
  convergence_results <- do.call(rbind, lapply(all_labels, function(lbl) {
    lbl_data <- x[x$label == lbl, ]
    renditions <- unique(lbl_data$rendition)

    do.call(rbind, lapply(renditions, function(r) {
      rend_data <- lbl_data[lbl_data$rendition == r, ]
      rend_data <- rend_data[order(rend_data$.time), ]
      distances <- metric_distances(rend_data)

      row_result <- data.frame(
        label = lbl,
        rendition = r,
        n_matched_time = distances$n_time,
        stringsAsFactors = FALSE
      )

      # RMS distance
      if ("rms" %in% metrics) {
        transformed <- transform_distance(distances$rms_dist, reference_scales[["rms"]])
        row_result$rms_dist <- distances$rms_dist
        row_result$rms_mean_dist <- distances$rms_mean_dist
        row_result$rms_max_dist <- distances$rms_max_dist
        row_result$rms_normalized <- transformed$normalized
        row_result$rms_excess_normalized <- transformed$excess_normalized
        row_result$rms_similarity <- transformed$similarity
      }

      # Fr\u00e9chet distance
      if ("frechet" %in% metrics) {
        transformed <- transform_distance(distances$frechet_dist, reference_scales[["frechet"]])
        row_result$frechet_dist <- distances$frechet_dist
        row_result$frechet_normalized <- transformed$normalized
        row_result$frechet_excess_normalized <- transformed$excess_normalized
        row_result$frechet_similarity <- transformed$similarity
      }

      # DTW distance
      if ("dtw" %in% metrics) {
        transformed <- transform_distance(distances$dtw_dist, reference_scales[["dtw"]])
        row_result$dtw_dist <- distances$dtw_dist
        row_result$dtw_normalized <- transformed$normalized
        row_result$dtw_excess_normalized <- transformed$excess_normalized
        row_result$dtw_similarity <- transformed$similarity
      }

      # Correlation
      if ("correlation" %in% metrics) {
        row_result$correlation <- distances$correlation
      }

      row_result
    }))
  }))

  # ---- Summary table ----
  # Simplified: Keep only similarity metrics (mean + sd) for clean output
  # Per-rendition data retains all details for power users
  summary_df <- convergence_results |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(c(
        dplyr::ends_with("_similarity"),
        dplyr::any_of("correlation")
      ),
        list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) |>
    as.data.frame()

  # ---- Statistical tests (exclude reference label) ----
  tests <- NULL
  n_non_ref <- nrow(convergence_results[convergence_results$label != reference_label, ])

  if (stats && length(all_labels) > 1 && n_non_ref > 0) {
    test_data <- convergence_results[convergence_results$label != reference_label, ]

    tests <- list(kruskal = list(), posthoc = list())

    for (metric in c("rms", "frechet", "dtw", "correlation")) {
      metric_col <- if (metric == "correlation") "correlation" else paste0(metric, "_dist")

      if (metric_col %in% names(test_data) && !all(is.na(test_data[[metric_col]]))) {
        kw_test <- stats::kruskal.test(test_data[[metric_col]] ~ test_data$label)
        tests$kruskal[[metric]] <- kw_test

        pw_test <- stats::pairwise.wilcox.test(
          test_data[[metric_col]],
          test_data$label,
          p.adjust.method = "bonferroni",
          exact = FALSE
        )
        tests$posthoc[[metric]] <- pw_test
      }
    }

    if (verbose) {
      message("--- Summary (Reference label excluded from tests) ---")
      print(summary_df)
      message("\nKruskal-Wallis tests:")
      for (metric in names(tests$kruskal)) {
        message(sprintf(
          "  %s: chi-sq = %.2f, p = %.2e",
          metric, tests$kruskal[[metric]]$statistic, tests$kruskal[[metric]]$p.value
        ))
      }
    }
  } else if (verbose) {
    message("--- Summary ---")
    print(summary_df)
  }

  invisible(list(
    type = "convergence",
    dims = dims,
    reference_label = reference_label,
    trim_fraction = trim_fraction,
    min_coverage = min_coverage,
    time_digits = time_digits,
    reference_scale = reference_scale,
    reference_scales = reference_scales,
    similarity_baseline = similarity_baseline,
    similarity_scale_multiplier = similarity_scale_multiplier,
    reference_path = reference_path,
    metrics = metrics,
    interpolate_n = interpolate_n,
    similarity = convergence_results,
    summary = summary_df,
    tests = tests
  ))
}


#' @rdname trajectory_convergence
#' @export
trajectory_convergence.Sap <- function(x,
                                       segment_type = c("motifs", "syllables", "bouts", "segments"),
                                       dims = c("PC1", "PC2"),
                                       reference_label = NULL,
                                       labels = NULL,
                                       metrics = c("rms", "frechet", "dtw", "correlation"),
                                       interpolate_n = NULL,
                                       trim_fraction = 0.1,
                                       min_coverage = 0.5,
                                       time_digits = 6,
                                       similarity_baseline = c("reference", "zero"),
                                       similarity_scale_multiplier = 1,
                                       stats = TRUE,
                                       verbose = TRUE,
                                       ...) {
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  traj_embeds <- x$features[[feature_type]][["traj.embeds"]]
  if (is.null(traj_embeds)) {
    stop(sprintf(
      "Trajectory embeddings not found for %s. Run create_trajectory_matrix() first.",
      segment_type
    ))
  }

  missing_dims <- setdiff(dims, names(traj_embeds))
  if (length(missing_dims) > 0) {
    stop(sprintf("Dimensions not found: %s", paste(missing_dims, collapse = ", ")))
  }

  result <- trajectory_convergence.default(
    x = traj_embeds,
    dims = dims,
    reference_label = reference_label,
    labels = labels,
    metrics = metrics,
    interpolate_n = interpolate_n,
    trim_fraction = trim_fraction,
    min_coverage = min_coverage,
    time_digits = time_digits,
    similarity_baseline = similarity_baseline,
    similarity_scale_multiplier = similarity_scale_multiplier,
    stats = stats,
    verbose = verbose,
    ...
  )

  x$features[[feature_type]][["trajectory_convergence"]] <- result

  invisible(x)
}
#' Plot Trajectory Convergence Results
#'
#' @description
#' Creates multi-panel visualization of trajectory convergence metrics showing
#' how individual trial trajectories compare to a reference path across labels.
#'
#' @param x A list returned by \code{trajectory_convergence()} or a SAP object
#'   with trajectory convergence results
#' @param palette Character. Color palette name for ggplot2 (default: "Set1")
#' @param max_annotations Integer. Maximum number of pairwise significance
#'   annotations to display per panel (default: 10)
#' @param similarity_baseline Similarity transform to use for plotted distance
#'   metrics. \code{"result"} uses the baseline stored in \code{x};
#'   \code{"reference"} treats distances within reference-label variability as
#'   converged; \code{"zero"} uses \code{exp(-normalized_distance)}.
#' @param similarity_scale_multiplier Optional multiplier for metric-specific
#'   reference scales when plotting distance similarities. Defaults to the
#'   multiplier stored in \code{x}, or \code{1} when unavailable.
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables',
#'   'bouts', 'segments')
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' Creates one panel for each computed metric. Panels are ordered as Pointwise
#' Convergence, Shape Correlation, Path-Shape Convergence, and Timing-Adjusted
#' Convergence when those metrics are present. For small label sets (\u22646 labels),
#' uses violin + box plots. For larger label sets, uses mean \u00b1 SE trend lines.
#' Distance similarities are recalculated from raw distances using
#' metric-specific reference scales when available.
#'
#' Reference label is included in results and plots for context but excluded
#' from statistical test annotations (reference defines target, not convergence).
#'
#' @return
#' A patchwork object combining all metric panels (returned invisibly)
#'
#' @examples
#' \dontrun{
#' # Plot trajectory convergence
#' result <- trajectory_convergence(sap)
#' plot_trajectory_convergence(result)
#'
#' # From SAP object directly
#' plot_trajectory_convergence(sap,
#'   segment_type = "motifs"
#' )
#' }
#'
#' @seealso \code{\link{trajectory_convergence}} for analysis computation
#'
#' @rdname plot_trajectory_convergence
#' @export
plot_trajectory_convergence <- function(x, ...) {
  UseMethod("plot_trajectory_convergence")
}


#' @rdname plot_trajectory_convergence
#' @export
plot_trajectory_convergence.default <- function(x,
                                                palette = "Set1",
                                                max_annotations = 10,
                                                similarity_baseline = c("result", "reference", "zero"),
                                                similarity_scale_multiplier = NULL,
                                                ...) {
  result <- x

  # Validate input
  if (!is.list(result) || is.null(result$type) || result$type != "convergence") {
    stop("'x' must be a list returned by trajectory_convergence().")
  }

  ensure_pkgs("ggplot2", "patchwork")

  sim <- result$similarity
  tst <- result$tests
  dims <- result$dims
  metrics <- result$metrics
  ref_label <- result$reference_label
  similarity_baseline <- match.arg(similarity_baseline)
  plot_baseline <- if (similarity_baseline == "result") {
    result$similarity_baseline %||% "reference"
  } else {
    similarity_baseline
  }
  plot_scale_multiplier <- if (is.null(similarity_scale_multiplier)) {
    result$similarity_scale_multiplier %||% 1
  } else {
    similarity_scale_multiplier
  }
  if (!is.numeric(plot_scale_multiplier) ||
      length(plot_scale_multiplier) != 1 ||
      is.na(plot_scale_multiplier) ||
      plot_scale_multiplier <= 0) {
    stop("similarity_scale_multiplier must be NULL or a positive number")
  }

  metric_scale <- function(metric) {
    if (!is.null(result$reference_scales) && metric %in% names(result$reference_scales)) {
      scale <- result$reference_scales[[metric]]
      if (is.finite(scale) && scale > 0) {
        return(scale)
      }
    }

    dist_col <- paste0(metric, "_dist")
    if (dist_col %in% names(sim) && !is.null(ref_label)) {
      ref_dist <- sim[[dist_col]][as.character(sim$label) == as.character(ref_label)]
      scale <- stats::median(ref_dist, na.rm = TRUE)
      if (is.finite(scale) && scale > 0) {
        return(scale)
      }
    }

    NA_real_
  }

  similarity_col <- function(metric) {
    dist_col <- paste0(metric, "_dist")
    norm_col <- paste0(metric, "_normalized")
    plot_col <- paste0(metric, "_plot_similarity")
    scale <- metric_scale(metric)

    if (dist_col %in% names(sim) && is.finite(scale) && scale > 0) {
      normalized <- sim[[dist_col]] / (scale * plot_scale_multiplier)
      sim[[plot_col]] <<- if (plot_baseline == "reference") {
        exp(-pmax(normalized - 1, 0))
      } else {
        exp(-normalized)
      }
      plot_col
    } else if (norm_col %in% names(sim)) {
      sim[[plot_col]] <<- if (plot_baseline == "reference") {
        exp(-pmax(sim[[norm_col]] - 1, 0))
      } else {
        exp(-sim[[norm_col]])
      }
      plot_col
    } else {
      stored_col <- paste0(metric, "_similarity")
      if (stored_col %in% names(sim)) stored_col else plot_col
    }
  }

  labs_order <- sort_labels(unique(as.character(sim$label)))
  many_labs <- length(labs_order) > 6
  pal_map <- make_pal(labs_order, palette)

  panels <- list()

  # Helper to create individual metric panels
  plot_metric <- function(metric_name, y_col, title, y_label = title) {
    if (!metric_name %in% metrics || !y_col %in% names(sim)) {
      return(NULL)
    }

    if (all(is.na(sim[[y_col]]))) {
      return(NULL)
    }

    kw <- if (!is.null(tst) && metric_name %in% names(tst$kruskal)) {
      fmt_p(tst$kruskal[[metric_name]]$p.value)
    } else {
      NULL
    }

    if (many_labs) {
      trend_panel(sim, "label", y_col, title, kw, labs_order, y_label = y_label)
    } else {
      p <- panel(sim, "label", y_col, title, kw, pal_map, labs_order, y_label = y_label)

      if (!is.null(tst) && metric_name %in% names(tst$posthoc)) {
        p <- add_brackets(p,
          brackets(sim[[y_col]], tst$posthoc[[metric_name]], labs_order, max_annotations)
        )
      }

      p
    }
  }

  # ---- Panel: RMS Similarity ----
  rms_col <- similarity_col("rms")
  if (!is.null(plot_metric("rms", rms_col, "Pointwise Convergence"))) {
    panels$p_rms <- plot_metric(
      "rms",
      rms_col,
      "Pointwise Convergence",
      "RMS similarity"
    )
  }

  # ---- Panel: Correlation ----
  if (!is.null(plot_metric("correlation", "correlation", "Shape Correlation"))) {
    panels$p_cor <- plot_metric(
      "correlation",
      "correlation",
      "Shape Correlation",
      "Correlation similarity"
    )
  }

  # ---- Panel: Fr\u00e9chet Similarity ----
  frechet_col <- similarity_col("frechet")
  if (!is.null(plot_metric("frechet", frechet_col, "Path-Shape Convergence"))) {
    panels$p_frechet <- plot_metric(
      "frechet",
      frechet_col,
      "Path-Shape Convergence",
      "Fr\u00e9chet similarity"
    )
  }

  # ---- Panel: DTW Similarity ----
  dtw_col <- similarity_col("dtw")
  if (!is.null(plot_metric("dtw", dtw_col, "Timing-Adjusted Convergence"))) {
    panels$p_dtw <- plot_metric(
      "dtw",
      dtw_col,
      "Timing-Adjusted Convergence",
      "DTW similarity"
    )
  }

  if (length(panels) == 0) {
    stop("No valid metrics to plot.")
  }

  # Combine panels
  combined <- Reduce("+", panels) +
    patchwork::plot_annotation(
      title = sprintf("Trajectory Convergence to Reference: %s", ref_label),
      subtitle = paste(
        "Dimensions:", paste(dims, collapse = " + "),
        "| Similarity baseline:", plot_baseline
      ),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      )
    )

  withCallingHandlers(
    print(combined),
    warning = function(w) {
      if (grepl("annotation$theme is not a valid theme", w$message, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  invisible(combined)
}


#' @rdname plot_trajectory_convergence
#' @export
plot_trajectory_convergence.Sap <- function(x,
                                            segment_type = c("motifs", "syllables", "bouts", "segments"),
                                            palette = "Set1",
                                            max_annotations = 10,
                                            similarity_baseline = c("result", "reference", "zero"),
                                            similarity_scale_multiplier = NULL,
                                            ...) {
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  result <- x$features[[feature_type]][["trajectory_convergence"]]
  if (is.null(result)) {
    stop(sprintf(
      "No trajectory convergence results found for %s. Run trajectory_convergence() first.",
      segment_type
    ))
  }

  plot_trajectory_convergence.default(
    x = result,
    palette = palette,
    max_annotations = max_annotations,
    similarity_baseline = similarity_baseline,
    similarity_scale_multiplier = similarity_scale_multiplier,
    ...
  )
}
