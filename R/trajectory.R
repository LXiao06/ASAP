# Song Trajectory Analysis
# Update date : Feb. 7, 2026

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

  # Combine results
  spc_list <- spc_list[!sapply(spc_list, is.null)]
  if (length(spc_list) == 0) {
    stop("No valid spectrograms were generated")
  }
  spectrogram_matrix <- do.call(cbind, spc_list) |> t()

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
  window_start_time <- start + (0:(max_windows - 1)) * step_size

  # Create sliding window data frame
  sliding_windows <- data.frame(
    filename = x$filename[i],
    day_post_hatch = x$day_post_hatch[i],
    label = x$label[i],
    rendition = i,
    selec = seq_along(window_start_time),
    start_time = window_start_time,
    end_time = window_start_time + window_size,
    .time = window_start_time - start
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
#' # Use filtered data directly with trajectory_variability
#' result <- trajectory_variability(traj_clean, dims = c("PC1", "PC2"))
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


# Trajectory Variability
# Update date : May 1, 2026

#' Trajectory Variability Analysis
#'
#' @description
#' Quantifies trajectory variability across experimental conditions using three
#' complementary metrics computed in PCA or UMAP embedding space.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or SAP object
#' @param dims Character vector of dimension columns to use (default: c("PC1", "PC2"))
#' @param labels Optional character vector of labels to include
#' @param max_pairs Maximum number of pairwise comparisons per label (default: 5000)
#' @param seed Random seed for reproducible pair sampling (default: 222)
#' @param palette Color palette name for plotting (default: "Set1")
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
#' @return A list (returned invisibly) with the following elements:
#' \itemize{
#'   \item \code{pairwise}: Data frame of pairwise distances (label, pair_id, mean_dist)
#'   \item \code{dispersion}: Data frame of centroid dispersion (label, rendition, dispersion)
#'   \item \code{path_length}: Data frame of path lengths (label, rendition, path_length)
#'   \item \code{summary}: Summary table with mean and SD for each metric per label
#'   \item \code{tests}: List of statistical test results
#' }
#'
#' A summary plot is printed as a side effect.
#'
#' @examples
#' \dontrun{
#' # From SAP object using PC dimensions
#' result <- trajectory_variability(sap)
#'
#' # Using UMAP dimensions
#' result <- trajectory_variability(sap, dims = c("UMAP1", "UMAP2"))
#'
#' # From trajectory embeddings data frame directly
#' result <- trajectory_variability(sap$features$motif$traj.embeds,
#'   dims = c("PC1", "PC2")
#' )
#'
#' # Access results
#' result$summary # summary table
#' result$tests # statistical tests
#' result$dispersion # per-rendition dispersion values
#' }
#'
#' @rdname trajectory_variability
#' @export
trajectory_variability <- function(x, ...) {
  UseMethod("trajectory_variability")
}


#' @rdname trajectory_variability
#' @export
trajectory_variability.default <- function(x,
                                           dims = c("PC1", "PC2"),
                                           labels = NULL,
                                           max_pairs = 5000,
                                           seed = 222,
                                           palette = "Set1",
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

  ensure_pkgs("ggplot2", "dplyr", "patchwork")

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

  if (verbose) {
    message("\n=== Trajectory Variability Analysis ===")
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

    # Compute mean distance for each pair
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
    message(sprintf("\nKruskal-Wallis tests:"))
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

  # ==== Plot ====
  format_p <- function(p) {
    if (p < 0.001) {
      sprintf("p = %.1e", p)
    } else {
      sprintf("p = %.3f", p)
    }
  }

  p1 <- ggplot2::ggplot(
    pairwise_results,
    ggplot2::aes(
      x = .data$label, y = .data$mean_dist,
      fill = .data$label
    )
  ) +
    ggplot2::geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::labs(
      title = "Mean Pairwise Distance",
      subtitle = format_p(test_pw$p.value),
      y = "Distance", x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  p2 <- ggplot2::ggplot(
    dispersion_results,
    ggplot2::aes(
      x = .data$label, y = .data$dispersion,
      fill = .data$label
    )
  ) +
    ggplot2::geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::labs(
      title = "Centroid Dispersion",
      subtitle = format_p(test_disp$p.value),
      y = "Mean Distance to Centroid", x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  p3 <- ggplot2::ggplot(
    path_length_results,
    ggplot2::aes(
      x = .data$label, y = .data$path_length,
      fill = .data$label
    )
  ) +
    ggplot2::geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::labs(
      title = "Trajectory Path Length",
      subtitle = format_p(test_pl$p.value),
      y = "Path Length", x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )

  combined_plot <- (p1 + p2 + p3) +
    patchwork::plot_annotation(
      title = "Trajectory Variability Comparison",
      subtitle = paste("Dimensions:", paste(dims, collapse = " + "))
    )

  print(combined_plot)

  # Return results
  invisible(list(
    pairwise = pairwise_results,
    dispersion = dispersion_results,
    path_length = path_length_results,
    summary = summary_table,
    tests = tests
  ))
}


#' @rdname trajectory_variability
#' @export
trajectory_variability.Sap <- function(x,
                                       segment_type = c(
                                         "motifs", "syllables",
                                         "bouts", "segments"
                                       ),
                                       dims = c("PC1", "PC2"),
                                       labels = NULL,
                                       max_pairs = 5000,
                                       seed = 222,
                                       palette = "Set1",
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

  result <- trajectory_variability.default(
    x = traj_embeds,
    dims = dims,
    labels = labels,
    max_pairs = max_pairs,
    seed = seed,
    palette = palette,
    verbose = verbose,
    ...
  )

  invisible(result)
}


# Trajectory Width Variability
# Update date : May 1, 2026

#' Trajectory Width Variability Analysis
#'
#' @description
#' Quantifies rendition-to-rendition trajectory "width" by measuring residual
#' spread around each label's mean trajectory and decomposing that spread into
#' orthogonal and parallel components relative to the local trajectory tangent.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or SAP object
#' @param dims Character vector of dimension columns to use (default: c("PC1", "PC2"))
#' @param trim_fraction Trim fraction for robust mean trajectory estimation
#'   (default: 0.1)
#' @param min_coverage Minimum fraction of renditions that must cover a time step
#'   for it to contribute to the reference trajectory (default: 0.5)
#' @param labels Optional character vector of labels to include
#' @param palette Color palette name for plotting (default: "Set1")
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables',
#'   'bouts', 'segments')
#' @param ... Additional arguments
#'
#' @details
#' For each label, the function builds a robust mean trajectory in the requested
#' dimensions, estimates a local tangent vector at each time step, and decomposes
#' each rendition's residual into:
#' \describe{
#'   \item{Total RMS Residual}{Overall deviation from the label-specific mean trajectory}
#'   \item{Orthogonal RMS Residual}{Deviation perpendicular to the local tangent;
#'     a direct measure of trajectory width / jitter around the backbone}
#'   \item{Parallel RMS Residual}{Deviation along the local tangent; often reflects
#'     phase or advance-lag variability}
#' }
#'
#' @return A list (returned invisibly) with the following elements:
#' \itemize{
#'   \item \code{width}: Per-rendition width metrics
#'   \item \code{summary}: Summary table with mean and SD for each metric per label
#'   \item \code{mean_trajectories}: Label-specific mean trajectories
#'   \item \code{tangent_vectors}: Label-specific unit tangent vectors
#'   \item \code{tests}: Kruskal-Wallis and pairwise Wilcoxon tests when multiple
#'     labels are present
#' }
#'
#' @examples
#' \dontrun{
#' result <- trajectory_width_variability(sap, dims = c("PC1", "PC2"))
#' result$summary
#' result$width
#' }
#'
#' @export
trajectory_width_variability <- function(x, ...) {
  UseMethod("trajectory_width_variability")
}


#' @rdname trajectory_width_variability
#' @export
trajectory_width_variability.default <- function(x,
                                                 dims = c("PC1", "PC2"),
                                                 trim_fraction = 0.1,
                                                 min_coverage = 0.5,
                                                 labels = NULL,
                                                 palette = "Set1",
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

  ensure_pkgs("ggplot2", "dplyr", "patchwork")

  if (!is.null(labels)) {
    missing_labels <- setdiff(labels, unique(x$label))
    if (length(missing_labels) > 0) {
      stop(sprintf("Labels not found: %s", paste(missing_labels, collapse = ", ")))
    }
    x <- x[x$label %in% labels, ]
  }

  all_labels <- unique(x$label)
  if (length(all_labels) == 0) stop("No labels available after filtering")

  if (verbose) {
    message("\n=== Trajectory Width Variability Analysis ===")
    message(sprintf("Dimensions    : %s", paste(dims, collapse = ", ")))
    message(sprintf("Trim fraction : %.0f%% each tail", trim_fraction * 100))
    message(sprintf(
      "Min coverage  : %.0f%% of renditions per time step",
      min_coverage * 100
    ))
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
  format_p <- function(p) {
    if (is.null(p) || is.na(p)) {
      return("p = NA")
    }
    if (p < 0.001) sprintf("p = %.1e", p) else sprintf("p = %.3f", p)
  }

  get_pairwise_p <- function(pmat, g1, g2) {
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

  build_pairwise_annotations <- function(values, posthoc_obj, labels_in_order) {
    comps <- utils::combn(labels_in_order, 2, simplify = FALSE)
    y_range <- range(values, na.rm = TRUE)
    span <- diff(y_range)
    if (!is.finite(span) || span == 0) span <- max(abs(y_range), na.rm = TRUE)
    if (!is.finite(span) || span == 0) span <- 1
    base_y <- y_range[2] + 0.08 * span
    step_y <- 0.12 * span

    ann <- do.call(rbind, lapply(seq_along(comps), function(i) {
      comp <- comps[[i]]
      p_val <- get_pairwise_p(posthoc_obj$p.value, comp[1], comp[2])
      data.frame(
        x1 = match(comp[1], labels_in_order),
        x2 = match(comp[2], labels_in_order),
        y = base_y + (i - 1) * step_y,
        label = format_p(p_val),
        stringsAsFactors = FALSE
      )
    }))

    ann$y_text <- ann$y + 0.025 * span
    ann$y_tip <- ann$y - 0.02 * span
    ann$y_max <- max(ann$y_text) + 0.06 * span
    ann
  }

  if (length(unique(width_results$label)) > 1) {
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
      kruskal = list(
        total = test_total,
        orthogonal = test_orth,
        parallel = test_parallel
      ),
      posthoc = list(
        total = posthoc_total,
        orthogonal = posthoc_orth,
        parallel = posthoc_parallel
      )
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

  p_total <- ggplot2::ggplot(
    width_results,
    ggplot2::aes(x = .data$label, y = .data$total_rms, fill = .data$label)
  ) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::labs(
      title = "Total RMS Residual",
      subtitle = if (is.null(tests)) NULL else format_p(tests$kruskal$total$p.value),
      y = "Distance to Mean Trajectory",
      x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  p_orth <- ggplot2::ggplot(
    width_results,
    ggplot2::aes(x = .data$label, y = .data$orthogonal_rms, fill = .data$label)
  ) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::labs(
      title = "Orthogonal RMS Residual",
      subtitle = if (is.null(tests)) NULL else format_p(tests$kruskal$orthogonal$p.value),
      y = "Width Around Trajectory Backbone",
      x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  p_parallel <- ggplot2::ggplot(
    width_results,
    ggplot2::aes(x = .data$label, y = .data$parallel_rms, fill = .data$label)
  ) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::labs(
      title = "Parallel RMS Residual",
      subtitle = if (is.null(tests)) NULL else format_p(tests$kruskal$parallel$p.value),
      y = "Along-Trajectory Variability",
      x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )

  combined_plot <- (p_total + p_orth + p_parallel) +
    patchwork::plot_annotation(
      title = "Trajectory Width Variability Comparison",
      subtitle = paste("Dimensions:", paste(dims, collapse = " + "))
    )

  print(combined_plot)

  invisible(list(
    width = width_results,
    summary = summary_df,
    mean_trajectories = mean_trajectories,
    tangent_vectors = tangent_vectors,
    tests = tests
  ))
}


#' @rdname trajectory_width_variability
#' @export
trajectory_width_variability.Sap <- function(x,
                                             segment_type = c(
                                               "motifs", "syllables",
                                               "bouts", "segments"
                                             ),
                                             dims = c("PC1", "PC2"),
                                             trim_fraction = 0.1,
                                             min_coverage = 0.5,
                                             labels = NULL,
                                             palette = "Set1",
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

  result <- trajectory_width_variability.default(
    x = traj_embeds,
    dims = dims,
    trim_fraction = trim_fraction,
    min_coverage = min_coverage,
    labels = labels,
    palette = palette,
    verbose = verbose,
    ...
  )

  invisible(result)
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
#' @param palette Color palette name for plotting (default: "Set1")
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
#'   \item \code{occupancy}: Per-rendition occupancy metrics
#'   \item \code{summary}: Summary table with mean and SD for each metric per label
#'   \item \code{annotated_points}: Original data with occupancy annotations
#'   \item \code{bin_counts}: Shared-grid counts per label and bin
#'   \item \code{grid_info}: Grid settings and peripheral threshold metadata
#'   \item \code{tests}: Kruskal-Wallis and pairwise Wilcoxon tests
#' }
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
                                              palette = "Set1",
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

  ensure_pkgs("ggplot2", "dplyr", "patchwork", "RANN")

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

  format_p <- function(p) {
    if (is.null(p) || is.na(p)) {
      return("p = NA")
    }
    if (p < 0.001) sprintf("p = %.1e", p) else sprintf("p = %.3f", p)
  }

  get_pairwise_p <- function(pmat, g1, g2) {
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

  build_pairwise_annotations <- function(values, posthoc_obj, labels_in_order) {
    comps <- utils::combn(labels_in_order, 2, simplify = FALSE)
    y_range <- range(values, na.rm = TRUE)
    span <- diff(y_range)
    if (!is.finite(span) || span == 0) span <- max(abs(y_range), na.rm = TRUE)
    if (!is.finite(span) || span == 0) span <- 1
    base_y <- y_range[2] + 0.08 * span
    step_y <- 0.12 * span

    ann <- do.call(rbind, lapply(seq_along(comps), function(i) {
      comp <- comps[[i]]
      p_val <- get_pairwise_p(posthoc_obj$p.value, comp[1], comp[2])
      data.frame(
        x1 = match(comp[1], labels_in_order),
        x2 = match(comp[2], labels_in_order),
        y = base_y + (i - 1) * step_y,
        label = format_p(p_val),
        stringsAsFactors = FALSE
      )
    }))

    ann$y_text <- ann$y + 0.025 * span
    ann$y_tip <- ann$y - 0.02 * span
    ann$y_max <- max(ann$y_text) + 0.06 * span
    ann
  }

  tests <- NULL
  if (length(unique(occupancy_results$label)) > 1) {
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

  label_order <- unique(as.character(occupancy_results$label))
  pair_ann_occ <- if (is.null(tests)) {
    NULL
  } else {
    build_pairwise_annotations(
      occupancy_results$occupied_fraction,
      tests$posthoc$occupied_fraction,
      label_order
    )
  }
  pair_ann_ent <- if (is.null(tests)) {
    NULL
  } else {
    build_pairwise_annotations(
      occupancy_results$occupancy_entropy_norm,
      tests$posthoc$entropy,
      label_order
    )
  }
  pair_ann_per <- if (is.null(tests)) {
    NULL
  } else {
    build_pairwise_annotations(
      occupancy_results$peripheral_fraction,
      tests$posthoc$peripheral_fraction,
      label_order
    )
  }
  pair_ann_knn <- if (is.null(tests)) {
    NULL
  } else {
    build_pairwise_annotations(
      occupancy_results$knn_dispersion,
      tests$posthoc$knn_dispersion,
      label_order
    )
  }

  add_pairwise_layers <- function(plot_obj, ann_df) {
    if (is.null(ann_df) || nrow(ann_df) == 0) {
      return(plot_obj)
    }
    plot_obj +
      ggplot2::geom_segment(
        data = ann_df,
        ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, yend = .data$y),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_segment(
        data = ann_df,
        ggplot2::aes(x = .data$x1, xend = .data$x1, y = .data$y_tip, yend = .data$y),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_segment(
        data = ann_df,
        ggplot2::aes(x = .data$x2, xend = .data$x2, y = .data$y_tip, yend = .data$y),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = ann_df,
        ggplot2::aes(x = (.data$x1 + .data$x2) / 2, y = .data$y_text, label = .data$label),
        inherit.aes = FALSE,
        size = 3
      ) +
      ggplot2::coord_cartesian(ylim = c(NA, ann_df$y_max[1]))
  }

  p_occ <- ggplot2::ggplot(
    occupancy_results,
    ggplot2::aes(x = .data$label, y = .data$occupied_fraction, fill = .data$label)
  ) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.35, size = 0.9) +
    ggplot2::labs(
      title = "Occupied Fraction",
      subtitle = if (is.null(tests)) NULL else format_p(tests$kruskal$occupied_fraction$p.value),
      y = "Visited Grid Fraction",
      x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  p_occ <- add_pairwise_layers(p_occ, pair_ann_occ)

  p_ent <- ggplot2::ggplot(
    occupancy_results,
    ggplot2::aes(x = .data$label, y = .data$occupancy_entropy_norm, fill = .data$label)
  ) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.35, size = 0.9) +
    ggplot2::labs(
      title = "Occupancy Entropy",
      subtitle = if (is.null(tests)) NULL else format_p(tests$kruskal$entropy$p.value),
      y = "Normalized Entropy",
      x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  p_ent <- add_pairwise_layers(p_ent, pair_ann_ent)

  p_per <- ggplot2::ggplot(
    occupancy_results,
    ggplot2::aes(x = .data$label, y = .data$peripheral_fraction, fill = .data$label)
  ) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.35, size = 0.9) +
    ggplot2::labs(
      title = "Peripheral Fraction",
      subtitle = if (is.null(tests)) NULL else format_p(tests$kruskal$peripheral_fraction$p.value),
      y = "Fraction in Sparse UMAP Bins",
      x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  p_per <- add_pairwise_layers(p_per, pair_ann_per)

  p_knn <- ggplot2::ggplot(
    occupancy_results,
    ggplot2::aes(x = .data$label, y = .data$knn_dispersion, fill = .data$label)
  ) +
    ggplot2::geom_violin(alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.35, size = 0.9) +
    ggplot2::labs(
      title = "Same-Label kNN Dispersion",
      subtitle = if (is.null(tests)) NULL else format_p(tests$kruskal$knn_dispersion$p.value),
      y = "Mean kNN Distance",
      x = NULL
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  p_knn <- add_pairwise_layers(p_knn, pair_ann_knn)

  combined_plot <- (p_occ + p_ent + p_per + p_knn) +
    patchwork::plot_annotation(
      title = "Trajectory UMAP Occupancy Comparison",
      subtitle = paste("Dimensions:", paste(dims, collapse = " + "))
    )

  print(combined_plot)

  invisible(list(
    occupancy = occupancy_results,
    summary = summary_df,
    annotated_points = x,
    bin_counts = bin_counts,
    grid_info = list(
      dims = dims,
      grid_n = grid_n,
      total_bins = total_bins,
      x_breaks = x_breaks,
      y_breaks = y_breaks,
      bin_area = bin_area,
      peripheral_quantile = peripheral_quantile,
      peripheral_cut = peripheral_cut,
      peripheral_bins = peripheral_bins
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
                                          palette = "Set1",
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
    palette = palette,
    verbose = verbose,
    ...
  )

  invisible(result)
}
