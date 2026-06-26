# Trajectory Similarity Analysis
# Update date : Jun. 23, 2026

#' Trajectory Similarity to Reference Path
#'
#' @description
#' Measures trajectory similarity by quantifying how
#' individual trial PC trajectories compare to a reference path (e.g.,
#' adult/crystallized song). Uses complementary distance metrics (RMS,
#' Fr\'echet, DTW) with unified variability scaling and correlation analysis.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or SAP object
#' @param dims Character vector of dimension columns to use (default: c("PC1", "PC2"))
#' @param reference_label Character. Label to use as reference trajectory
#'   (default: NULL, uses last sorted label)
#' @param labels Optional character vector of labels to include
#' @param metrics Character vector of metrics to compute, or \code{"all"} for all
#'   four metrics (default: c("rms", "correlation"))
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
#' **Similarity Metrics:**
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
#' similarity target).
#'
#' @return
#' For default method: A list (returned invisibly) with:
#' \itemize{
#'   \item \code{type}: Character string \code{"similarity"}
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
#'   \code{x$features[[feature_type]][["trajectory_similarity"]]}
#'
#' @examples
#' \dontrun{
#' # Compute similarity from trajectory embeddings
#' result <- trajectory_similarity(sap$features$motif$traj.embeds,
#'   dims = c("PC1", "PC2")
#' )
#'
#' # From SAP object
#' sap <- trajectory_similarity(sap)
#'
#' # Custom reference and metrics
#' result <- trajectory_similarity(sap,
#'   reference_label = "Adult",
#'   metrics = c("rms", "correlation")
#' )
#'
#' # All four metrics (RMS, correlation, Frechet, DTW)
#' result_all <- trajectory_similarity(sap,
#'   metrics = "all"
#' )
#'
#' # Access results
#' result$summary # per-label statistics
#' result$tests # statistical tests
#' }
#'
#' @seealso \code{\link{plot_trajectory_similarity}} for visualization,
#'   \code{\link{trajectory_path_deviation}} for width-based trajectory analysis
#'
#' @rdname trajectory_similarity
#' @keywords internal
#' @export
trajectory_similarity <- function(x, ...) {
  UseMethod("trajectory_similarity")
}


#' @rdname trajectory_similarity
#' @export
trajectory_similarity.default <- function(x,
                                           dims = c("PC1", "PC2"),
                                           reference_label = NULL,
                                           labels = NULL,
                                           metrics = c("rms", "correlation"),
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
  if (identical(metrics, "all")) {
    metrics <- valid_metrics
  }
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
    message("\n=== Trajectory Similarity Analysis ===")
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
  # ---- Helper functions for distance metrics ----
  # TODO: extract to standalone helper
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

  # TODO: extract to standalone helper
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

  # TODO: extract to standalone helper
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

  # ---- Compute similarity metrics for all trials ----
  similarity_results <- do.call(rbind, lapply(all_labels, function(lbl) {
    lbl_data <- x[x$label == lbl, ]
    renditions <- unique(lbl_data$rendition)

    do.call(rbind, lapply(renditions, function(r) {
      rend_data <- lbl_data[lbl_data$rendition == r, ]
      rend_data <- rend_data[order(rend_data$.time), ]
      distances <- metric_distances(rend_data)

      row_result <- data.frame(
        label = lbl,
        rendition = r,
        .source_row = if (".source_row" %in% names(x) && nrow(rend_data) > 0) unique(rend_data$.source_row)[1] else NA_integer_,
        n_matched_time = distances$n_time,
        stringsAsFactors = FALSE
      )

      # RMS distance
      if ("rms" %in% metrics) {
        transformed <- transform_distance(distances$rms_dist, reference_scales[["rms"]],
          similarity_scale_multiplier, similarity_baseline)
        row_result$rms_dist <- distances$rms_dist
        row_result$rms_mean_dist <- distances$rms_mean_dist
        row_result$rms_max_dist <- distances$rms_max_dist
        row_result$rms_normalized <- transformed$normalized
        row_result$rms_excess_normalized <- transformed$excess_normalized
        row_result$rms_similarity <- transformed$similarity
      }

      # Fr\u00e9chet distance
      if ("frechet" %in% metrics) {
        transformed <- transform_distance(distances$frechet_dist, reference_scales[["frechet"]],
          similarity_scale_multiplier, similarity_baseline)
        row_result$frechet_dist <- distances$frechet_dist
        row_result$frechet_normalized <- transformed$normalized
        row_result$frechet_excess_normalized <- transformed$excess_normalized
        row_result$frechet_similarity <- transformed$similarity
      }

      # DTW distance
      if ("dtw" %in% metrics) {
        transformed <- transform_distance(distances$dtw_dist, reference_scales[["dtw"]],
          similarity_scale_multiplier, similarity_baseline)
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
  summary_df <- similarity_results |>
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
  n_non_ref <- nrow(similarity_results[similarity_results$label != reference_label, ])

  if (stats && length(all_labels) > 1 && n_non_ref > 0) {
    test_data <- similarity_results[similarity_results$label != reference_label, ]

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
    type = "similarity",
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
    similarity = similarity_results,
    summary = summary_df,
    tests = tests
  ))
}


#' @rdname trajectory_similarity
#' @export
trajectory_similarity.Sap <- function(x,
                                       segment_type = c("motifs", "syllables", "bouts", "segments"),
                                       dims = c("PC1", "PC2"),
                                       reference_label = NULL,
                                       labels = NULL,
                                       metrics = c("rms", "correlation"),
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

  result <- trajectory_similarity.default(
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

  x$features[[feature_type]][["trajectory_similarity"]] <- result

  invisible(x)
}
#' Interpolate Trajectory Path
#'
#' @description
#' Internal helper used by \code{trajectory_similarity.default()}.
#' Interpolates a trajectory path to a fixed number of equally spaced points.
#'
#' @param path_data Data frame with .time and dimension columns
#' @param dims Character vector of dimension column names
#' @param n_interp Integer number of interpolation points, or NULL to return matrix as-is
#'
#' @keywords internal
#' @noRd
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


#' Fr\u00e9chet Distance Between Trajectories
#'
#' @description
#' Internal helper used by \code{trajectory_similarity.default()}.
#' Computes the discrete Fr\u00e9chet distance between two trajectories.
#'
#' @param P First trajectory matrix (n x d)
#' @param Q Second trajectory matrix (m x d)
#'
#' @keywords internal
#' @noRd
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


#' Transform Distance to Similarity Score
#'
#' @description
#' Internal helper used by \code{trajectory_similarity.default()}.
#' Converts a distance value to a similarity score using metric-specific
#' scaling and baseline.
#'
#' @param distance Numeric distance value
#' @param scale Numeric reference scale
#' @param similarity_scale_multiplier Numeric multiplier applied to scale
#' @param similarity_baseline Character "reference" or "zero"
#'
#' @keywords internal
#' @noRd
transform_distance <- function(distance, scale,
                                similarity_scale_multiplier,
                                similarity_baseline) {
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


#' Plot Trajectory Similarity Results
#'
#' @description
#' Creates multi-panel visualization of trajectory similarity metrics showing
#' how individual trial trajectories compare to a reference path across labels.
#'
#' @param x A list returned by \code{trajectory_similarity()} or a SAP object
#'   with trajectory similarity results
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
#' Similarity, Shape Correlation, Path-Shape Similarity, and Timing-Adjusted
#' Similarity when those metrics are present. For small label sets (\eqn{\le}6 labels),
#' uses violin + box plots. For larger label sets, uses mean \eqn{\pm} SE trend lines.
#' Distance similarities are recalculated from raw distances using
#' metric-specific reference scales when available.
#'
#' Reference label is included in results and plots for context but excluded
#' from statistical test annotations (reference defines target, not similarity).
#'
#' @return
#' A patchwork object combining all metric panels (returned invisibly)
#'
#' @examples
#' \dontrun{
#' # Plot trajectory similarity
#' result <- trajectory_similarity(sap)
#' plot_trajectory_similarity(result)
#'
#' # From SAP object directly
#' plot_trajectory_similarity(sap,
#'   segment_type = "motifs"
#' )
#'
#' # With all four metrics
#' result_all <- trajectory_similarity(sap, metrics = "all")
#' plot_trajectory_similarity(result_all)
#' }
#'
#' @seealso \code{\link{trajectory_similarity}} for analysis computation
#'
#' @rdname plot_trajectory_similarity
#' @export
plot_trajectory_similarity <- function(x, ...) {
  UseMethod("plot_trajectory_similarity")
}


#' @rdname plot_trajectory_similarity
#' @export
plot_trajectory_similarity.default <- function(x,
                                                palette = "Set1",
                                                max_annotations = 10,
                                                similarity_baseline = c("result", "reference", "zero"),
                                                similarity_scale_multiplier = NULL,
                                                ...) {
  result <- x

  # Validate input
  if (!is.list(result) || is.null(result$type) || result$type != "similarity") {
    stop("'x' must be a list returned by trajectory_similarity().")
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

  # TODO: extract to standalone helper
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

  # TODO: extract to standalone helper
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
  # TODO: extract to standalone helper
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
  if (!is.null(plot_metric("rms", rms_col, "Pointwise Similarity"))) {
    panels$p_rms <- plot_metric(
      "rms",
      rms_col,
      "Pointwise Similarity",
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
  if (!is.null(plot_metric("frechet", frechet_col, "Path-Shape Similarity"))) {
    panels$p_frechet <- plot_metric(
      "frechet",
      frechet_col,
      "Path-Shape Similarity",
      "Fr\u00e9chet similarity"
    )
  }

  # ---- Panel: DTW Similarity ----
  dtw_col <- similarity_col("dtw")
  if (!is.null(plot_metric("dtw", dtw_col, "Timing-Adjusted Similarity"))) {
    panels$p_dtw <- plot_metric(
      "dtw",
      dtw_col,
      "Timing-Adjusted Similarity",
      "DTW similarity"
    )
  }

  if (length(panels) == 0) {
    stop("No valid metrics to plot.")
  }

  # Combine panels
  combined <- Reduce("+", panels) +
    patchwork::plot_annotation(
      title = sprintf("Trajectory Similarity to Reference: %s", ref_label),
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


#' @rdname plot_trajectory_similarity
#' @export
plot_trajectory_similarity.Sap <- function(x,
                                            segment_type = c("motifs", "syllables", "bouts", "segments"),
                                            palette = "Set1",
                                            max_annotations = 10,
                                            similarity_baseline = c("result", "reference", "zero"),
                                            similarity_scale_multiplier = NULL,
                                            ...) {
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  result <- x$features[[feature_type]][["trajectory_similarity"]]
  if (is.null(result)) {
    stop(sprintf(
      "No trajectory similarity results found for %s. Run trajectory_similarity() first.",
      segment_type
    ))
  }

  plot_trajectory_similarity.default(
    x = result,
    palette = palette,
    max_annotations = max_annotations,
    similarity_baseline = similarity_baseline,
    similarity_scale_multiplier = similarity_scale_multiplier,
    ...
  )
}
