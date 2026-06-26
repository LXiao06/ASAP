# Trajectory Variability Analysis
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
#' @keywords internal
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

  # Attach .source_row from original data
  if (".source_row" %in% names(x)) {
    source_lookup <- unique(x[, c("label", "rendition", ".source_row")])
    dispersion_results <- merge(dispersion_results, source_lookup,
      by = c("label", "rendition"), all.x = TRUE)
  } else {
    dispersion_results$.source_row <- NA_integer_
  }

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
        label = lbl, rendition = r,
        .source_row = if (".source_row" %in% names(x) && nrow(rend_data) > 0) unique(rend_data$.source_row)[1] else NA_integer_,
        path_length = path_len,
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
      rend_data <- rend_data[order(rend_data$.time), ]
      keep_cols <- c("label", "rendition", ".time", dims)
      if (".source_row" %in% names(rend_data)) keep_cols <- c(keep_cols, ".source_row")
      rend_data <- rend_data[, keep_cols]
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
        .source_row = if (".source_row" %in% names(x) && nrow(rend_data) > 0) unique(rend_data$.source_row)[1] else NA_integer_,
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
        .source_row = if (".source_row" %in% names(x) && nrow(rend_data) > 0) unique(rend_data$.source_row)[1] else NA_integer_,
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
#'   \item{\code{"dispersion"}}{3 panels: Mean Pairwise Distance \eqn{\cdot} Centroid
#'     Dispersion \eqn{\cdot} Path Length}
#'   \item{\code{"path_deviation"}}{3 panels: Total RMS \eqn{\cdot} Orthogonal RMS \eqn{\cdot}
#'     Parallel RMS}
#'   \item{\code{"umap_occupancy"}}{4 panels: Occupied Fraction \eqn{\cdot} Occupancy
#'     Entropy \eqn{\cdot} Peripheral Fraction \eqn{\cdot} kNN Dispersion}
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
