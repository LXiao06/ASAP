# Trajectory Discrimination Analysis
# Update date : Sep. 21, 2026

# Suppress R CMD check notes about internal variables and ggplot aesthetics
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "across", "all_of", ".data", "rendition", ".time", "label",
    "variability", "variability_se", "ci_lower", "ci_upper",
    "d_prime", "centroid_dist", "bhattacharyya", "silhouette",
    "pair", "label_1", "label_2",
    "p_fdr", "is_significant", "p_value", "mean_variability",
    "auc_variability", "mean_d_prime", "peak_d_prime",
    "mean_bhattacharyya", "peak_bhattacharyya",
    "mean_silhouette", "peak_silhouette",
    "p_value_d_prime", "p_value_dist",
    "p_value_bhattacharyya", "p_value_silhouette",
    "metric_val", "ymin", "ymax",
    "xmin", "xmax", "y_sig"
  ))
}

# -------------------------------------------------------------------------
# Trajectory Timecourse Analysis
# -------------------------------------------------------------------------

#' Trajectory Timecourse Analysis
#'
#' @description
#' Quantifies continuous, time-resolved trajectory variability along the
#' temporal progression of song renditions for each experimental condition or
#' social context (e.g. female-directed vs. undirected song). Supports
#' unbalanced sample sizes, multiple groups (>2), and pointwise statistical
#' comparisons across time steps with FDR correction.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or a SAP object.
#' @param dims Character vector of dimension columns to use (default: \code{c("PC1", "PC2")}).
#' @param labels Optional character vector of labels to include. If \code{NULL} (default),
#'   all unique labels in \code{x$label} are used.
#' @param metric Character. Method used to quantify dispersion at each time step:
#'   \itemize{
#'     \item \code{"centroid_dist"} (default): Mean Euclidean distance of each rendition
#'       to the label's centroid coordinate at time \eqn{t}. Provides per-rendition deviations.
#'     \item \code{"trace"}: Total variance at time \eqn{t}, computed as the sum of
#'       variances across selected dimensions.
#'     \item \code{"mad"}: Median absolute deviation from the median coordinate at time \eqn{t}.
#'   }
#' @param min_coverage Numeric in \code{(0, 1]}. Minimum fraction of renditions within
#'   a label that must cover a time step for it to be analyzed (default: \code{0.5}).
#' @param time_digits Integer. Number of decimal places used to bin \code{.time}
#'   to ensure exact temporal alignment across renditions (default: \code{6}).
#' @param smooth_span Optional numeric. LOESS smoothing span to fit a smoothed
#'   variability curve across time (default: \code{NULL}).
#' @param stats Logical. If \code{TRUE} (default) and \eqn{\ge 2} labels are present,
#'   runs pointwise non-parametric tests (Wilcoxon rank-sum for 2 groups; Kruskal-Wallis
#'   and post-hoc pairwise Wilcoxon for >2 groups) on per-rendition deviations across time,
#'   with Benjamini-Hochberg FDR correction.
#' @param verbose Logical. Whether to print progress messages (default: \code{TRUE}).
#' @param segment_type For SAP objects: Type of segments to analyze
#'   (\code{'motifs'}, \code{'syllables'}, \code{'bouts'}, \code{'segments'}).
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' Unlike scalar variability measures that collapse an entire rendition into a single
#' summary number, \code{trajectory_timecourse()} tracks how acoustic exploration and
#' motor precision fluctuate dynamically across the syllable or motif timeline.
#'
#' For each label and each aligned time step \eqn{t}:
#' \describe{
#'   \item{Centroid Distance}{
#'     Computes the label's centroid \eqn{\bar{\mathbf{X}}_k(t)} and calculates
#'     the Euclidean distance for each rendition \eqn{d_{k,r}(t) = \|\mathbf{X}_{k,r}(t) - \bar{\mathbf{X}}_k(t)\|}.
#'     Mean variability is \eqn{V_k(t) = \frac{1}{N_k(t)}\sum d_{k,r}(t)}, with standard error
#'     \eqn{SE_k(t) = \text{sd}(d_{k,r}(t)) / \sqrt{N_k(t)}} and 95\% confidence intervals.
#'   }
#'   \item{Trace / Total Variance}{
#'     Computes the sum of sample variances across dimensions: \eqn{\sum_d \text{Var}(\mathbf{X}_{k,d}(t))}.
#'   }
#'   \item{Area Under Curve (AUC)}{
#'     Total trajectory variability integrated over time using the trapezoidal rule,
#'     providing an overall measure of motor exploration across the rendition.
#'   }
#' }
#'
#' @return
#' For the default method, a list of class \code{c("trajectory_timecourse", "list")} containing:
#' \itemize{
#'   \item \code{timecourse}: Data frame with columns \code{label}, \code{.time}, \code{n_renditions},
#'     \code{variability}, \code{variability_se}, \code{ci_lower}, and \code{ci_upper}
#'     (and \code{variability_smoothed} if \code{smooth_span} is supplied).
#'   \item \code{rendition_deviations}: Data frame of per-rendition deviations at each time step
#'     (\code{NULL} if \code{metric = "trace"}).
#'   \item \code{summary}: Summary table per label with \code{n_renditions}, \code{mean_variability},
#'     \code{sd_variability}, \code{peak_time}, \code{peak_variability}, and \code{auc_variability}.
#'   \item \code{tests}: Data frame of pointwise comparisons across labels, including raw
#'     and Benjamini-Hochberg FDR-adjusted p-values (\code{NULL} if \code{stats = FALSE} or only 1 label).
#'   \item \code{metric}: Character string specifying the metric used.
#'   \item \code{dims}: Dimensions analyzed.
#'   \item \code{type}: Character string \code{"timecourse"}, used by \code{\link{plot_trajectory_discrimination}()}.
#' }
#' For SAP objects, the updated SAP object with results stored in
#' \code{x$features[[feature_type]][["trajectory_timecourse"]]} (returned invisibly).
#'
#' @examples
#' \dontrun{
#' # From trajectory embeddings data frame
#' tc <- trajectory_timecourse(sap$features$motif$traj.embeds, dims = c("PC1", "PC2"))
#' tc$summary
#'
#' # Plot results
#' plot_trajectory_discrimination(tc)
#'
#' # From SAP object directly
#' sap <- trajectory_timecourse(sap, segment_type = "motifs")
#' }
#'
#' @export
trajectory_timecourse <- function(x, ...) {
  UseMethod("trajectory_timecourse")
}


#' @rdname trajectory_timecourse
#' @export
trajectory_timecourse.default <- function(x,
                                          dims = c("PC1", "PC2"),
                                          labels = NULL,
                                          metric = c("centroid_dist", "trace", "mad"),
                                          min_coverage = 0.5,
                                          time_digits = 6,
                                          smooth_span = NULL,
                                          stats = TRUE,
                                          verbose = TRUE,
                                          ...) {
  metric <- match.arg(metric)
  ensure_pkgs("dplyr")

  # ---- Input Validation ----
  if (!is.data.frame(x)) {
    stop("'x' must be a data frame containing trajectory embeddings.")
  }
  required_cols <- c("label", "rendition", ".time")
  missing_req <- setdiff(required_cols, names(x))
  if (length(missing_req) > 0) {
    stop(sprintf(
      "Missing required columns in 'x': %s",
      paste(missing_req, collapse = ", ")
    ))
  }
  missing_dims <- setdiff(dims, names(x))
  if (length(missing_dims) > 0) {
    stop(sprintf(
      "Dimensions not found in 'x': %s",
      paste(missing_dims, collapse = ", ")
    ))
  }
  if (!is.numeric(min_coverage) || min_coverage <= 0 || min_coverage > 1) {
    stop("'min_coverage' must be a single numeric value in (0, 1].")
  }

  # ---- Filter labels ----
  if (!is.null(labels)) {
    x <- x[x$label %in% labels, , drop = FALSE]
    if (nrow(x) == 0) {
      stop("No data remaining after filtering by specified 'labels'.")
    }
  }

  # Ensure label is character
  x$label <- as.character(x$label)
  present_labels <- sort_labels(unique(x$label))
  n_labels <- length(present_labels)

  if (n_labels == 0) {
    stop("No labels found in 'x'.")
  }

  if (verbose) {
    message(sprintf(
      "Calculating trajectory timecourse for %d label(s) across %s...",
      n_labels, paste(dims, collapse = ", ")
    ))
  }

  # ---- Bin time data ----
  binned <- bin_trajectory_time_data(x, dims = dims, time_digits = time_digits)

  # ---- Compute total renditions per label for coverage filtering ----
  total_n_per_label <- tapply(binned$rendition, binned$label, function(r) length(unique(r)))

  # ---- Calculate pointwise timecourse per label ----
  timecourse_list <- list()
  rendition_dev_list <- list()

  for (lbl in present_labels) {
    df_lbl <- binned[binned$label == lbl, , drop = FALSE]
    n_total <- total_n_per_label[[lbl]]
    min_renditions <- ceiling(n_total * min_coverage)

    # Group by time step
    time_groups <- split(df_lbl, df_lbl$.time)

    records <- list()
    lbl_rend_devs <- list()

    for (t_str in names(time_groups)) {
      sub_t <- time_groups[[t_str]]
      n_cur <- nrow(sub_t)
      if (n_cur < min_renditions || n_cur < 2) next

      t_val <- as.numeric(t_str)
      coords <- as.matrix(sub_t[, dims, drop = FALSE])

      if (metric == "centroid_dist") {
        centroid <- colMeans(coords)
        devs <- sqrt(rowSums(sweep(coords, 2, centroid, "-")^2))
        v_mean <- mean(devs)
        v_sd <- stats::sd(devs)
        v_se <- v_sd / sqrt(n_cur)

        records[[length(records) + 1]] <- data.frame(
          label = lbl,
          .time = t_val,
          n_renditions = n_cur,
          variability = v_mean,
          variability_se = v_se,
          ci_lower = max(0, v_mean - 1.96 * v_se),
          ci_upper = v_mean + 1.96 * v_se,
          stringsAsFactors = FALSE
        )

        lbl_rend_devs[[length(lbl_rend_devs) + 1]] <- data.frame(
          label = lbl,
          rendition = sub_t$rendition,
          .time = t_val,
          deviation = devs,
          stringsAsFactors = FALSE
        )
      } else if (metric == "trace") {
        variances <- apply(coords, 2, stats::var)
        v_tot <- sum(variances)
        # Jackknife SE for total variance
        jk_vars <- numeric(n_cur)
        for (i in seq_len(n_cur)) {
          jk_vars[i] <- sum(apply(coords[-i, , drop = FALSE], 2, stats::var))
        }
        jk_se <- sqrt(((n_cur - 1) / n_cur) * sum((jk_vars - mean(jk_vars))^2))

        records[[length(records) + 1]] <- data.frame(
          label = lbl,
          .time = t_val,
          n_renditions = n_cur,
          variability = v_tot,
          variability_se = jk_se,
          ci_lower = max(0, v_tot - 1.96 * jk_se),
          ci_upper = v_tot + 1.96 * jk_se,
          stringsAsFactors = FALSE
        )
      } else if (metric == "mad") {
        med_coord <- apply(coords, 2, stats::median)
        devs <- sqrt(rowSums(sweep(coords, 2, med_coord, "-")^2))
        v_mad <- stats::median(devs)
        v_se <- (1.4826 * stats::median(abs(devs - v_mad))) / sqrt(n_cur)

        records[[length(records) + 1]] <- data.frame(
          label = lbl,
          .time = t_val,
          n_renditions = n_cur,
          variability = v_mad,
          variability_se = v_se,
          ci_lower = max(0, v_mad - 1.96 * v_se),
          ci_upper = v_mad + 1.96 * v_se,
          stringsAsFactors = FALSE
        )

        lbl_rend_devs[[length(lbl_rend_devs) + 1]] <- data.frame(
          label = lbl,
          rendition = sub_t$rendition,
          .time = t_val,
          deviation = devs,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(records) > 0) {
      tc_lbl <- do.call(rbind, records)
      tc_lbl <- tc_lbl[order(tc_lbl$.time), ]

      # Optional LOESS smoothing
      if (!is.null(smooth_span) && nrow(tc_lbl) >= 4) {
        fit <- stats::loess(variability ~ .time, data = tc_lbl, span = smooth_span)
        tc_lbl$variability_smoothed <- stats::predict(fit)
      }

      timecourse_list[[lbl]] <- tc_lbl
    }
    if (length(lbl_rend_devs) > 0) {
      rendition_dev_list[[lbl]] <- do.call(rbind, lbl_rend_devs)
    }
  }

  if (length(timecourse_list) == 0) {
    stop("No time points satisfied min_coverage for any label.")
  }

  timecourse_df <- do.call(rbind, timecourse_list)
  rownames(timecourse_df) <- NULL

  rend_devs_df <- if (length(rendition_dev_list) > 0) {
    out_dev <- do.call(rbind, rendition_dev_list)
    rownames(out_dev) <- NULL
    out_dev
  } else {
    NULL
  }

  # ---- Summary Table per Label ----
  summary_rows <- list()
  for (lbl in unique(timecourse_df$label)) {
    sub_df <- timecourse_df[timecourse_df$label == lbl, ]
    peak_idx <- which.max(sub_df$variability)
    times <- sub_df$.time
    vals <- sub_df$variability

    # Trapezoidal rule for AUC
    auc_val <- if (length(times) >= 2) {
      dt <- diff(times)
      sum(dt * (vals[-1] + vals[-length(vals)]) / 2)
    } else {
      NA_real_
    }

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      label = lbl,
      n_renditions = total_n_per_label[[lbl]],
      n_timepoints = nrow(sub_df),
      mean_variability = mean(vals),
      sd_variability = stats::sd(vals),
      peak_time = times[peak_idx],
      peak_variability = vals[peak_idx],
      auc_variability = auc_val,
      stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, summary_rows)

  # ---- Pointwise Statistical Comparisons across Labels ----
  tests_df <- NULL
  if (stats && n_labels >= 2 && !is.null(rend_devs_df)) {
    # Perform pairwise comparisons for all pairs of labels
    label_pairs <- utils::combn(present_labels, 2, simplify = FALSE)
    pair_tests <- list()

    for (pair in label_pairs) {
      lbl1 <- pair[1]
      lbl2 <- pair[2]
      pair_name <- paste0(lbl1, "_vs_", lbl2)

      d1 <- rend_devs_df[rend_devs_df$label == lbl1, ]
      d2 <- rend_devs_df[rend_devs_df$label == lbl2, ]
      common_times <- intersect(unique(d1$.time), unique(d2$.time))

      if (length(common_times) == 0) next

      time_records <- list()
      for (t_val in sort(common_times)) {
        v1 <- d1$deviation[d1$.time == t_val]
        v2 <- d2$deviation[d2$.time == t_val]

        if (length(v1) >= 2 && length(v2) >= 2) {
          wt <- suppressWarnings(stats::wilcox.test(v1, v2))
          m1 <- mean(v1)
          m2 <- mean(v2)

          time_records[[length(time_records) + 1]] <- data.frame(
            pair = pair_name,
            label_1 = lbl1,
            label_2 = lbl2,
            .time = t_val,
            mean_1 = m1,
            mean_2 = m2,
            diff = m1 - m2,
            ratio = if (m2 > 0) m1 / m2 else NA_real_,
            p_value = wt$p.value,
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(time_records) > 0) {
        pair_df <- do.call(rbind, time_records)
        pair_df$p_fdr <- stats::p.adjust(pair_df$p_value, method = "BH")
        pair_df$is_significant <- pair_df$p_fdr < 0.05
        pair_tests[[length(pair_tests) + 1]] <- pair_df
      }
    }

    if (length(pair_tests) > 0) {
      tests_df <- do.call(rbind, pair_tests)
      rownames(tests_df) <- NULL
    }
  }

  result <- list(
    timecourse = timecourse_df,
    rendition_deviations = rend_devs_df,
    summary = summary_df,
    tests = tests_df,
    metric = metric,
    dims = dims,
    type = "timecourse"
  )

  class(result) <- c("trajectory_timecourse", "list")
  result
}


#' @rdname trajectory_timecourse
#' @export
trajectory_timecourse.Sap <- function(x,
                                      segment_type = c("motifs", "syllables", "bouts", "segments"),
                                      dims = c("PC1", "PC2"),
                                      labels = NULL,
                                      metric = c("centroid_dist", "trace", "mad"),
                                      min_coverage = 0.5,
                                      time_digits = 6,
                                      smooth_span = NULL,
                                      stats = TRUE,
                                      verbose = TRUE,
                                      ...) {
  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  traj_embeds <- x$features[[feature_type]][["traj.embeds"]]
  if (is.null(traj_embeds)) {
    stop(sprintf(
      "No trajectory embeddings found for '%s'. Run create_trajectory_matrix() first.",
      segment_type
    ))
  }

  # Borrow missing dimensions from traj_mat if available
  missing_dims <- setdiff(dims, names(traj_embeds))
  if (length(missing_dims) > 0) {
    traj_mat <- x$features[[feature_type]][["traj_mat"]]
    dims_in_mat <- if (!is.null(traj_mat)) intersect(missing_dims, names(traj_mat)) else character(0)
    if (length(dims_in_mat) > 0 && nrow(traj_mat) == nrow(traj_embeds)) {
      if (verbose) {
        message(sprintf(
          "Dimensions [%s] not in traj.embeds; borrowing from traj_mat (run_pca() output).",
          paste(dims_in_mat, collapse = ", ")
        ))
      }
      traj_embeds[dims_in_mat] <- traj_mat[dims_in_mat]
    }
  }

  result <- trajectory_timecourse.default(
    x = traj_embeds,
    dims = dims,
    labels = labels,
    metric = metric,
    min_coverage = min_coverage,
    time_digits = time_digits,
    smooth_span = smooth_span,
    stats = stats,
    verbose = verbose,
    ...
  )

  x$features[[feature_type]][["trajectory_timecourse"]] <- result
  invisible(x)
}


# -------------------------------------------------------------------------
# Trajectory Separation Analysis
# -------------------------------------------------------------------------

#' Trajectory Separation and Discriminability Analysis
#'
#' @description
#' Quantifies the geometric separation, probabilistic divergence, and statistical
#' discriminability between song trajectories recorded across different conditions,
#' experimental treatments, or social contexts (e.g. directed vs. undirected song)
#' in PCA or UMAP embedding space.
#'
#' @param x An object to analyze: a trajectory embeddings data frame or a SAP object.
#' @param dims Character vector of dimension columns to use (default: \code{c("PC1", "PC2")}).
#' @param labels Optional character vector of labels to include. If \code{NULL} (default),
#'   all unique labels in \code{x$label} are used.
#' @param metrics Character vector of separation metrics to compute (default computes all four):
#'   \itemize{
#'     \item \code{"d_prime"}: Fisher criterion / d-prime between group centroids,
#'       normalized by pooled standard deviation: \eqn{\|\bar{\mathbf{X}}_A(t) - \bar{\mathbf{X}}_B(t)\| / s_{\text{pooled}}(t)}.
#'     \item \code{"centroid_dist"}: Euclidean distance between group centroids at each time step.
#'     \item \code{"bhattacharyya"}: Bhattacharyya distance measuring both mean divergence
#'       and covariance contraction / difference in shape:
#'       \eqn{D_B(t) = D_{\text{mean}}(t) + D_{\text{cov}}(t)}. Uses ridge regularization
#'       to ensure stability with small sample sizes.
#'     \item \code{"silhouette"}: Average silhouette score at each time step, measuring
#'       whether song renditions cluster into distinct behavioral categories.
#'   }
#' @param min_coverage Numeric in \code{(0, 1]}. Minimum fraction of renditions within
#'   each label that must cover a time step for it to be compared (default: \code{0.5}).
#' @param time_digits Integer. Number of decimal places used to bin \code{.time} (default: \code{6}).
#' @param stats Logical. If \code{TRUE} (default), runs rendition-level permutation tests
#'   by randomly shuffling label identities to generate empirical null distributions
#'   and p-values for trajectory separation.
#' @param n_perm Integer. Number of permutations for statistical testing (default: \code{500}).
#' @param seed Integer. Random seed for reproducible permutation testing (default: \code{222}).
#' @param verbose Logical. Whether to print progress messages (default: \code{TRUE}).
#' @param segment_type For SAP objects: Type of segments to analyze
#'   (\code{'motifs'}, \code{'syllables'}, \code{'bouts'}, \code{'segments'}).
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' \code{trajectory_separation()} evaluates whether and when trajectories diverge between
#' conditions in embedding space. For every unique pair of labels \eqn{(A, B)}:
#' \describe{
#'   \item{Pointwise Centroid Distance}{
#'     \eqn{D_{AB}(t) = \|\bar{\mathbf{X}}_A(t) - \bar{\mathbf{X}}_B(t)\|}. Measures the absolute
#'     spatial divergence between the mean acoustic paths at time \eqn{t}.
#'   }
#'   \item{Pointwise d-prime (\eqn{d'})}{
#'     Standardized separation ratio:
#'     \deqn{d'_{AB}(t) = \frac{\|\bar{\mathbf{X}}_A(t) - \bar{\mathbf{X}}_B(t)\|}{s_{\text{pooled}}(t) + \epsilon}}
#'     where \eqn{s_{\text{pooled}}(t)} is the pooled sample standard deviation across renditions in
#'     both groups at time \eqn{t}. High \eqn{d'} indicates that between-group separation far exceeds
#'     within-group rendition variability.
#'   }
#'   \item{Bhattacharyya Distance}{
#'     Measures probabilistic divergence between Gaussian approximations of the two trajectory
#'     clouds:
#'     \deqn{D_B = \frac{1}{8}(\boldsymbol{\mu}_A - \boldsymbol{\mu}_B)^T \boldsymbol{\Sigma}_{\text{reg}}^{-1} (\boldsymbol{\mu}_A - \boldsymbol{\mu}_B) + \frac{1}{2}\ln\left(\frac{\det \boldsymbol{\Sigma}_{\text{reg}}}{\sqrt{\det\boldsymbol{\Sigma}_{A, \text{reg}} \det\boldsymbol{\Sigma}_{B, \text{reg}}}}\right)}
#'     The second term explicitly penalizes variance contraction (differences in cloud spread),
#'     capturing the hallmark variability reduction of courtship song even when mean paths align.
#'   }
#'   \item{Silhouette Score}{
#'     Computes \eqn{s(i, t) = \frac{b(i, t) - a(i, t)}{\max(a(i, t), b(i, t))}} for each rendition
#'     at time \eqn{t}. Values near \eqn{+1} indicate distinct non-overlapping clusters; values near
#'     \eqn{0} indicate completely overlapping distributions.
#'   }
#'   \item{Rendition-Level Permutation Test}{
#'     Shuffles label assignments across entire renditions (preserving within-rendition temporal
#'     autocorrelation) and recalculates the mean separation metrics.
#'     The empirical p-value is the proportion of permutations yielding separation greater than
#'     or equal to the observed value.
#'   }
#' }
#'
#' @return
#' For the default method, a list of class \code{c("trajectory_separation", "list")} containing:
#' \itemize{
#'   \item \code{pointwise}: Data frame of separation metrics at each time step for every label pair:
#'     \code{pair}, \code{label_1}, \code{label_2}, \code{.time}, \code{centroid_dist}, \code{d_prime},
#'     \code{bhattacharyya} (if computed), and \code{silhouette} (if computed).
#'   \item \code{pairwise_summary}: Summary table per pair with mean and peak values for each computed metric,
#'     along with permutation test p-values (\code{p_value_d_prime}, \code{p_value_dist}, etc.).
#'   \item \code{omnibus}: Summary statistics across all pairs if >2 labels are present.
#'   \item \code{metrics}: Character vector of metrics computed.
#'   \item \code{dims}: Dimensions analyzed.
#'   \item \code{type}: Character string \code{"separation"}, used by \code{\link{plot_trajectory_discrimination}()}.
#' }
#' For SAP objects, the updated SAP object with results stored in
#' \code{x$features[[feature_type]][["trajectory_separation"]]} (returned invisibly).
#'
#' @examples
#' \dontrun{
#' # From trajectory embeddings data frame
#' sep <- trajectory_separation(sap$features$motif$traj.embeds, dims = c("PC1", "PC2"))
#' sep$pairwise_summary
#'
#' # Plot results
#' plot_trajectory_discrimination(sep, score_col = "d_prime")
#' plot_trajectory_discrimination(sep, score_col = "bhattacharyya")
#'
#' # From SAP object directly
#' sap <- trajectory_separation(sap, segment_type = "motifs")
#' }
#'
#' @export
trajectory_separation <- function(x, ...) {
  UseMethod("trajectory_separation")
}


#' @rdname trajectory_separation
#' @export
trajectory_separation.default <- function(x,
                                          dims = c("PC1", "PC2"),
                                          labels = NULL,
                                          metrics = c("d_prime", "centroid_dist", "bhattacharyya", "silhouette"),
                                          min_coverage = 0.5,
                                          time_digits = 6,
                                          stats = TRUE,
                                          n_perm = 500,
                                          seed = 222,
                                          verbose = TRUE,
                                          ...) {
  ensure_pkgs("dplyr")

  valid_metrics <- c("d_prime", "centroid_dist", "bhattacharyya", "silhouette")
  matched_metrics <- match.arg(metrics, valid_metrics, several.ok = TRUE)

  # ---- Input Validation ----
  if (!is.data.frame(x)) {
    stop("'x' must be a data frame containing trajectory embeddings.")
  }
  required_cols <- c("label", "rendition", ".time")
  missing_req <- setdiff(required_cols, names(x))
  if (length(missing_req) > 0) {
    stop(sprintf(
      "Missing required columns in 'x': %s",
      paste(missing_req, collapse = ", ")
    ))
  }
  missing_dims <- setdiff(dims, names(x))
  if (length(missing_dims) > 0) {
    stop(sprintf(
      "Dimensions not found in 'x': %s",
      paste(missing_dims, collapse = ", ")
    ))
  }
  if (!is.numeric(min_coverage) || min_coverage <= 0 || min_coverage > 1) {
    stop("'min_coverage' must be a single numeric value in (0, 1].")
  }

  # ---- Filter labels ----
  if (!is.null(labels)) {
    x <- x[x$label %in% labels, , drop = FALSE]
    if (nrow(x) == 0) {
      stop("No data remaining after filtering by specified 'labels'.")
    }
  }

  x$label <- as.character(x$label)
  present_labels <- sort_labels(unique(x$label))
  n_labels <- length(present_labels)

  if (n_labels < 2) {
    stop("At least 2 unique labels are required for trajectory separation analysis.")
  }

  if (verbose) {
    message(sprintf(
      "Calculating trajectory separation across %d labels (%s) in %s...",
      n_labels, paste(present_labels, collapse = ", "), paste(dims, collapse = ", ")
    ))
  }

  # ---- Bin time data ----
  binned <- bin_trajectory_time_data(x, dims = dims, time_digits = time_digits)

  # Helper function to compute separation metrics for a single pair of labels
  compute_pair_metrics <- function(df_pair, lbl1, lbl2) {
    d1 <- df_pair[df_pair$label == lbl1, , drop = FALSE]
    d2 <- df_pair[df_pair$label == lbl2, , drop = FALSE]

    n_tot1 <- length(unique(d1$rendition))
    n_tot2 <- length(unique(d2$rendition))
    min_r1 <- ceiling(n_tot1 * min_coverage)
    min_r2 <- ceiling(n_tot2 * min_coverage)

    t1 <- split(d1, d1$.time)
    t2 <- split(d2, d2$.time)

    common_times <- intersect(names(t1), names(t2))
    if (length(common_times) == 0) return(NULL)

    d_dim <- length(dims)
    records <- list()

    for (t_str in common_times) {
      sub1 <- t1[[t_str]]
      sub2 <- t2[[t_str]]
      n1 <- nrow(sub1)
      n2 <- nrow(sub2)

      if (n1 < min_r1 || n2 < min_r2 || n1 < 2 || n2 < 2) next

      t_val <- as.numeric(t_str)
      c1 <- as.matrix(sub1[, dims, drop = FALSE])
      c2 <- as.matrix(sub2[, dims, drop = FALSE])

      mean1 <- colMeans(c1)
      mean2 <- colMeans(c2)

      row_res <- list(
        pair = paste0(lbl1, "_vs_", lbl2),
        label_1 = lbl1,
        label_2 = lbl2,
        .time = t_val
      )

      # 1. Centroid distance
      c_dist <- sqrt(sum((mean1 - mean2)^2))
      if ("centroid_dist" %in% matched_metrics) {
        row_res$centroid_dist <- c_dist
      }

      # 2. d-prime
      if ("d_prime" %in% matched_metrics) {
        var1 <- sum(apply(c1, 2, stats::var))
        var2 <- sum(apply(c2, 2, stats::var))
        df_pool <- (n1 - 1) + (n2 - 1)
        pooled_var <- if (df_pool > 0) ((n1 - 1) * var1 + (n2 - 1) * var2) / df_pool else (var1 + var2) / 2
        s_pool <- sqrt(max(0, pooled_var))
        dp <- c_dist / (s_pool + 1e-8)
        row_res$d_prime <- dp
      }

      # 3. Bhattacharyya distance
      if ("bhattacharyya" %in% matched_metrics) {
        cov1 <- stats::cov(c1)
        cov2 <- stats::cov(c2)
        cov_pooled <- (cov1 + cov2) / 2

        # Ridge regularization to prevent singularity when N <= D
        reg_scale <- max(mean(diag(cov_pooled)), 1e-6) * 1e-4
        diag_reg <- diag(reg_scale, d_dim)
        c1_reg <- cov1 + diag_reg
        c2_reg <- cov2 + diag_reg
        pool_reg <- cov_pooled + diag_reg

        diff_m <- matrix(mean1 - mean2, ncol = 1)
        inv_pool <- tryCatch(solve(pool_reg), error = function(e) diag(1 / (diag(pool_reg) + 1e-6), d_dim))
        d_mean <- as.numeric(0.125 * t(diff_m) %*% inv_pool %*% diff_m)

        det_p <- max(det(pool_reg), 1e-12)
        det_1 <- max(det(c1_reg), 1e-12)
        det_2 <- max(det(c2_reg), 1e-12)
        d_cov <- 0.5 * log(det_p / sqrt(det_1 * det_2))

        b_dist <- max(0, d_mean) + max(0, d_cov)
        row_res$bhattacharyya <- b_dist
      }

      # 4. Silhouette score
      if ("silhouette" %in% matched_metrics) {
        comb_coords <- rbind(c1, c2)
        comb_dists <- as.matrix(stats::dist(comb_coords))

        idx1 <- seq_len(n1)
        idx2 <- seq(n1 + 1, n1 + n2)

        # Intra-cluster distances (excluding self)
        a1 <- rowSums(comb_dists[idx1, idx1, drop = FALSE]) / (n1 - 1)
        b1 <- rowMeans(comb_dists[idx1, idx2, drop = FALSE])
        s1 <- (b1 - a1) / pmax(a1, b1, 1e-8)

        a2 <- rowSums(comb_dists[idx2, idx2, drop = FALSE]) / (n2 - 1)
        b2 <- rowMeans(comb_dists[idx2, idx1, drop = FALSE])
        s2 <- (b2 - a2) / pmax(a2, b2, 1e-8)

        sil_val <- mean(c(s1, s2), na.rm = TRUE)
        row_res$silhouette <- sil_val
      }

      records[[length(records) + 1]] <- as.data.frame(row_res, stringsAsFactors = FALSE)
    }

    if (length(records) == 0) return(NULL)
    out <- do.call(rbind, records)
    out[order(out$.time), ]
  }

  label_pairs <- utils::combn(present_labels, 2, simplify = FALSE)
  pointwise_list <- list()
  summary_list <- list()

  for (pair in label_pairs) {
    lbl1 <- pair[1]
    lbl2 <- pair[2]
    pair_name <- paste0(lbl1, "_vs_", lbl2)

    df_pair <- binned[binned$label %in% c(lbl1, lbl2), , drop = FALSE]
    pw <- compute_pair_metrics(df_pair, lbl1, lbl2)

    if (is.null(pw) || nrow(pw) == 0) {
      if (verbose) {
        warning(sprintf("No common time points satisfied coverage for pair '%s'.", pair_name))
      }
      next
    }

    pointwise_list[[pair_name]] <- pw

    # Pairwise summary building
    sum_entry <- list(
      pair = pair_name,
      label_1 = lbl1,
      label_2 = lbl2
    )

    if ("centroid_dist" %in% matched_metrics) {
      max_c_idx <- which.max(pw$centroid_dist)
      sum_entry$mean_centroid_dist <- mean(pw$centroid_dist)
      sum_entry$peak_centroid_dist <- pw$centroid_dist[max_c_idx]
      sum_entry$peak_centroid_time <- pw$.time[max_c_idx]
    }
    if ("d_prime" %in% matched_metrics) {
      max_d_idx <- which.max(pw$d_prime)
      sum_entry$mean_d_prime <- mean(pw$d_prime)
      sum_entry$peak_d_prime <- pw$d_prime[max_d_idx]
      sum_entry$peak_d_prime_time <- pw$.time[max_d_idx]
    }
    if ("bhattacharyya" %in% matched_metrics) {
      max_b_idx <- which.max(pw$bhattacharyya)
      sum_entry$mean_bhattacharyya <- mean(pw$bhattacharyya)
      sum_entry$peak_bhattacharyya <- pw$bhattacharyya[max_b_idx]
      sum_entry$peak_bhattacharyya_time <- pw$.time[max_b_idx]
    }
    if ("silhouette" %in% matched_metrics) {
      max_s_idx <- which.max(pw$silhouette)
      sum_entry$mean_silhouette <- mean(pw$silhouette)
      sum_entry$peak_silhouette <- pw$silhouette[max_s_idx]
      sum_entry$peak_silhouette_time <- pw$.time[max_s_idx]
    }

    # Permutation test
    if (stats && n_perm > 0) {
      if (!is.null(seed)) set.seed(seed)

      rend_table <- unique(df_pair[, c("label", "rendition")])
      rend_labels <- rend_table$label
      n_rends <- nrow(rend_table)

      perm_dp <- numeric(n_perm)
      perm_cd <- numeric(n_perm)
      perm_bh <- numeric(n_perm)
      perm_sl <- numeric(n_perm)

      for (p in seq_len(n_perm)) {
        shuffled_labels <- sample(rend_labels)
        label_map <- stats::setNames(shuffled_labels, rend_table$rendition)

        df_perm <- df_pair
        df_perm$label <- label_map[as.character(df_perm$rendition)]

        pw_perm <- compute_pair_metrics(df_perm, lbl1, lbl2)
        if (!is.null(pw_perm) && nrow(pw_perm) > 0) {
          if ("d_prime" %in% matched_metrics) perm_dp[p] <- mean(pw_perm$d_prime)
          if ("centroid_dist" %in% matched_metrics) perm_cd[p] <- mean(pw_perm$centroid_dist)
          if ("bhattacharyya" %in% matched_metrics) perm_bh[p] <- mean(pw_perm$bhattacharyya)
          if ("silhouette" %in% matched_metrics) perm_sl[p] <- mean(pw_perm$silhouette)
        }
      }

      if ("d_prime" %in% matched_metrics) {
        sum_entry$p_value_d_prime <- (sum(perm_dp >= sum_entry$mean_d_prime) + 1) / (n_perm + 1)
      }
      if ("centroid_dist" %in% matched_metrics) {
        sum_entry$p_value_dist <- (sum(perm_cd >= sum_entry$mean_centroid_dist) + 1) / (n_perm + 1)
      }
      if ("bhattacharyya" %in% matched_metrics) {
        sum_entry$p_value_bhattacharyya <- (sum(perm_bh >= sum_entry$mean_bhattacharyya) + 1) / (n_perm + 1)
      }
      if ("silhouette" %in% matched_metrics) {
        sum_entry$p_value_silhouette <- (sum(perm_sl >= sum_entry$mean_silhouette) + 1) / (n_perm + 1)
      }
    }

    summary_list[[length(summary_list) + 1]] <- as.data.frame(sum_entry, stringsAsFactors = FALSE)
  }

  if (length(pointwise_list) == 0) {
    stop("No pairs had sufficient common coverage to calculate separation.")
  }

  pointwise_df <- do.call(rbind, pointwise_list)
  rownames(pointwise_df) <- NULL

  pairwise_summary_df <- do.call(rbind, summary_list)
  rownames(pairwise_summary_df) <- NULL

  # Omnibus summary across all pairs
  omnibus_summary <- data.frame(
    n_labels = n_labels,
    n_pairs = nrow(pairwise_summary_df),
    stringsAsFactors = FALSE
  )
  if ("d_prime" %in% matched_metrics) {
    omnibus_summary$mean_pairwise_d_prime <- mean(pairwise_summary_df$mean_d_prime)
    omnibus_summary$max_pairwise_d_prime <- max(pairwise_summary_df$mean_d_prime)
  }
  if ("centroid_dist" %in% matched_metrics) {
    omnibus_summary$mean_pairwise_dist <- mean(pairwise_summary_df$mean_centroid_dist)
  }
  if ("bhattacharyya" %in% matched_metrics) {
    omnibus_summary$mean_pairwise_bhattacharyya <- mean(pairwise_summary_df$mean_bhattacharyya)
  }
  if ("silhouette" %in% matched_metrics) {
    omnibus_summary$mean_pairwise_silhouette <- mean(pairwise_summary_df$mean_silhouette)
  }

  result <- list(
    pointwise = pointwise_df,
    pairwise_summary = pairwise_summary_df,
    omnibus = omnibus_summary,
    metrics = matched_metrics,
    dims = dims,
    type = "separation"
  )

  class(result) <- c("trajectory_separation", "list")
  result
}


#' @rdname trajectory_separation
#' @export
trajectory_separation.Sap <- function(x,
                                      segment_type = c("motifs", "syllables", "bouts", "segments"),
                                      dims = c("PC1", "PC2"),
                                      labels = NULL,
                                      metrics = c("d_prime", "centroid_dist", "bhattacharyya", "silhouette"),
                                      min_coverage = 0.5,
                                      time_digits = 6,
                                      stats = TRUE,
                                      n_perm = 500,
                                      seed = 222,
                                      verbose = TRUE,
                                      ...) {
  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  traj_embeds <- x$features[[feature_type]][["traj.embeds"]]
  if (is.null(traj_embeds)) {
    stop(sprintf(
      "No trajectory embeddings found for '%s'. Run create_trajectory_matrix() first.",
      segment_type
    ))
  }

  # Borrow missing dimensions from traj_mat if available
  missing_dims <- setdiff(dims, names(traj_embeds))
  if (length(missing_dims) > 0) {
    traj_mat <- x$features[[feature_type]][["traj_mat"]]
    dims_in_mat <- if (!is.null(traj_mat)) intersect(missing_dims, names(traj_mat)) else character(0)
    if (length(dims_in_mat) > 0 && nrow(traj_mat) == nrow(traj_embeds)) {
      if (verbose) {
        message(sprintf(
          "Dimensions [%s] not in traj.embeds; borrowing from traj_mat (run_pca() output).",
          paste(dims_in_mat, collapse = ", ")
        ))
      }
      traj_embeds[dims_in_mat] <- traj_mat[dims_in_mat]
    }
  }

  result <- trajectory_separation.default(
    x = traj_embeds,
    dims = dims,
    labels = labels,
    metrics = metrics,
    min_coverage = min_coverage,
    time_digits = time_digits,
    stats = stats,
    n_perm = n_perm,
    seed = seed,
    verbose = verbose,
    ...
  )

  x$features[[feature_type]][["trajectory_separation"]] <- result
  invisible(x)
}


# -------------------------------------------------------------------------
# Unified Discrimination Plotting
# -------------------------------------------------------------------------

#' Plot Trajectory Discrimination Results
#'
#' @description
#' Visualizes results produced by \code{\link{trajectory_timecourse}} or
#' \code{\link{trajectory_separation}}, displaying time-resolved variability
#' ribbons, divergence windows, and pairwise discriminability metrics.
#'
#' @param x A list returned by \code{trajectory_timecourse()} or
#'   \code{trajectory_separation()} (for the default method); or a SAP object
#'   (for the Sap method).
#' @param score_col Optional character specifying which separation metric to plot
#'   when visualizing \code{trajectory_separation()} results. One of \code{"d_prime"},
#'   \code{"centroid_dist"}, \code{"bhattacharyya"}, or \code{"silhouette"}.
#'   If \code{NULL} (default), uses the first available metric in the results.
#' @param palette RColorBrewer palette name (default: \code{"Set1"}).
#' @param segment_type For SAP objects: Type of segments to visualize
#'   (\code{'motifs'}, \code{'syllables'}, \code{'bouts'}, \code{'segments'}).
#' @param data_type For SAP objects: Which result to plot (\code{'timecourse'} or \code{'separation'}).
#' @param alpha_ribbon Numeric. Transparency of confidence interval ribbons (default: \code{0.25}).
#' @param ... Additional arguments passed to specific methods.
#'
#' @return An assembled \pkg{patchwork} object, returned invisibly.
#'
#' @examples
#' \dontrun{
#' # Visualizing variability timecourse
#' tc <- trajectory_timecourse(sap$features$motif$traj.embeds, dims = c("PC1", "PC2"))
#' plot_trajectory_discrimination(tc)
#'
#' # Visualizing trajectory separation (d-prime or Bhattacharyya distance)
#' sep <- trajectory_separation(sap$features$motif$traj.embeds, dims = c("PC1", "PC2"))
#' plot_trajectory_discrimination(sep, score_col = "d_prime")
#' plot_trajectory_discrimination(sep, score_col = "bhattacharyya")
#' }
#'
#' @export
plot_trajectory_discrimination <- function(x, ...) {
  UseMethod("plot_trajectory_discrimination")
}


#' @rdname plot_trajectory_discrimination
#' @export
plot_trajectory_discrimination.default <- function(x,
                                                   score_col = NULL,
                                                   palette = "Set1",
                                                   alpha_ribbon = 0.25,
                                                   ...) {
  result <- x
  if (!is.list(result) || is.null(result$type)) {
    stop("'x' must be a list returned by trajectory_timecourse() or trajectory_separation().")
  }

  ensure_pkgs("ggplot2", "patchwork")

  if (result$type == "timecourse") {
    tc <- result$timecourse
    summ <- result$summary
    tst <- result$tests

    labs <- sort_labels(unique(tc$label))
    pal_map <- make_pal(labs, palette)

    # Panel A: Timecourse ribbon plot
    p1 <- ggplot2::ggplot(tc, ggplot2::aes(x = .data$.time, y = .data$variability, color = .data$label, fill = .data$label)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
        alpha = alpha_ribbon,
        color = NA
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_fill_manual(values = pal_map) +
      ggplot2::labs(
        title = "Trajectory Variability Timecourse",
        subtitle = sprintf("Metric: %s | Dimensions: %s", result$metric, paste(result$dims, collapse = ", ")),
        x = "Time (seconds)",
        y = "Variability",
        color = "Condition",
        fill = "Condition"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "top"
      )

    # If significant pointwise intervals exist, highlight them
    if (!is.null(tst) && any(tst$is_significant)) {
      sig_times <- tst$.time[tst$is_significant]
      y_min <- min(tc$ci_lower, na.rm = TRUE)
      y_max <- max(tc$ci_upper, na.rm = TRUE)
      y_rug <- y_min - 0.05 * (y_max - y_min)

      p1 <- p1 +
        ggplot2::geom_point(
          data = data.frame(.time = sig_times, y_sig = y_rug),
          ggplot2::aes(x = .data$.time, y = .data$y_sig),
          inherit.aes = FALSE,
          shape = 124,
          size = 3,
          color = "red3"
        ) +
        ggplot2::labs(caption = "Red tick marks indicate time points with FDR-adjusted p < 0.05")
    }

    # Panel B: Summary bar plot of AUC variability
    p2 <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$label, y = .data$auc_variability, fill = .data$label)) +
      ggplot2::geom_col(width = 0.6, alpha = 0.85, show.legend = FALSE) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$auc_variability, ymax = .data$auc_variability),
        width = 0.2
      ) +
      ggplot2::scale_fill_manual(values = pal_map) +
      ggplot2::labs(
        title = "Total Integrated Variability (AUC)",
        x = "Condition",
        y = "Area Under Curve"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

    combined <- patchwork::wrap_plots(p1, p2, heights = c(2.2, 1))
    return(combined)

  } else if (result$type == "separation") {
    pw <- result$pointwise
    summ <- result$pairwise_summary

    pairs <- unique(pw$pair)
    pal_map <- make_pal(pairs, palette)

    # Determine metric column to plot
    metric_candidates <- c("d_prime", "centroid_dist", "bhattacharyya", "silhouette")
    available_metrics <- intersect(metric_candidates, names(pw))

    if (is.null(score_col)) {
      score_col <- available_metrics[1]
    } else if (!score_col %in% available_metrics) {
      stop(sprintf(
        "Requested score_col '%s' not found in separation results. Available: %s",
        score_col, paste(available_metrics, collapse = ", ")
      ))
    }

    # Label styling by metric
    metric_labels <- list(
      d_prime = list(name = "d-prime", title = "Pointwise Trajectory Discriminability (d')", summary_col = "mean_d_prime", p_col = "p_value_d_prime"),
      centroid_dist = list(name = "Centroid Distance", title = "Pointwise Centroid Distance", summary_col = "mean_centroid_dist", p_col = "p_value_dist"),
      bhattacharyya = list(name = "Bhattacharyya Distance", title = "Pointwise Bhattacharyya Divergence", summary_col = "mean_bhattacharyya", p_col = "p_value_bhattacharyya"),
      silhouette = list(name = "Silhouette Score", title = "Pointwise Trajectory Silhouette Score", summary_col = "mean_silhouette", p_col = "p_value_silhouette")
    )
    meta <- metric_labels[[score_col]]

    # Panel A: Pointwise profile
    p1 <- ggplot2::ggplot(pw, ggplot2::aes(x = .data$.time, y = .data[[score_col]], color = .data$pair)) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::labs(
        title = meta$title,
        subtitle = sprintf("Dimensions: %s", paste(result$dims, collapse = ", ")),
        x = "Time (seconds)",
        y = meta$name,
        color = "Pair"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "top"
      )

    # Panel B: Mean score per pair with permutation p-values
    p_col_name <- meta$p_col
    summ$p_label <- if (p_col_name %in% names(summ) && !anyNA(summ[[p_col_name]])) {
      vapply(summ[[p_col_name]], fmt_p, character(1))
    } else {
      ""
    }

    sum_val_col <- meta$summary_col
    p2 <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$pair, y = .data[[sum_val_col]], fill = .data$pair)) +
      ggplot2::geom_col(width = 0.6, alpha = 0.85, show.legend = FALSE) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$p_label),
        vjust = -0.5,
        size = 3.5
      ) +
      ggplot2::scale_fill_manual(values = pal_map) +
      ggplot2::labs(
        title = sprintf("Mean %s", meta$name),
        x = "Pair",
        y = meta$name
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

    combined <- patchwork::wrap_plots(p1, p2, heights = c(2.2, 1))
    return(combined)
  } else {
    stop(sprintf("Unknown discrimination result type '%s'.", result$type))
  }
}


#' @rdname plot_trajectory_discrimination
#' @export
plot_trajectory_discrimination.Sap <- function(x,
                                               segment_type = c("motifs", "syllables", "bouts", "segments"),
                                               data_type = c("timecourse", "separation"),
                                               score_col = NULL,
                                               palette = "Set1",
                                               alpha_ribbon = 0.25,
                                               ...) {
  segment_type <- match.arg(segment_type)
  data_type <- match.arg(data_type)
  feature_type <- sub("s$", "", segment_type)

  target_name <- if (data_type == "timecourse") "trajectory_timecourse" else "trajectory_separation"
  result <- x$features[[feature_type]][[target_name]]

  if (is.null(result)) {
    stop(sprintf(
      "No %s results found for '%s'. Run %s() first.",
      data_type, segment_type, target_name
    ))
  }

  plot_trajectory_discrimination.default(
    x = result,
    score_col = score_col,
    palette = palette,
    alpha_ribbon = alpha_ribbon,
    ...
  )
}
