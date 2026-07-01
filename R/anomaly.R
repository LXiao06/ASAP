# Anomaly Detection
# Update date : Jul. 01, 2026
# Identifies time points with abnormal variability, maturation drops, or
# multivariate outliers using Mahalanobis distance.

# Global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
  "n_samples", "is_anomalous",
  "label_numeric", "method", "flagged",
  "total_flags",
  "orthogonal_rms", "parallel_rms", "total_rms",
  "orthogonal_rms_scaled", "parallel_rms_scaled", "total_rms_scaled",
  "mean_orthogonal_rms", "mean_parallel_rms", "mean_total_rms",
  "mean_orthogonal_rms_scaled", "mean_parallel_rms_scaled", "mean_total_rms_scaled",
  "mean_rms_similarity", "mean_correlation", "rms_similarity", "correlation",
  "trajectory_score", "trajectory_score_signed", "trajectory_metric", "anomaly_score",
  "sd_rms_similarity", "sd_correlation",
  "metric_name", "label_cat_metric",
  "dispersion_cv", "path_length_cv", "lof_score", "row_group"
))


# Detection -------------------------------------------------------------------

#' Detect Anomalous Time Labels
#'
#' @description
#' Identifies labels with abnormal variability, low maturation scores,
#' or unexpected deviations from trajectory trends. Useful for detecting
#' developmental anomalies or time points requiring further investigation.
#'
#' @param sap SAP object with computed maturation scores.
#' @param segment_type Type of segments: \code{"motifs"}, \code{"syllables"},
#'   \code{"bouts"}, or \code{"segments"}.
#' @param methods Vector of detection methods to use:
#'   \itemize{
#'     \item \code{"trajectory"}: Fits a loess trend to each metric and
#'       z-scores the residuals. Flags labels with |z-score| > extreme_threshold.
#'     \item \code{"multivariate"}: Mahalanobis distance on the loess residual
#'       z-scores. Falls back to raw feature means when trajectory is not run.
#'       Requires MASS package.
#'     \item \code{"lof"}: Local Outlier Factor (non-parametric). Requires
#'       the \code{dbscan} package.
#'   }
#' @param extreme_threshold Z-score threshold for detecting extreme values
#'   (default: 2). Flags when |z-score| > extreme_threshold.
#' @param isolation_quantile Chi-squared quantile for multivariate threshold
#'   (default: 0.99). Increase to flag fewer labels; decrease to flag more.
#' @param min_sample_size Minimum samples per label to include (default: 3).
#' @param lof_k Integer. Number of nearest neighbours for LOF. Defaults to
#'   \code{max(3, min(10, floor(n/4)))} where \code{n} is the number of labels.
#' @param cov_regularization Numeric. Regularization constant for covariance
#'   shrinkage (default: 0.5).
#' @param verbose Logical. Print progress messages (default: TRUE).
#'
#' @return A list of class \code{"label_anomaly"} with:
#'   \itemize{
#'     \item \code{label_scores}: Data frame with anomaly scores per label.
#'     \item \code{anomalous_labels}: Vector of flagged labels.
#'     \item \code{method_flags}: Data frame with per-method flags
#'       (\code{trajectory_flag}, \code{multivariate_flag}, \code{lof_flag},
#'       \code{total_flags}, \code{is_anomalous}).
#'     \item \code{parameters}: List of detection parameters used.
#'   }
#'
#' @examples
#' \dontrun{
#' # Detect anomalous labels using all three methods (default)
#' anomaly_results <- detect_anomalous_labels(sap, segment_type = "motifs")
#'
#' # View flagged labels
#' print(anomaly_results$anomalous_labels)
#'
#' # View detailed scores
#' View(anomaly_results$label_scores)
#'
#' # Trajectory method only
#' anomaly_results <- detect_anomalous_labels(
#'   sap,
#'   segment_type = "motifs",
#'   methods = "trajectory"
#' )
#'
#' # Custom LOF neighbourhood size
#' anomaly_results <- detect_anomalous_labels(
#'   sap,
#'   segment_type = "motifs",
#'   lof_k = 5
#' )
#' }
#'
#' @seealso \code{\link{plot_label_anomaly}} for visualising results.
#'
#' @export
detect_anomalous_labels <- function(sap,
                                    segment_type = c("motifs", "syllables",
                                                     "bouts", "segments"),
                                    methods = c("trajectory", "multivariate",
                                                "lof"),
                                    extreme_threshold  = 2,
                                    isolation_quantile = 0.99,
                                    min_sample_size    = 3,
                                    lof_k              = NULL,
                                    cov_regularization = 0.5,
                                    verbose            = TRUE) {

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  # Validate feature slot
  if (is.null(sap$features[[feature_type]])) {
    stop(sprintf("No features found for '%s' in SAP object.", feature_type))
  }

  scores <- sap$features[[feature_type]][["maturation_scores"]]
  if (is.null(scores)) {
    stop(sprintf(
      "Maturation scores not found for %s. Run trajectory_maturation() first.",
      segment_type
    ))
  }

  # Get supplementary data if available
  dispersion <- sap$features[[feature_type]][["dispersion_comparison"]]
  path_dev   <- sap$features[[feature_type]][["trajectory_path_deviation"]]

  if (verbose) {
    message("\n=== Detecting Anomalous Labels ===")
    message(sprintf("Segment type       : %s", segment_type))
    message(sprintf("Methods            : %s", paste(methods, collapse = ", ")))
    message(sprintf("Min samples/label  : %d", min_sample_size))
    message("")
  }

  # Aggregate metrics by label
  label_stats <- scores |>
    dplyr::group_by(label) |>
    dplyr::summarise(
      n_samples           = dplyr::n(),
      mean_rms_similarity = mean(rms_similarity, na.rm = TRUE),
      sd_rms_similarity   = stats::sd(rms_similarity, na.rm = TRUE),
      mean_correlation = if ("correlation" %in% names(scores)) {
        mean(correlation, na.rm = TRUE)
      } else {
        NA_real_
      },
      sd_correlation = if ("correlation" %in% names(scores)) {
        stats::sd(correlation, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::filter(n_samples >= min_sample_size)

  # Merge dispersion metrics
  if (!is.null(dispersion)) {
    label_stats <- label_stats |>
      dplyr::left_join(
        dispersion |> dplyr::select(label, dispersion_cv, path_length_cv),
        by = "label"
      )
  }

  # Merge path deviation metrics
  if (!is.null(path_dev) && !is.null(path_dev$width)) {
    has_scaled_orth  <- "orthogonal_rms_scaled" %in% names(path_dev$width)
    has_scaled_par   <- "parallel_rms_scaled"   %in% names(path_dev$width)
    has_total        <- "total_rms"             %in% names(path_dev$width)
    has_total_scaled <- "total_rms_scaled"      %in% names(path_dev$width)

    path_metrics <- path_dev$width |>
      dplyr::group_by(label) |>
      dplyr::summarise(
        mean_orthogonal_rms = mean(orthogonal_rms, na.rm = TRUE),
        mean_parallel_rms   = mean(parallel_rms,   na.rm = TRUE),
        sd_orthogonal_rms   = stats::sd(orthogonal_rms, na.rm = TRUE),
        sd_parallel_rms     = stats::sd(parallel_rms,   na.rm = TRUE),
        .groups = "drop"
      )

    if (has_scaled_orth) {
      temp <- path_dev$width |>
        dplyr::group_by(label) |>
        dplyr::summarise(
          mean_orthogonal_rms_scaled = mean(orthogonal_rms_scaled, na.rm = TRUE),
          .groups = "drop"
        )
      path_metrics <- dplyr::left_join(path_metrics, temp, by = "label")
    }

    if (has_scaled_par) {
      temp <- path_dev$width |>
        dplyr::group_by(label) |>
        dplyr::summarise(
          mean_parallel_rms_scaled = mean(parallel_rms_scaled, na.rm = TRUE),
          .groups = "drop"
        )
      path_metrics <- dplyr::left_join(path_metrics, temp, by = "label")
    }

    if (has_total) {
      temp <- path_dev$width |>
        dplyr::group_by(label) |>
        dplyr::summarise(
          mean_total_rms = mean(total_rms, na.rm = TRUE),
          .groups = "drop"
        )
      path_metrics <- dplyr::left_join(path_metrics, temp, by = "label")
    }

    if (has_total_scaled) {
      temp <- path_dev$width |>
        dplyr::group_by(label) |>
        dplyr::summarise(
          mean_total_rms_scaled = mean(total_rms_scaled, na.rm = TRUE),
          .groups = "drop"
        )
      path_metrics <- dplyr::left_join(path_metrics, temp, by = "label")
    }

    label_stats <- label_stats |>
      dplyr::left_join(path_metrics, by = "label")
  }

  # Initialise per-method flag table
  flag_df <- data.frame(
    label             = label_stats$label,
    trajectory_flag   = FALSE,
    multivariate_flag = FALSE,
    lof_flag          = FALSE,
    stringsAsFactors  = FALSE
  )

  # Loess residual z-scores shared between trajectory, multivariate, and LOF
  all_zscores_for_mv <- NULL


  # Method 1: Trajectory Extremes (Similarity + Variability) ----------------

  if ("trajectory" %in% methods) {
    # Similarity metrics (expect to increase; drops are anomalies)
    similarity_metrics <- intersect(
      c("mean_rms_similarity", "mean_correlation"),
      names(label_stats)
    )

    # Variability metrics (expect to decrease; increases are anomalies)
    variability_metrics <- intersect(
      c("mean_orthogonal_rms_scaled", "mean_parallel_rms_scaled",
        "mean_total_rms_scaled",
        "mean_orthogonal_rms", "mean_parallel_rms", "mean_total_rms"),
      names(label_stats)
    )

    all_metrics <- c(similarity_metrics, variability_metrics)

    if (length(all_metrics) > 0) {
      label_stats   <- dplyr::arrange(label_stats, label)
      label_numeric <- as.numeric(as.character(label_stats$label))

      if (nrow(label_stats) > 5 && all(!is.na(label_numeric))) {
        # Retain only metrics with at least 3 non-NA values
        valid_metrics <- all_metrics[sapply(all_metrics, function(m) {
          sum(!is.na(label_stats[[m]])) >= 3
        })]

        if (length(valid_metrics) > 0) {
          all_zscores <- sapply(valid_metrics, function(m) {
            vals     <- label_stats[[m]]
            fit_data <- data.frame(vals = vals, label_numeric = label_numeric)
            fit_data <- fit_data[stats::complete.cases(fit_data), ]

            if (nrow(fit_data) < 3) {
              return(rep(NA_real_, nrow(label_stats)))
            }

            smooth_fit <- stats::loess(vals ~ label_numeric,
                                       data = fit_data, span = 0.3)

            pred <- rep(NA_real_, nrow(label_stats))
            pred[stats::complete.cases(
              data.frame(vals, label_numeric)
            )] <- stats::predict(smooth_fit)

            resid <- vals - pred

            if (sum(!is.na(resid)) >= 3) {
              zscore <- (resid - mean(resid, na.rm = TRUE)) /
                stats::sd(resid, na.rm = TRUE)
            } else {
              zscore <- rep(NA_real_, length(resid))
            }

            # Flip sign for similarity so drops produce positive z-scores,
            # consistent with variability where increases = positive
            if (m %in% similarity_metrics) zscore <- -zscore

            zscore
          })

          # Share residual z-score matrix with multivariate and LOF
          all_zscores_for_mv <- all_zscores

          # Max absolute z-score across all metrics per label
          label_stats$trajectory_score <- apply(
            abs(all_zscores), 1, max, na.rm = TRUE
          )

          # Signed score for plotting (shows direction of worst deviation)
          label_stats$trajectory_score_signed <- apply(
            all_zscores, 1, function(row) row[which.max(abs(row))]
          )

          # Which metric drove the max deviation
          label_stats$trajectory_metric <- apply(
            all_zscores, 1, function(row) {
              if (all(is.na(row))) return(NA_character_)
              valid_metrics[which.max(abs(row))]
            }
          )

          flag_df$trajectory_flag <-
            !is.na(label_stats$trajectory_score) &
            label_stats$trajectory_score > extreme_threshold

          # Per-metric flags for the heatmap
          for (metric in valid_metrics) {
            metric_col    <- all_zscores[, metric]
            flag_col_name <- paste0(metric, "_flag")
            flag_df[[flag_col_name]] <-
              !is.na(metric_col) & abs(metric_col) > extreme_threshold
            label_stats[[paste0(metric, "_zscore")]] <- metric_col
          }

          if (verbose) {
            message(sprintf(
              "Trajectory extremes: %d labels flagged (|z-score| > %.1f, using %d metrics)",
              sum(flag_df$trajectory_flag, na.rm = TRUE),
              extreme_threshold,
              length(valid_metrics)
            ))
            message(sprintf(
              "  Similarity metrics: %s",
              ifelse(
                length(similarity_metrics) > 0,
                paste(similarity_metrics, collapse = ", "),
                "none"
              )
            ))
            message(sprintf(
              "  Variability metrics: %s",
              ifelse(
                length(variability_metrics) > 0,
                paste(variability_metrics, collapse = ", "),
                "none"
              )
            ))
          }
        } else if (verbose) {
          message("Trajectory extremes: No valid metrics with sufficient data")
        }
      } else if (verbose) {
        message(
          "Trajectory extremes: Too few labels for trend fitting or labels not numeric"
        )
      }
    } else if (verbose) {
      message("Trajectory extremes: No trajectory metrics available")
    }
  }


  # Method 2: Multivariate (Mahalanobis on loess residuals) -----------------
  #
  # Operates on the loess residual z-scores from the trajectory method rather
  # than raw per-label means. Removing the developmental trend makes the
  # residuals closer to i.i.d., better satisfying the multivariate-normal
  # assumption behind the chi-squared threshold. Scaled and unscaled versions
  # of the same metric are collapsed to the scaled version to prevent
  # near-singular covariance estimation.

  if ("multivariate" %in% methods) {
    if (!is.null(all_zscores_for_mv) && ncol(all_zscores_for_mv) >= 2) {
      zscore_mat <- all_zscores_for_mv

      # Prefer scaled when both scaled and unscaled versions exist
      drop_unscaled <- character(0)
      scale_pairs <- list(
        c("mean_orthogonal_rms_scaled", "mean_orthogonal_rms"),
        c("mean_parallel_rms_scaled",   "mean_parallel_rms"),
        c("mean_total_rms_scaled",       "mean_total_rms")
      )
      for (pair in scale_pairs) {
        if (all(pair %in% colnames(zscore_mat))) {
          drop_unscaled <- c(drop_unscaled, pair[2])
        }
      }
      if (length(drop_unscaled) > 0) {
        zscore_mat <- zscore_mat[
          , !colnames(zscore_mat) %in% drop_unscaled, drop = FALSE
        ]
      }

      complete_idx   <- stats::complete.cases(zscore_mat)
      features_clean <- zscore_mat[complete_idx, , drop = FALSE]
      p              <- ncol(features_clean)

      # Need at least 2p samples for reliable covariance estimation
      if (nrow(features_clean) >= max(5, 2 * p)) {
        cov_result <- tryCatch(
          MASS::cov.mcd(features_clean),
          error = function(e) NULL
        )

        if (!is.null(cov_result)) {
          center  <- cov_result$center
          cov_mat <- cov_result$cov
        } else {
          center  <- colMeans(features_clean)
          cov_mat <- stats::cov(features_clean)
        }

        # Regularize covariance for invertibility and scale stability
        diag_val <- mean(diag(cov_mat))
        if (is.na(diag_val) || diag_val == 0) diag_val <- 1.0
        diag(cov_mat) <- diag(cov_mat) + cov_regularization * diag_val

        md <- stats::mahalanobis(features_clean, center, cov_mat)
        label_stats$anomaly_score <- NA_real_
        label_stats$anomaly_score[complete_idx] <- md

        chi_threshold <- stats::qchisq(isolation_quantile, df = p)
        flag_df$multivariate_flag <-
          !is.na(label_stats$anomaly_score) &
          label_stats$anomaly_score > chi_threshold

        if (verbose) {
          message(sprintf(
            "Multivariate: %d labels flagged (Mahalanobis on residuals, %d metrics)",
            sum(flag_df$multivariate_flag, na.rm = TRUE), p
          ))
        }
      } else if (verbose) {
        message("Multivariate: Too few labels for reliable covariance estimation")
      }
    } else {
      # Fall back to raw feature means when trajectory method was not run
      feature_cols <- intersect(
        c("mean_rms_similarity", "mean_correlation",
          "mean_orthogonal_rms_scaled", "mean_parallel_rms_scaled",
          "mean_total_rms_scaled",
          "mean_orthogonal_rms", "mean_parallel_rms", "mean_total_rms"),
        names(label_stats)
      )

      # Prefer scaled when both exist
      for (pair in list(
        c("mean_orthogonal_rms_scaled", "mean_orthogonal_rms"),
        c("mean_parallel_rms_scaled",   "mean_parallel_rms"),
        c("mean_total_rms_scaled",       "mean_total_rms")
      )) {
        if (all(pair %in% feature_cols)) {
          feature_cols <- setdiff(feature_cols, pair[2])
        }
      }

      if (length(feature_cols) >= 2) {
        features       <- label_stats[, feature_cols, drop = FALSE]
        complete_idx   <- stats::complete.cases(features)
        features_clean <- features[complete_idx, , drop = FALSE]
        p              <- length(feature_cols)

        if (nrow(features_clean) >= max(5, 2 * p)) {
          cov_result <- tryCatch(
            MASS::cov.mcd(features_clean),
            error = function(e) NULL
          )

          if (!is.null(cov_result)) {
            center  <- cov_result$center
            cov_mat <- cov_result$cov
          } else {
            center  <- colMeans(features_clean)
            cov_mat <- stats::cov(features_clean)
          }

          # Regularize covariance for invertibility and scale stability
          diag_val <- mean(diag(cov_mat))
          if (is.na(diag_val) || diag_val == 0) diag_val <- 1.0
          diag(cov_mat) <- diag(cov_mat) + cov_regularization * diag_val

          md <- stats::mahalanobis(features_clean, center, cov_mat)
          label_stats$anomaly_score <- NA_real_
          label_stats$anomaly_score[complete_idx] <- md

          chi_threshold <- stats::qchisq(isolation_quantile, df = p)
          flag_df$multivariate_flag <-
            !is.na(label_stats$anomaly_score) &
            label_stats$anomaly_score > chi_threshold

          if (verbose) {
            message(sprintf(
              "Multivariate: %d labels flagged (Mahalanobis on raw features, %d features)",
              sum(flag_df$multivariate_flag, na.rm = TRUE), p
            ))
          }
        } else if (verbose) {
          message("Multivariate: Too few labels for reliable covariance estimation")
        }
      } else if (verbose) {
        message("Multivariate: Not enough features available (need >= 2)")
      }
    }
  }


  # Method 3: LOF (Local Outlier Factor) ------------------------------------
  #
  # Non-parametric density-based method -- no multivariate normality required.
  # Detects labels whose feature neighbourhood is significantly less dense
  # than their neighbours (LOF score >> 1). Uses loess residual z-scores when
  # available (consistent with the multivariate method) and falls back to raw
  # per-label means otherwise.

  if ("lof" %in% methods) {
    ensure_pkgs("dbscan")

    # Prefer residual matrix; fall back to raw feature means
    lof_mat <- if (!is.null(all_zscores_for_mv) &&
                   ncol(all_zscores_for_mv) >= 2) {
      all_zscores_for_mv
    } else {
      feat_cols <- intersect(
        c("mean_rms_similarity", "mean_correlation",
          "mean_orthogonal_rms_scaled", "mean_parallel_rms_scaled",
          "mean_total_rms_scaled",
          "mean_orthogonal_rms", "mean_parallel_rms", "mean_total_rms"),
        names(label_stats)
      )
      for (pair in list(
        c("mean_orthogonal_rms_scaled", "mean_orthogonal_rms"),
        c("mean_parallel_rms_scaled",   "mean_parallel_rms"),
        c("mean_total_rms_scaled",       "mean_total_rms")
      )) {
        if (all(pair %in% feat_cols)) feat_cols <- setdiff(feat_cols, pair[2])
      }
      if (length(feat_cols) >= 2) {
        as.matrix(label_stats[, feat_cols, drop = FALSE])
      } else {
        NULL
      }
    }

    if (!is.null(lof_mat) && ncol(lof_mat) >= 2) {
      complete_idx_lof <- stats::complete.cases(lof_mat)
      lof_mat_clean    <- lof_mat[complete_idx_lof, , drop = FALSE]
      n_clean          <- nrow(lof_mat_clean)

      if (n_clean >= 5) {
        # k defaults to ~25% of labels, clamped to [3, 10]
        k_lof <- if (!is.null(lof_k)) {
          as.integer(lof_k)
        } else {
          max(3L, min(10L, as.integer(floor(n_clean / 4))))
        }
        k_lof <- min(k_lof, n_clean - 1L)

        lof_scores_clean <- dbscan::lof(lof_mat_clean, minPts = k_lof)

        label_stats$lof_score <- NA_real_
        label_stats$lof_score[complete_idx_lof] <- lof_scores_clean

        # LOF > 1.5 is the conventional heuristic for outliers
        lof_threshold    <- 1.5
        flag_df$lof_flag <-
          !is.na(label_stats$lof_score) &
          label_stats$lof_score > lof_threshold

        if (verbose) {
          message(sprintf(
            "LOF: %d labels flagged (LOF > %.1f, k = %d)",
            sum(flag_df$lof_flag, na.rm = TRUE), lof_threshold, k_lof
          ))
        }
      } else if (verbose) {
        message("LOF: Too few labels for reliable LOF estimation (need >= 5)")
      }
    } else if (verbose) {
      message("LOF: Not enough features available (need >= 2)")
    }
  }


  # Combine flags -----------------------------------------------------------

  flag_col_names       <- grep("_flag$", names(flag_df), value = TRUE)
  flag_df$total_flags  <- rowSums(flag_df[, flag_col_names, drop = FALSE])
  flag_df$is_anomalous <- flag_df$total_flags > 0

  label_scores <- dplyr::left_join(label_stats, flag_df, by = "label")

  anomalous_labels <- label_scores |>
    dplyr::filter(is_anomalous) |>
    dplyr::pull(label)

  if (verbose) {
    message(sprintf("\nTotal anomalous labels: %d", length(anomalous_labels)))
    if (length(anomalous_labels) > 0) {
      message(sprintf("Labels: %s", paste(anomalous_labels, collapse = ", ")))
    }
    message("")
  }

  result <- list(
    label_scores     = label_scores,
    anomalous_labels = anomalous_labels,
    method_flags     = flag_df,
    parameters       = list(
      segment_type       = segment_type,
      methods            = methods,
      extreme_threshold  = extreme_threshold,
      isolation_quantile = isolation_quantile,
      min_sample_size    = min_sample_size,
      lof_k              = lof_k,
      cov_regularization = cov_regularization
    )
  )

  class(result) <- c("label_anomaly", class(result))
  return(result)
}


# Plotting --------------------------------------------------------------------

#' Plot Label Anomaly Detection Results
#'
#' @description
#' Visualises label anomaly detection results, highlighting time points flagged
#' as having abnormal variability or maturation patterns. When
#' \code{heatmap_only = FALSE} (default), a multi-panel figure of individual
#' metric traces is returned (heatmap excluded). When \code{heatmap_only = TRUE},
#' a standalone summary heatmap is returned instead.
#'
#' @param detection_result Output from \code{\link{detect_anomalous_labels}}.
#' @param sap Original SAP object (reserved for future enhancements).
#' @param highlight_color Color for anomalous labels (default: \code{"red"}).
#' @param heatmap_only Logical. If \code{TRUE}, returns only the flag heatmap
#'   as a standalone summary figure. If \code{FALSE} (default), returns a
#'   multi-panel figure of individual metric traces without the heatmap.
#'
#' @return A \code{ggplot} or \code{patchwork} object.
#'
#' @examples
#' \dontrun{
#' anomaly_results <- detect_anomalous_labels(sap, segment_type = "motifs")
#'
#' # Multi-panel metric traces
#' plot_label_anomaly(anomaly_results, sap)
#'
#' # Standalone heatmap summary
#' plot_label_anomaly(anomaly_results, sap, heatmap_only = TRUE)
#' }
#'
#' @seealso \code{\link{detect_anomalous_labels}}
#'
#' @export
plot_label_anomaly <- function(detection_result,
                               sap             = NULL,
                               highlight_color = "red",
                               heatmap_only    = FALSE) {

  if (!inherits(detection_result, "label_anomaly")) {
    stop("detection_result must be output from detect_anomalous_labels()")
  }

  ensure_pkgs("ggplot2", "patchwork", "scales")

  label_scores <- detection_result$label_scores
  label_col    <- names(label_scores)[1]

  # Add numeric label column for continuous x-axis
  plot_data <- label_scores |>
    dplyr::mutate(
      label_numeric = as.numeric(as.character(!!rlang::sym(label_col)))
    )

  # Heatmap-only path: build and return standalone summary figure
  if (heatmap_only) {
    p_heat <- build_anomaly_heatmap(
      plot_data        = plot_data,
      detection_result = detection_result,
      label_col        = label_col,
      highlight_color  = highlight_color
    )
    if (is.null(p_heat)) stop("No flags were generated to plot in heatmap.")

    methods_used    <- paste(detection_result$parameters$methods, collapse = ", ")
    total_anomalous <- length(detection_result$anomalous_labels)
    p_heat <- p_heat + ggplot2::labs(
      title    = "Anomalous Label Detection Summary Heatmap",
      subtitle = sprintf(
        "Methods: %s  |  Anomalous: %d / %d labels",
        methods_used, total_anomalous, nrow(label_scores)
      )
    )
    return(p_heat)
  }

  # Multi-panel metric traces (heatmap excluded from combined figure)
  plots <- list()


  # Plot 1: RMS Similarity ---------------------------------------------------

  if ("mean_rms_similarity" %in% names(plot_data)) {
    plot_data$label_cat_metric <- make_metric_anomaly_factor(
      plot_data, "mean_rms_similarity"
    )

    p1 <- ggplot2::ggplot(
      plot_data, ggplot2::aes(x = label_numeric, y = mean_rms_similarity)
    ) +
      ggplot2::geom_smooth(
        method = "loess", span = 0.3, se = TRUE,
        color = "gray50", fill = "gray80", alpha = 0.2, linewidth = 0.5
      ) +
      ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(color = label_cat_metric, size = label_cat_metric)
      ) +
      ggplot2::scale_color_manual(
        values = c("Normal" = "steelblue", "Anomalous" = highlight_color)
      ) +
      ggplot2::scale_size_manual(
        values = c("Normal" = 2, "Anomalous" = 4), guide = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      ggplot2::labs(
        title = "RMS Similarity",
        x     = label_col,
        y     = "Mean RMS Similarity",
        color = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")
    plots <- c(plots, list(p1))
  }


  # Plot 2: Correlation -------------------------------------------------------

  if ("mean_correlation" %in% names(plot_data) &&
      !all(is.na(plot_data$mean_correlation))) {
    plot_data$label_cat_metric <- make_metric_anomaly_factor(
      plot_data, "mean_correlation"
    )

    p2 <- ggplot2::ggplot(
      plot_data, ggplot2::aes(x = label_numeric, y = mean_correlation)
    ) +
      ggplot2::geom_smooth(
        method = "loess", span = 0.3, se = TRUE,
        color = "gray50", fill = "gray80", alpha = 0.2, linewidth = 0.5
      ) +
      ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(color = label_cat_metric, size = label_cat_metric)
      ) +
      ggplot2::scale_color_manual(
        values = c("Normal" = "steelblue", "Anomalous" = highlight_color)
      ) +
      ggplot2::scale_size_manual(
        values = c("Normal" = 2, "Anomalous" = 4), guide = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      ggplot2::labs(
        title = "Correlation",
        x     = label_col,
        y     = "Mean Correlation",
        color = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")
    plots <- c(plots, list(p2))
  }


  # Plot 3: Orthogonal RMS ---------------------------------------------------

  orth_col <- if ("mean_orthogonal_rms_scaled" %in% names(plot_data)) {
    "mean_orthogonal_rms_scaled"
  } else if ("mean_orthogonal_rms" %in% names(plot_data)) {
    "mean_orthogonal_rms"
  } else {
    NULL
  }

  if (!is.null(orth_col)) {
    plot_data$label_cat_metric <- make_metric_anomaly_factor(
      plot_data, orth_col
    )

    p3 <- ggplot2::ggplot(
      plot_data, ggplot2::aes(x = label_numeric, y = .data[[orth_col]])
    ) +
      ggplot2::geom_smooth(
        method = "loess", span = 0.3, se = TRUE,
        color = "gray50", fill = "gray80", alpha = 0.2, linewidth = 0.5
      ) +
      ggplot2::geom_line(color = "black", linewidth = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(color = label_cat_metric, size = label_cat_metric)
      ) +
      ggplot2::scale_color_manual(
        values = c("Normal" = "black", "Anomalous" = highlight_color)
      ) +
      ggplot2::scale_size_manual(
        values = c("Normal" = 2, "Anomalous" = 4), guide = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      ggplot2::labs(
        title = ifelse(
          orth_col == "mean_orthogonal_rms_scaled",
          "Orthogonal RMS (Scaled)", "Orthogonal RMS"
        ),
        x     = label_col,
        y     = "Mean Orthogonal RMS",
        color = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")
    plots <- c(plots, list(p3))
  }


  # Plot 4: Parallel RMS -----------------------------------------------------

  par_col <- if ("mean_parallel_rms_scaled" %in% names(plot_data)) {
    "mean_parallel_rms_scaled"
  } else if ("mean_parallel_rms" %in% names(plot_data)) {
    "mean_parallel_rms"
  } else {
    NULL
  }

  if (!is.null(par_col)) {
    plot_data$label_cat_metric <- make_metric_anomaly_factor(
      plot_data, par_col
    )

    p4 <- ggplot2::ggplot(
      plot_data, ggplot2::aes(x = label_numeric, y = .data[[par_col]])
    ) +
      ggplot2::geom_smooth(
        method = "loess", span = 0.3, se = TRUE,
        color = "gray50", fill = "gray80", alpha = 0.2, linewidth = 0.5
      ) +
      ggplot2::geom_line(color = "black", linewidth = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(color = label_cat_metric, size = label_cat_metric)
      ) +
      ggplot2::scale_color_manual(
        values = c("Normal" = "black", "Anomalous" = highlight_color)
      ) +
      ggplot2::scale_size_manual(
        values = c("Normal" = 2, "Anomalous" = 4), guide = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      ggplot2::labs(
        title = ifelse(
          par_col == "mean_parallel_rms_scaled",
          "Parallel RMS (Scaled)", "Parallel RMS"
        ),
        x     = label_col,
        y     = "Mean Parallel RMS",
        color = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")
    plots <- c(plots, list(p4))
  }


  # Plot 5: Trajectory Deviation Score ---------------------------------------
  # Shows the signed max-residual z-score across all metrics. Threshold lines
  # indicate the extreme_threshold used for flagging.

  if ("trajectory_score_signed" %in% names(plot_data)) {
    traj_vals <- detection_result$method_flags$trajectory_flag
    traj_vals[is.na(traj_vals)] <- FALSE
    plot_data$label_cat_metric <- factor(
      traj_vals,
      levels = c(FALSE, TRUE),
      labels = c("Normal", "Anomalous")
    )

    p5 <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = label_numeric, y = trajectory_score_signed)
    ) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5
      ) +
      ggplot2::geom_hline(
        yintercept = c(
          -detection_result$parameters$extreme_threshold,
           detection_result$parameters$extreme_threshold
        ),
        linetype = "dashed", color = highlight_color, alpha = 0.3
      ) +
      ggplot2::geom_line(color = "purple", linewidth = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(color = label_cat_metric, size = label_cat_metric)
      ) +
      ggplot2::scale_color_manual(
        values = c("Normal" = "purple", "Anomalous" = highlight_color)
      ) +
      ggplot2::scale_size_manual(
        values = c("Normal" = 2, "Anomalous" = 4), guide = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      ggplot2::labs(
        title    = "Trajectory Deviation Score",
        subtitle = sprintf(
          "Max loess residual z-score across all metrics  |  threshold \u00b1%.1f",
          detection_result$parameters$extreme_threshold
        ),
        x     = label_col,
        y     = "Signed Z-score",
        color = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")
    plots <- c(plots, list(p5))
  }


  # Plot 6: LOF Score --------------------------------------------------------
  # LOF operates on the same residual space as Mahalanobis, so surfacing the
  # raw LOF score gives complementary non-parametric information.

  if ("lof_score" %in% names(plot_data) && !all(is.na(plot_data$lof_score))) {
    lof_vals <- detection_result$method_flags$lof_flag
    lof_vals[is.na(lof_vals)] <- FALSE
    plot_data$label_cat_metric <- factor(
      lof_vals,
      levels = c(FALSE, TRUE),
      labels = c("Normal", "Anomalous")
    )

    p6 <- ggplot2::ggplot(
      plot_data, ggplot2::aes(x = label_numeric, y = lof_score)
    ) +
      ggplot2::geom_hline(
        yintercept = 1, linetype = "dotted", color = "gray60", alpha = 0.6
      ) +
      ggplot2::geom_hline(
        yintercept = 1.5, linetype = "dashed", color = highlight_color, alpha = 0.3
      ) +
      ggplot2::geom_line(color = "#e07b39", linewidth = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(color = label_cat_metric, size = label_cat_metric)
      ) +
      ggplot2::scale_color_manual(
        values = c("Normal" = "#e07b39", "Anomalous" = highlight_color)
      ) +
      ggplot2::scale_size_manual(
        values = c("Normal" = 2, "Anomalous" = 4), guide = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      ggplot2::labs(
        title    = "LOF Score (Local Outlier Factor)",
        subtitle = "Threshold > 1.5 (dashed)  |  score = 1 means same density as neighbours",
        x        = label_col,
        y        = "LOF Score",
        color    = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")
    plots <- c(plots, list(p6))
  }


  # Combine and annotate multi-panel figure
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = 2)

    methods_used    <- paste(detection_result$parameters$methods, collapse = ", ")
    total_anomalous <- length(detection_result$anomalous_labels)

    combined <- combined +
      patchwork::plot_annotation(
        title    = "Anomalous Label Detection Results",
        subtitle = sprintf(
          "Methods: %s  |  Anomalous: %d / %d labels",
          methods_used, total_anomalous, nrow(label_scores)
        ),
        theme = ggplot2::theme(
          plot.title    = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 10)
        )
      )

    return(combined)
  } else {
    return(plots)
  }
}


# Internal helpers ------------------------------------------------------------

#' Make Metric-specific Anomaly Factor
#'
#' @description
#' Internal helper used by \code{plot_label_anomaly()}. Creates a two-level
#' factor (\code{"Normal"} / \code{"Anomalous"}) from a per-metric flag
#' column in the plot data.
#'
#' @param data Data frame containing flag columns.
#' @param metric_name Name of the metric; the function looks for a column
#'   named \code{<metric_name>_flag}.
#'
#' @return Factor with levels \code{"Normal"} and \code{"Anomalous"}.
#' @keywords internal
#' @noRd
make_metric_anomaly_factor <- function(data, metric_name) {
  flag_col <- paste0(metric_name, "_flag")
  if (flag_col %in% names(data)) {
    values <- data[[flag_col]]
    values[is.na(values)] <- FALSE
    factor(values, levels = c(FALSE, TRUE), labels = c("Normal", "Anomalous"))
  } else {
    factor(rep("Normal", nrow(data)), levels = c("Normal", "Anomalous"))
  }
}


#' Build Anomaly Flag Heatmap
#'
#' @description
#' Internal helper used by \code{plot_label_anomaly()}. Builds the method-level
#' flag heatmap showing per-metric trajectory flags (top rows) and method-level
#' summary flags (bottom rows).
#'
#' @param plot_data Data frame with label scores and flag columns.
#' @param detection_result Output from \code{detect_anomalous_labels()}.
#' @param label_col Name of the label column.
#' @param highlight_color Color for anomalous tiles.
#'
#' @return A \code{ggplot} object, or \code{NULL} when no flag columns exist.
#' @keywords internal
#' @noRd
build_anomaly_heatmap <- function(plot_data, detection_result,
                                  label_col, highlight_color) {

  ensure_pkgs("tidyr")

  # Per-metric trajectory flag columns
  metric_flag_cols <- grep(
    "^mean_(rms_similarity|correlation|orthogonal_rms|parallel_rms|total_rms).*_flag$",
    names(plot_data), value = TRUE
  )

  # Method-level summary flag columns
  method_flag_cols <- character(0)
  for (mflag in c("trajectory_flag", "multivariate_flag", "lof_flag")) {
    if (mflag %in% names(detection_result$method_flags)) {
      plot_data[[mflag]] <- detection_result$method_flags[[mflag]]
      method_flag_cols   <- c(method_flag_cols, mflag)
    }
  }

  all_flag_cols <- c(metric_flag_cols, method_flag_cols)
  if (length(all_flag_cols) == 0) return(NULL)

  heatmap_data <- plot_data |>
    dplyr::select(dplyr::all_of(c(label_col, "label_numeric", all_flag_cols))) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(all_flag_cols),
      names_to  = "metric_name",
      values_to = "flagged"
    ) |>
    dplyr::mutate(
      metric_name = dplyr::case_when(
        metric_name == "trajectory_flag"   ~ "[method] trajectory",
        metric_name == "multivariate_flag" ~ "[method] multivariate",
        metric_name == "lof_flag"          ~ "[method] LOF",
        TRUE                               ~ metric_name
      ),
      metric_name = gsub("^mean_",  "", metric_name),
      metric_name = gsub("_flag$",  "", metric_name),
      metric_name = gsub("_scaled", " (scaled)", metric_name),
      # Method rows sort to the bottom; metric rows sort alphabetically on top
      row_group   = ifelse(
        grepl("^\\[method\\]", metric_name), "b_method", "a_metric"
      ),
      metric_name = factor(
        metric_name,
        levels = rev(unique(metric_name[order(row_group, metric_name)]))
      ),
      flagged = as.character(flagged)
    )

  ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(x = label_numeric, y = metric_name, fill = flagged)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values   = c("FALSE" = "gray92", "TRUE" = highlight_color),
      na.value = "gray92",
      labels   = c("FALSE" = "Normal", "TRUE" = "Anomalous"),
      name     = NULL
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::labs(
      title    = "Anomaly Flag Heatmap",
      subtitle = "Top rows: per-metric trajectory flags  |  Bottom rows: method-level summaries",
      x        = label_col,
      y        = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}
