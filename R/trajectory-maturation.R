# Trajectory Maturation Analysis
# Update date : Jun. 25, 2026

#' Compute Trajectory Maturation Scores
#'
#' @description
#' Calculate composite developmental maturity scores by combining trajectory
#' similarity (similarity to reference) and within-rendition variability
#' (consistency). Higher scores indicate more mature, crystallized trajectories.
#'
#' @param x An object to analyze: a data.frame with similarity and variability
#'   results, or a SAP object
#' @param similarity_metric Character. Which similarity metric to use:
#'   "rms" (default), "frechet", "dtw", or "correlation"
#' @param variability_metric Character. Which variability metric to use:
#'   "dispersion" (default), "orthogonal_rms", or "parallel_rms"
#' @param score_type Character vector. Which scores to compute:
#'   "maturation" (default), "stability", or "both"
#' @param invert_variability Logical. If TRUE (default), higher variability
#'   reduces the score
#' @param epsilon Numeric. Small constant to avoid division by zero in
#'   stability_index (default: 0.01)
#' @param scale_method Character. How to scale variability: "minmax" (default),
#'   "zscore", or "none"
#' @param verbose Logical. Print progress messages (default: TRUE)
#' @param segment_type For SAP objects: Type of segments ("motifs", "syllables",
#'   "bouts", or "segments")
#' @param ... Additional arguments
#'
#' @details
#' **Score Types:**
#' \itemize{
#'   \item Maturation score: similarity × (1 - scaled_variability). Combines
#'     convergence to reference with consistency into a unified developmental
#'     score
#'   \item Stability index: similarity / (scaled_variability + epsilon). Ratio
#'     form that shows relative stability independent of scale
#' }
#'
#' **Prerequisites:**
#' Must first run \code{\link{trajectory_similarity}} and either
#' \code{\link{trajectory_dispersion}} or \code{\link{trajectory_path_deviation}}.
#'
#' **Variability Metrics:**
#' \itemize{
#'   \item \code{dispersion}: Mean distance to centroid trajectory (from
#'     trajectory_dispersion)
#'   \item \code{orthogonal_rms}: Perpendicular deviation from mean trajectory
#'     (from trajectory_path_deviation)
#'   \item \code{parallel_rms}: Along-trajectory timing/phase deviation (from
#'     trajectory_path_deviation)
#' }
#'
#' @return
#' For data.frame: A data.frame with original data plus score columns.
#'
#' For SAP objects: Updated object with scores stored in
#'   \code{sap$features[[feature_type]]$maturation_scores}
#'
#' @examples
#' \dontrun{
#' # After running trajectory_similarity and trajectory_dispersion
#' sap <- trajectory_maturation(sap,
#'   segment_type = "motifs",
#'   similarity_metric = "rms",
#'   score_type = "both"
#' )
#'
#' # Access results
#' scores <- sap$features$motif$maturation_scores
#' }
#'
#' @seealso \code{\link{plot_trajectory_maturation}} for visualization,
#'   \code{\link{trajectory_similarity}} for similarity computation,
#'   \code{\link{trajectory_dispersion}} for dispersion-based variability
#'
#' @rdname trajectory_maturation
#' @export
trajectory_maturation <- function(x, ...) {
  UseMethod("trajectory_maturation")
}


#' @rdname trajectory_maturation
#' @export
trajectory_maturation.default <- function(
    x,                      # x is data.frame
    similarity_metric = c("rms", "frechet", "dtw", "correlation"),
    variability_metric = c("dispersion", "orthogonal_rms", "parallel_rms"),
    score_type = c("maturation", "stability", "both"),
    invert_variability = TRUE,
    epsilon = 0.01,
    scale_method = c("minmax", "zscore", "none"),
    verbose = TRUE,
    ...) {

  # Validate inputs
  if (!is.data.frame(x)) stop("Input must be a data.frame")

  similarity_metric <- match.arg(similarity_metric)
  variability_metric <- match.arg(variability_metric)
  scale_method <- match.arg(scale_method)

  if ("both" %in% score_type) {
    score_type <- c("maturation", "stability")
  } else {
    score_type <- match.arg(score_type, several.ok = TRUE)
  }

  # Check required columns
  required_cols <- c("label", "rendition")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Map metric names to column names
  sim_col <- paste0(similarity_metric, "_similarity")
  if (similarity_metric == "correlation") sim_col <- "correlation"

  if (!sim_col %in% names(x)) {
    stop(sprintf(
      "Similarity metric '%s' not found. Available columns: %s",
      sim_col,
      paste(names(x), collapse = ", ")
    ))
  }

  # Try both with and without variability_ prefix
  var_col <- variability_metric
  var_col_prefixed <- paste0("variability_", variability_metric)

  if (var_col_prefixed %in% names(x)) {
    var_col <- var_col_prefixed
  } else if (!var_col %in% names(x)) {
    stop(sprintf(
      "Variability metric '%s' not found. Available columns: %s",
      variability_metric,
      paste(names(x), collapse = ", ")
    ))
  }

  if (verbose) {
    message("\n=== Computing Trajectory Maturation Scores ===")
    message(sprintf("Similarity metric  : %s", similarity_metric))
    message(sprintf("Variability metric : %s", variability_metric))
    message(sprintf("Score types        : %s", paste(score_type, collapse = ", ")))
    message(sprintf("Scaling method     : %s\n", scale_method))
  }

  # Extract metrics
  similarity <- x[[sim_col]]
  variability <- x[[var_col]]

  # Scale variability
  variability_scaled <- scale_variability(
    variability = variability,
    method = scale_method
  )

  result <- x

  # Compute requested scores
  if ("maturation" %in% score_type) {
    # Maturation uses inverted variability (high variability = bad)
    variability_inverted <- if (invert_variability) {
      1 - variability_scaled
    } else {
      variability_scaled
    }
    result$maturation_score <- similarity * variability_inverted
    if (verbose) {
      message(sprintf("Maturation score: mean = %.3f, sd = %.3f",
                      mean(result$maturation_score, na.rm = TRUE),
                      sd(result$maturation_score, na.rm = TRUE)))
    }
  }

  if ("stability" %in% score_type) {
    # Stability uses raw scaled variability (not inverted)
    result$stability_index <- similarity / (variability_scaled + epsilon)
    if (verbose) {
      message(sprintf("Stability index: mean = %.3f, sd = %.3f",
                      mean(result$stability_index, na.rm = TRUE),
                      sd(result$stability_index, na.rm = TRUE)))
    }
  }

  # Store metadata
  attr(result, "similarity_metric") <- similarity_metric
  attr(result, "variability_metric") <- variability_metric
  attr(result, "score_type") <- score_type
  attr(result, "scale_method") <- scale_method

  if (verbose) message("Done.\n")

  return(result)
}


#' @rdname trajectory_maturation
#' @export
trajectory_maturation.Sap <- function(
    x,                      # x is SAP object
    segment_type = c("motifs", "syllables", "bouts", "segments"),
    similarity_metric = c("rms", "frechet", "dtw", "correlation"),
    variability_metric = c("dispersion", "orthogonal_rms", "parallel_rms"),
    score_type = c("maturation", "stability", "both"),
    invert_variability = TRUE,
    epsilon = 0.01,
    scale_method = c("minmax", "zscore", "none"),
    verbose = TRUE,
    ...) {

  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  # Get similarity results
  sim_result <- x$features[[feature_type]][["trajectory_similarity"]]
  if (is.null(sim_result)) {
    stop(sprintf(
      "No trajectory similarity results found for %s. Run trajectory_similarity() first.",
      segment_type
    ))
  }

  # Get variability results
  var_result <- x$features[[feature_type]][["trajectory_dispersion"]]
  var_result_path_dev <- x$features[[feature_type]][["trajectory_path_deviation"]]

  if (is.null(var_result) && is.null(var_result_path_dev)) {
    stop(sprintf(
      paste("No trajectory variability results found for %s.",
            "Run trajectory_dispersion() or trajectory_path_deviation() first."),
      segment_type
    ))
  }

  # Merge similarity and variability data
  sim_sim <- sim_result$similarity
  var_sim <- data.frame(label = character(), rendition = character(),
                        stringsAsFactors = FALSE)

  # Extract dispersion metrics if available
  if (!is.null(var_result) && var_result$type == "dispersion") {
    disp_data <- var_result$dispersion[, c("label", "rendition", "dispersion")]
    names(disp_data)[names(disp_data) == "dispersion"] <- "variability_dispersion"

    if (nrow(var_sim) == 0) {
      var_sim <- disp_data
    } else {
      var_sim <- merge(var_sim, disp_data, by = c("label", "rendition"), all = TRUE)
    }
  }

  # Extract path deviation metrics if available
  if (!is.null(var_result_path_dev) && var_result_path_dev$type == "path_deviation") {
    width_data <- var_result_path_dev$width[, c("label", "rendition", "orthogonal_rms", "parallel_rms")]
    names(width_data)[names(width_data) == "orthogonal_rms"] <- "variability_orthogonal_rms"
    names(width_data)[names(width_data) == "parallel_rms"] <- "variability_parallel_rms"

    if (nrow(var_sim) == 0) {
      var_sim <- width_data
    } else {
      var_sim <- merge(var_sim, width_data, by = c("label", "rendition"), all = TRUE)
    }
  }

  if (nrow(var_sim) == 0) {
    stop("Could not extract variability metrics from results")
  }

  merged_data <- merge(sim_sim, var_sim,
                       by = c("label", "rendition"),
                       suffixes = c("_sim", "_var"))

  if (nrow(merged_data) == 0) {
    stop("No matching renditions found between similarity and variability results")
  }

  # Compute scores
  scores <- trajectory_maturation.default(
    x = merged_data,
    similarity_metric = similarity_metric,
    variability_metric = variability_metric,
    score_type = score_type,
    invert_variability = invert_variability,
    epsilon = epsilon,
    scale_method = scale_method,
    verbose = verbose,
    ...
  )

  # Store in Sap object
  x$features[[feature_type]][["maturation_scores"]] <- scores
  x$misc$last_modified <- Sys.time()

  invisible(x)
}


#' Scale Variability Values
#'
#' @description
#' Internal helper used by \code{trajectory_maturation.default()} to scale
#' variability values using minmax, zscore, or no scaling.
#'
#' @param variability Numeric vector of variability values
#' @param method Character. Scaling method: "minmax", "zscore", or "none"
#'
#' @keywords internal
#' @noRd
scale_variability <- function(variability, method = c("minmax", "zscore", "none")) {
  method <- match.arg(method)

  switch(method,
    minmax = {
      min_val <- min(variability, na.rm = TRUE)
      max_val <- max(variability, na.rm = TRUE)
      if (max_val - min_val > 0) {
        (variability - min_val) / (max_val - min_val)
      } else {
        rep(0, length(variability))
      }
    },
    zscore = {
      (variability - mean(variability, na.rm = TRUE)) /
        sd(variability, na.rm = TRUE)
    },
    none = variability
  )
}


# Plotting Functions ------------------------------------------------------

#' Plot Trajectory Maturation Scores
#'
#' @description
#' Visualize developmental maturation scores across labels or age groups.
#' Uses ASAP's standard trajectory plotting style (violin + boxplot panels).
#'
#' @param x An object to plot: a data.frame with maturation scores or a SAP object
#' @param score_cols Character vector. Which score columns to plot
#'   (default: auto-detect from available scores)
#' @param palette Character. RColorBrewer palette name (default: "Set1")
#' @param max_annotations Numeric. Maximum number of significance brackets
#'   to show (default: 10)
#' @param segment_type For SAP objects: Type of segments ("motifs", "syllables",
#'   "bouts", or "segments")
#' @param ... Additional arguments
#'
#' @details
#' This function uses ASAP's existing plotting infrastructure from trajectory.R:
#' \itemize{
#'   \item \code{sort_labels()}: Intelligent numeric-aware label sorting
#'   \item \code{make_pal()}: RColorBrewer palette generation
#'   \item \code{panel()}: Violin + boxplot panel creation
#'   \item \code{trend_panel()}: Line plot for >6 groups
#'   \item \code{fmt_p()}: P-value formatting
#'   \item \code{brackets()}: Significance bracket computation
#'   \item \code{add_brackets()}: Significance annotation overlay
#' }
#'
#' The resulting plots match the style of other ASAP trajectory visualizations.
#'
#' @return A ggplot2 or patchwork object
#'
#' @examples
#' \dontrun{
#' # Plot all computed scores
#' plot_trajectory_maturation(sap, segment_type = "motifs")
#'
#' # Plot specific score
#' plot_trajectory_maturation(scores_df, score_cols = "maturation_score")
#' }
#'
#' @seealso \code{\link{trajectory_maturation}} for score computation
#'
#' @rdname plot_trajectory_maturation
#' @export
plot_trajectory_maturation <- function(x, ...) {
  UseMethod("plot_trajectory_maturation")
}


#' @rdname plot_trajectory_maturation
#' @export
plot_trajectory_maturation.default <- function(
    x,                      # x is data.frame
    score_cols = NULL,
    palette = "Set1",
    max_annotations = 10,
    ...) {

  if (!is.data.frame(x)) stop("Input must be a data.frame")

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'patchwork' are required for plotting")
  }

  # Auto-detect score columns if not specified
  if (is.null(score_cols)) {
    score_cols <- grep("^(maturation_score|stability_index)$",
                       names(x), value = TRUE)
  }

  if (length(score_cols) == 0) {
    stop("No maturation score columns found. Run trajectory_maturation() first.")
  }

  # Validate columns exist
  missing_cols <- setdiff(score_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf("Score columns not found: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!"label" %in% names(x)) {
    stop("'label' column required for plotting")
  }

  # Sort labels using ASAP's sort_labels function
  labs_order <- sort_labels(unique(as.character(x$label)))
  many_labs <- length(labs_order) > 6

  # Create color palette using ASAP's make_pal function
  pal_map <- make_pal(labs_order, palette)

  # Compute statistical tests
  tests <- NULL
  if (length(labs_order) > 1) {
    tests <- lapply(score_cols, function(score_col) {
      kw <- stats::kruskal.test(x[[score_col]] ~ x$label)

      posthoc <- NULL
      if (kw$p.value < 0.05 && !many_labs) {
        pw <- stats::pairwise.wilcox.test(
          x[[score_col]],
          x$label,
          p.adjust.method = "bonferroni"
        )
        posthoc <- pw$p.value
      }

      list(kruskal = list(p.value = kw$p.value), posthoc = posthoc)
    })
    names(tests) <- score_cols
  }

  # Create panels for each score
  panels <- lapply(score_cols, function(score_col) {
    y_title <- format_score_title(score_col)
    kw_p <- if (!is.null(tests)) fmt_p(tests[[score_col]]$kruskal$p.value) else NULL

    if (many_labs) {
      # Use trend panel for many groups
      p <- trend_panel(x, "label", score_col, y_title, kw_p, labs_order)
    } else {
      # Use standard panel (violin + boxplot)
      p <- panel(x, "label", score_col, y_title, kw_p, pal_map, labs_order)

      # Add significance brackets
      if (!is.null(tests) && !is.null(tests[[score_col]]$posthoc)) {
        brackets_data <- brackets(x[[score_col]], tests[[score_col]]$posthoc,
                                  labs_order, max_annotations)
        p <- add_brackets(p, brackets_data)
      }
    }

    p
  })

  # Combine panels
  if (length(panels) == 1) {
    combined <- panels[[1]]
  } else {
    combined <- Reduce("+", panels) +
      patchwork::plot_layout(ncol = min(2, length(panels)))
  }

  # Add overall title with metadata
  sim_metric <- attr(x, "similarity_metric")
  var_metric <- attr(x, "variability_metric")

  title <- "Trajectory Maturation Scores"
  subtitle <- if (!is.null(sim_metric) && !is.null(var_metric)) {
    sprintf("Similarity: %s | Variability: %s", sim_metric, var_metric)
  } else {
    NULL
  }

  combined <- combined +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      )
    )

  combined
}


#' @rdname plot_trajectory_maturation
#' @export
plot_trajectory_maturation.Sap <- function(
    x,                      # x is SAP object
    segment_type = c("motifs", "syllables", "bouts", "segments"),
    score_cols = NULL,
    palette = "Set1",
    max_annotations = 10,
    ...) {

  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)

  scores <- x$features[[feature_type]][["maturation_scores"]]
  if (is.null(scores)) {
    stop(sprintf(
      "No maturation scores found for %s. Run trajectory_maturation() first.",
      segment_type
    ))
  }

  plot_trajectory_maturation.default(
    x = scores,
    score_cols = score_cols,
    palette = palette,
    max_annotations = max_annotations,
    ...
  )
}


#' Format Score Column Names as Titles
#'
#' @description
#' Internal helper used by \code{plot_trajectory_maturation.default()} to
#' convert score column names to readable titles.
#'
#' @param col_name Character. Column name to format
#'
#' @keywords internal
#' @noRd
format_score_title <- function(col_name) {
  title <- gsub("_", " ", col_name)
  tools::toTitleCase(title)
}
