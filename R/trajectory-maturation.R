# Trajectory Maturation Analysis
# Update date : Jun. 29, 2026

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
#' @param normalize_variability Character. How to normalize variability for
#'   cross-animal comparison: "none" (default), or "reference" (normalize to
#'   reference label). When "reference", variability is divided by the mean
#'   variability at the reference label, making values interpretable as
#'   "X times more variable than reference"
#' @param reference_label Character. Label to use as normalization reference.
#'   If NULL (default), uses the last label (typically adult). Only used when
#'   normalize_variability = "reference"
#' @param norm_epsilon Numeric. Small constant added to reference variability
#'   to avoid division by zero (default: 1e-6)
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
#' For data.frame: A data.frame with original data plus score columns:
#' \itemize{
#'   \item \code{variability_raw}: Original variability metric value (renamed from input)
#'   \item \code{variability_scaled}: Scaled variability (for ML)
#'   \item \code{variability_normalized}: Normalized to reference label (if normalize_variability = "reference").
#'     Values represent "X times the reference variability". For example, 2.0 means
#'     twice as variable as the reference label (typically adult)
#'   \item \code{maturation_score}: Computed maturation score (if requested)
#'   \item \code{stability_index}: Computed stability index (if requested)
#' }
#' Metadata attributes include: variability_metric, scale_method, similarity_metric,
#' score_type, invert_variability, normalize_variability, reference_label
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
#' # With reference-based normalization (for cross-animal ML)
#' sap <- trajectory_maturation(sap,
#'   segment_type = "motifs",
#'   normalize_variability = "reference",
#'   reference_label = "90"  # Adult as reference
#' )
#'
#' # Access results
#' scores <- sap$features$motif$maturation_scores
#' # variability_normalized = 2.0 means "2x adult variability"
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
    x, # x is data.frame
    similarity_metric = c("rms", "frechet", "dtw", "correlation"),
    variability_metric = c("dispersion", "orthogonal_rms", "parallel_rms"),
    score_type = c("maturation", "stability", "both"),
    invert_variability = TRUE,
    epsilon = 0.01,
    normalize_variability = c("none", "reference"),
    reference_label = NULL,
    norm_epsilon = 1e-6,
    scale_method = c("minmax", "zscore", "none"),
    verbose = TRUE,
    ...) {
  # Validate inputs
  if (!is.data.frame(x)) stop("Input must be a data.frame")

  similarity_metric <- match.arg(similarity_metric)
  variability_metric <- match.arg(variability_metric)
  scale_method <- match.arg(scale_method)
  normalize_variability <- match.arg(normalize_variability)

  if ("both" %in% score_type) {
    score_type <- c("maturation", "stability")
  } else {
    score_type <- match.arg(score_type, several.ok = TRUE)
  }

  # Check required columns
  required_cols <- c("label", "rendition")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
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
    message(sprintf("Scaling method     : %s", scale_method))
    if (normalize_variability == "reference") {
      ref_label_display <- if (is.null(reference_label)) "last label (auto)" else reference_label
      message(sprintf("Normalize to       : %s", ref_label_display))
    }
    message("")
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

  # Store raw and scaled variability for ML/analysis
  result$variability_raw <- variability
  result$variability_scaled <- variability_scaled

  # Normalize variability to reference label (for cross-animal comparison)
  if (normalize_variability == "reference") {
    # Determine reference label
    if (is.null(reference_label)) {
      # Use last label (typically adult)
      all_labels <- unique(x$label)
      reference_label <- sort_labels(all_labels)[length(all_labels)]
      if (verbose) {
        message(sprintf("Using last label as reference: %s", reference_label))
      }
    }
    
    # Check reference label exists
    if (!reference_label %in% x$label) {
      stop(sprintf(
        "Reference label '%s' not found in data. Available labels: %s",
        reference_label,
        paste(unique(x$label), collapse = ", ")
      ))
    }
    
    # Identify all variability columns to normalize
    # Look for columns with "variability_" prefix (from Sap merge)
    # Or raw metric names like "dispersion", "orthogonal_rms", "parallel_rms"
    variability_cols <- grep("^variability_", names(x), value = TRUE)
    
    # Also check for raw metric columns
    raw_metrics <- c("dispersion", "orthogonal_rms", "parallel_rms", "total_rms")
    raw_metrics_present <- raw_metrics[raw_metrics %in% names(x)]
    
    # Combine all variability columns to normalize
    all_var_cols <- unique(c(variability_cols, raw_metrics_present))
    
    if (verbose && length(all_var_cols) > 0) {
      message(sprintf("Normalizing %d variability metric(s):", length(all_var_cols)))
      for (col in all_var_cols) {
        message(sprintf("  - %s", col))
      }
    }
    
    # Normalize each variability metric
    ref_mask <- x$label == reference_label
    
    for (col in all_var_cols) {
      # Extract values
      values <- x[[col]]
      
      # Scale first (if not already scaled)
      values_scaled <- scale_variability(values, method = scale_method)
      
      # Compute reference mean
      ref_mean <- mean(values_scaled[ref_mask], na.rm = TRUE)
      
      # Normalize
      normalized_col_name <- paste0(col, "_normalized")
      result[[normalized_col_name]] <- values_scaled / (ref_mean + norm_epsilon)
      
      if (verbose) {
        message(sprintf(
          "  %s -> %s (ref: %.4f)",
          col, normalized_col_name, ref_mean
        ))
      }
    }
    
    # Also store the primary normalized variability
    # (the one selected for maturation scoring)
    result$variability_normalized <- variability_scaled / 
      (mean(variability_scaled[ref_mask], na.rm = TRUE) + norm_epsilon)
    
    if (verbose) {
      message(sprintf(
        "Primary variability_normalized range: [%.2f, %.2f]",
        min(result$variability_normalized, na.rm = TRUE),
        max(result$variability_normalized, na.rm = TRUE)
      ))
    }
  } else {
    # No normalization - just copy scaled values
    result$variability_normalized <- variability_scaled
  }

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
      message(sprintf(
        "Maturation score: mean = %.3f, sd = %.3f",
        mean(result$maturation_score, na.rm = TRUE),
        sd(result$maturation_score, na.rm = TRUE)
      ))
    }
  }

  if ("stability" %in% score_type) {
    # Stability uses raw scaled variability (not inverted)
    result$stability_index <- similarity / (variability_scaled + epsilon)
    if (verbose) {
      message(sprintf(
        "Stability index: mean = %.3f, sd = %.3f",
        mean(result$stability_index, na.rm = TRUE),
        sd(result$stability_index, na.rm = TRUE)
      ))
    }
  }

  # Store metadata
  attr(result, "similarity_metric") <- similarity_metric
  attr(result, "variability_metric") <- variability_metric
  attr(result, "score_type") <- score_type
  attr(result, "scale_method") <- scale_method
  attr(result, "invert_variability") <- invert_variability
  attr(result, "normalize_variability") <- normalize_variability
  if (normalize_variability == "reference") {
    attr(result, "reference_label") <- reference_label
  }

  if (verbose) message("Done.\n")

  return(result)
}


#' Clean Duplicate .source_row Columns
#'
#' @description
#' Internal helper to clean up duplicate .source_row columns that can occur
#' during merge operations (e.g., .source_row.x, .source_row.y).
#'
#' @param df Data frame to clean
#'
#' @return Data frame with single .source_row column
#'
#' @keywords internal
#' @noRd
clean_source_row_columns <- function(df) {
  # If we already have a clean .source_row, we're good
  if (".source_row" %in% names(df)) {
    # Check if there are also .source_row.x or .source_row.y
    dup_cols <- grep("^\\.source_row\\.[xy]", names(df), value = TRUE)
    if (length(dup_cols) > 0) {
      # Remove the duplicates, keep the original
      df <- df[, !names(df) %in% dup_cols, drop = FALSE]
    }
    return(df)
  }
  
  # Look for .source_row.x, .source_row.y, etc.
  dup_cols <- grep("^\\.source_row\\.", names(df), value = TRUE)
  
  if (length(dup_cols) > 0) {
    # Create .source_row by taking first non-NA value for each row
    df$.source_row <- apply(df[, dup_cols, drop = FALSE], 1, function(row) {
      non_na <- row[!is.na(row)]
      if (length(non_na) > 0) non_na[1] else NA_integer_
    })
    
    # Drop the duplicate columns
    df <- df[, !names(df) %in% dup_cols, drop = FALSE]
  }
  
  return(df)
}


#' Deduplicate ref_* Columns After Merge
#'
#' @description
#' Internal helper to collapse any \code{ref_*.x} / \code{ref_*.y} pairs that
#' arise when \code{ref_day} or \code{ref_scale_*} columns are present in both
#' sides of a merge.  Values should be identical (same animal, same run), so
#' the \code{.x} value is kept and the \code{.y} duplicate is dropped.
#'
#' @param df Data frame to clean
#'
#' @return Data frame with one column per \code{ref_*} variable
#'
#' @keywords internal
#' @noRd
clean_ref_columns <- function(df) {
  # Find columns ending in .x whose base name starts with ref_
  x_cols <- grep("\\.x$", names(df), value = TRUE)
  ref_x_cols <- x_cols[grepl("^ref_", sub("\\.x$", "", x_cols))]

  for (xcol in ref_x_cols) {
    base <- sub("\\.x$", "", xcol)
    ycol <- paste0(base, ".y")
    # Coalesce: prefer .x value, fall back to .y when .x is NA
    if (ycol %in% names(df)) {
      df[[base]] <- ifelse(!is.na(df[[xcol]]), df[[xcol]], df[[ycol]])
      df[[xcol]] <- NULL
      df[[ycol]] <- NULL
    }
  }
  df
}


#' @rdname trajectory_maturation
#' @export
trajectory_maturation.Sap <- function(
    x, # x is SAP object
    segment_type = c("motifs", "syllables", "bouts", "segments"),
    similarity_metric = c("rms", "frechet", "dtw", "correlation"),
    variability_metric = c("dispersion", "orthogonal_rms", "parallel_rms"),
    score_type = c("maturation", "stability", "both"),
    invert_variability = TRUE,
    epsilon = 0.01,
    normalize_variability = c("none", "reference"),
    reference_label = NULL,
    norm_epsilon = 1e-6,
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
      paste(
        "No trajectory variability results found for %s.",
        "Run trajectory_dispersion() or trajectory_path_deviation() first."
      ),
      segment_type
    ))
  }

  # Merge similarity and variability data
  sim_sim <- sim_result$similarity
  var_sim <- data.frame(
    label = character(), rendition = character(),
    stringsAsFactors = FALSE
  )

  # Extract dispersion metrics if available
  if (!is.null(var_result) && var_result$type == "dispersion") {
    disp_cols <- c("label", "rendition", "dispersion")
    if (".source_row" %in% names(var_result$dispersion)) {
      disp_cols <- c(disp_cols, ".source_row")
    }
    # Pull ref_scale_* covariates; exclude ref_day — it already comes from sim_sim
    var_ref_scale_cols <- grep("^ref_scale_", names(var_result$dispersion), value = TRUE)
    disp_cols <- c(disp_cols, var_ref_scale_cols)
    disp_data <- var_result$dispersion[, disp_cols, drop = FALSE]
    names(disp_data)[names(disp_data) == "dispersion"] <- "variability_dispersion"

    if (nrow(var_sim) == 0) {
      var_sim <- disp_data
    } else {
      var_sim <- merge(var_sim, disp_data, by = c("label", "rendition"), all = TRUE)
    }
  }

  # Extract path deviation metrics if available
  if (!is.null(var_result_path_dev) && var_result_path_dev$type == "path_deviation") {
    width_cols <- c("label", "rendition", "orthogonal_rms", "parallel_rms")
    if (".source_row" %in% names(var_result_path_dev$width)) {
      width_cols <- c(width_cols, ".source_row")
    }
    # Pull ref_scale_* covariates; exclude ref_day — it already comes from sim_sim
    var_ref_scale_cols <- grep("^ref_scale_", names(var_result_path_dev$width), value = TRUE)
    width_cols <- c(width_cols, var_ref_scale_cols)
    width_data <- var_result_path_dev$width[, width_cols, drop = FALSE]
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

  # Determine merge keys — use .source_row only if present in both
  merge_keys <- c("label", "rendition")
  if (".source_row" %in% names(sim_sim) && ".source_row" %in% names(var_sim)) {
    merge_keys <- c(merge_keys, ".source_row")
  }

  merged_data <- merge(sim_sim, var_sim,
    by = merge_keys,
    suffixes = c("_sim", "_var")
  )

  if (nrow(merged_data) == 0) {
    stop("No matching renditions found between similarity and variability results")
  }

  # Clean up duplicate columns that can arise from merge operations
  merged_data <- clean_source_row_columns(merged_data)  # .source_row.x/.y
  merged_data <- clean_ref_columns(merged_data)          # ref_day.x/.y, ref_scale_*.x/.y

  # Compute scores
  scores <- trajectory_maturation.default(
    x = merged_data,
    similarity_metric = similarity_metric,
    variability_metric = variability_metric,
    score_type = score_type,
    invert_variability = invert_variability,
    epsilon = epsilon,
    normalize_variability = normalize_variability,
    reference_label = reference_label,
    norm_epsilon = norm_epsilon,
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
    x, # x is data.frame
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
      names(x),
      value = TRUE
    )
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
        brackets_data <- brackets(
          x[[score_col]], tests[[score_col]]$posthoc,
          labs_order, max_annotations
        )
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
    x, # x is SAP object
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

# Trajectory Trial Score Plot
# Update date : Jun. 26, 2026

#' Plot Individual Trial Scores
#'
#' @description
#' Creates a dot plot of per-trial similarity, maturation, or variability scores,
#' ordered by \code{.source_row} within each label/day. Helps visualise the
#' sequential progression of scores across individual renditions within and
#' across days. An optional running-mean trace (default window: 50 trials)
#' highlights the local trend.
#'
#' @param x An object to plot: a data frame with per-trial scores, or a SAP
#'   object with pre-computed trajectory results.
#' @param score_col Character. Column name of the score to plot
#'   (default: "rms_similarity").
#' @param labels Character vector of labels (days) to include. Default NULL
#'   shows the first and last label when labels are sorted numerically.
#' @param running_mean Logical. If TRUE (default), overlay a running-mean trace.
#' @param window_size Integer. Number of trials for the running-mean window
#'   (default: 50).
#' @param palette Character. RColorBrewer palette name (default: "Set1").
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables',
#'   'bouts', 'segments').
#' @param data_type For SAP objects: Which results to plot ('similarity',
#'   'maturation', 'dispersion', or 'path_deviation').
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' **Data source for SAP objects:**
#' \itemize{
#'   \item \code{data_type = "similarity"}: Uses \code{x$features[[feature_type]]$trajectory_similarity$similarity}.
#'     Common score columns: "rms_similarity", "correlation", "frechet_similarity",
#'     "dtw_similarity".
#'   \item \code{data_type = "maturation"}: Uses \code{x$features[[feature_type]]$maturation_scores}.
#'     Common score columns: "maturation_score", "stability_index".
#'   \item \code{data_type = "dispersion"}: Uses \code{x$features[[feature_type]]$trajectory_dispersion$dispersion}.
#'     Score column: "dispersion".
#'   \item \code{data_type = "path_deviation"}: Uses \code{x$features[[feature_type]]$trajectory_path_deviation$width}.
#'     Score columns: "orthogonal_rms", "parallel_rms".
#' }
#'
#' **Label selection:**
#' When \code{labels = NULL}, the function selects the first and last labels
#' from the sorted unique labels in the data (numeric-aware sorting). This
#' typically corresponds to the earliest and latest days in a developmental
#' series. Pass an explicit vector to customise.
#'
#' **Running mean:**
#' The running mean is computed with \code{stats::filter()} using a centred
#' window of \code{window_size} trials. It is only drawn when
#' \code{running_mean = TRUE} and there are at least \code{window_size}
#' trials in a label.
#'
#' @return
#' A ggplot2 object, printed as a side-effect and returned invisibly.
#'
#' @examples
#' \dontrun{
#' # Plot RMS similarity for first and last labels
#' plot_trajectory_trials(sap,
#'   data_type = "similarity",
#'   score_col = "rms_similarity"
#' )
#'
#' # Plot maturation scores for specific days without running mean
#' plot_trajectory_trials(sap,
#'   data_type = "maturation",
#'   score_col = "maturation_score",
#'   labels = c("60", "80", "100"),
#'   running_mean = FALSE
#' )
#'
#' # Plot dispersion variability
#' plot_trajectory_trials(sap,
#'   data_type = "dispersion",
#'   score_col = "dispersion"
#' )
#'
#' # Plot orthogonal RMS from path deviation
#' plot_trajectory_trials(sap,
#'   data_type = "path_deviation",
#'   score_col = "orthogonal_rms"
#' )
#'
#' # Plot parallel RMS from path deviation
#' plot_trajectory_trials(sap,
#'   data_type = "path_deviation",
#'   score_col = "parallel_rms",
#'   labels = c("60", "90")
#' )
#' }
#'
#' @seealso \code{\link{trajectory_similarity}} for similarity computation,
#'   \code{\link{trajectory_maturation}} for maturation score computation
#'
#' @rdname plot_trajectory_trials
#' @export
plot_trajectory_trials <- function(x, ...) {
  UseMethod("plot_trajectory_trials")
}


#' @rdname plot_trajectory_trials
#' @export
plot_trajectory_trials.default <- function(
    x, # x is data.frame with label, .source_row, score_col
    score_col = "rms_similarity",
    labels = NULL,
    running_mean = TRUE,
    window_size = 50,
    palette = "Set1",
    ...) {
  if (!is.data.frame(x)) stop("Input must be a data.frame")

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  # Validate score column
  if (!score_col %in% names(x)) {
    stop(sprintf(
      "Score column '%s' not found. Available columns: %s",
      score_col, paste(names(x), collapse = ", ")
    ))
  }

  # Sort labels numerically if possible
  labs <- sort_labels(unique(as.character(x$label)))

  # Default: first and last label
  if (is.null(labels)) {
    labels <- labs[c(1, length(labs))]
  } else {
    invalid <- setdiff(labels, labs)
    if (length(invalid) > 0) {
      stop(sprintf("Labels not found: %s", paste(invalid, collapse = ", ")))
    }
  }

  # Filter to requested labels and sort by .source_row within each label
  plot_data <- x[as.character(x$label) %in% labels, , drop = FALSE]
  plot_data <- plot_data[order(plot_data$label, plot_data$.source_row), , drop = FALSE]

  # Add trial index within each label
  plot_data$trial_idx <- stats::ave(seq_len(nrow(plot_data)),
    plot_data$label,
    FUN = seq_along
  )

  # Build colour palette and set factor order for facets
  # Preserve the input order if labels were specified, otherwise use numeric sort
  if (!is.null(labels)) {
    labs_order <- labels
  } else {
    labs_order <- sort_labels(unique(as.character(plot_data$label)))
  }
  
  # Set label as factor with correct order for faceting
  plot_data$label <- factor(plot_data$label, levels = labs_order)
  
  pal_map <- make_pal(labs_order, palette)

  # Running mean
  if (running_mean && window_size > 1) {
    running_list <- by(plot_data, plot_data$label, function(sub) {
      vals <- sub[[score_col]]
      if (length(vals) >= window_size) {
        rm <- as.numeric(stats::filter(vals, rep(1 / window_size, window_size),
          sides = 2
        ))
        data.frame(
          label = sub$label[1],
          trial_idx = sub$trial_idx,
          running_mean = rm,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    })
    running_df <- do.call(rbind, running_list)
  }

  # Build plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$trial_idx,
    y = .data[[score_col]],
    colour = .data$label
  )) +
    ggplot2::geom_point(size = 1, alpha = 0.5) +
    ggplot2::facet_wrap(~label, scales = "free_x", ncol = 1) +
    ggplot2::labs(
      title = sprintf("Trial-level %s", score_col),
      x = "Trial index (ordered by time sequence)",
      y = score_col,
      colour = "Label"
    ) +
    ggplot2::scale_colour_manual(values = pal_map) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "right"
    )

  # Add running mean trace
  if (running_mean && window_size > 1 && exists("running_df") &&
    !is.null(running_df) && nrow(running_df) > 0) {
    p <- p + ggplot2::geom_line(
      data = running_df,
      ggplot2::aes(x = .data$trial_idx, y = .data$running_mean),
      colour = "black", linewidth = 0.8, na.rm = TRUE
    )
  }

  print(p)
  invisible(p)
}


#' @rdname plot_trajectory_trials
#' @export
plot_trajectory_trials.Sap <- function(
    x, # x is SAP object
    segment_type = c("motifs", "syllables", "bouts", "segments"),
    data_type = c("similarity", "maturation", "dispersion", "path_deviation"),
    score_col = NULL,
    labels = NULL,
    running_mean = TRUE,
    window_size = 50,
    palette = "Set1",
    ...) {
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  segment_type <- match.arg(segment_type)
  data_type <- match.arg(data_type)
  feature_type <- sub("s$", "", segment_type)

  # Extract the appropriate results
  if (data_type == "similarity") {
    sim_result <- x$features[[feature_type]][["trajectory_similarity"]]
    if (is.null(sim_result)) {
      stop(sprintf(
        "No trajectory similarity results found for %s. Run trajectory_similarity() first.",
        segment_type
      ))
    }
    plot_data <- sim_result$similarity

    # Default score column for similarity
    if (is.null(score_col)) {
      score_col <- if ("rms_similarity" %in% names(plot_data)) {
        "rms_similarity"
      } else if ("correlation" %in% names(plot_data)) {
        "correlation"
      } else {
        names(plot_data)[!names(plot_data) %in% c("label", "rendition", ".source_row")][1]
      }
    }
  } else if (data_type == "maturation") {
    scores <- x$features[[feature_type]][["maturation_scores"]]
    if (is.null(scores)) {
      stop(sprintf(
        "No maturation scores found for %s. Run trajectory_maturation() first.",
        segment_type
      ))
    }
    plot_data <- scores

    # Default score column for maturation
    if (is.null(score_col)) {
      score_col <- if ("maturation_score" %in% names(plot_data)) {
        "maturation_score"
      } else if ("stability_index" %in% names(plot_data)) {
        "stability_index"
      } else {
        names(plot_data)[!names(plot_data) %in% c("label", "rendition", ".source_row")][1]
      }
    }
  } else if (data_type == "dispersion") {
    disp_result <- x$features[[feature_type]][["trajectory_dispersion"]]
    if (is.null(disp_result)) {
      stop(sprintf(
        "No trajectory dispersion results found for %s. Run trajectory_dispersion() first.",
        segment_type
      ))
    }
    plot_data <- disp_result$dispersion

    # Default score column for dispersion
    if (is.null(score_col)) {
      score_col <- "dispersion"
    }
  } else if (data_type == "path_deviation") {
    path_result <- x$features[[feature_type]][["trajectory_path_deviation"]]
    if (is.null(path_result)) {
      stop(sprintf(
        "No trajectory path deviation results found for %s. Run trajectory_path_deviation() first.",
        segment_type
      ))
    }
    plot_data <- path_result$width

    # Default score column for path_deviation
    if (is.null(score_col)) {
      score_col <- "orthogonal_rms"
    }
  }

  # Check that score_col actually exists
  if (!score_col %in% names(plot_data)) {
    stop(sprintf(
      "Score column '%s' not found in %s results. Available: %s",
      score_col, data_type, paste(names(plot_data), collapse = ", ")
    ))
  }

  plot_trajectory_trials.default(
    x = plot_data,
    score_col = score_col,
    labels = labels,
    running_mean = running_mean,
    window_size = window_size,
    palette = palette,
    ...
  )
}


