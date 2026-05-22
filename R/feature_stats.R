# Feature Statistics Summary and Analysis
# Update date: May 22, 2026

#' Feature Statistics Summary and Analysis
#'
#' @description
#' A generic function to calculate, print, and plot feature statistics.
#' The default method accepts a pre-computed segment statistics data frame,
#' while the \code{Sap} method extracts feature statistics from a \code{Sap} object
#' and forwards them to the default method.
#'
#' @param x An object containing statistics: either a data frame (e.g.
#'   \code{sap$features$motif$stats_fund_freq}) or a \code{Sap} object.
#' @param feature_name For default method: Character string used as the y-axis label
#'   and plot title (default: \code{"Feature"}).
#' @param labels Optional character vector of labels to filter the data.
#' @param run_anova Logical. Whether to perform a one-way ANOVA (default: \code{TRUE}).
#' @param plot Logical. Whether to generate the ggplot visualization (default: \code{TRUE}).
#' @param plot_type Character. Style of visual layer: \code{"both"} (default),
#'   \code{"boxplot"}, or \code{"violin"}.
#' @param palette Character. RColorBrewer palette name for fills/colors (default: \code{"Set1"}).
#' @param jitter Logical. Overlay individual rendition data points (default: \code{TRUE}).
#' @param point_alpha Numeric \code{[0, 1]}. Transparency of jittered points (default: \code{0.35}).
#' @param facet_segments Logical. Facet the plot by \code{segment_id} (default: \code{FALSE}).
#' @param ncol Integer. Number of columns for facets if \code{facet_segments = TRUE}.
#' @param feature For SAP objects: Feature to analyze, either \code{"fund_freq"} or
#'   \code{"wiener_entropy"}.
#' @param segment_type For SAP objects: Type of segments to extract, either \code{"motifs"},
#'   \code{"syllables"}, or \code{"segments"}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A list (invisibly) with components:
#' \describe{
#'   \item{\code{summary}}{A tibble of per-label summary statistics.}
#'   \item{\code{anova}}{ANOVA result tibble, or \code{NULL} if not run.}
#'   \item{\code{plot}}{The ggplot object, or \code{NULL} if \code{plot = FALSE}.}
#' }
#'
#' @seealso \code{\link{anova_analysis}}, \code{\link{refine_FF}},
#'          \code{\link{refine_sh}}
#'
#' @export
feature_stats <- function(x, ...) {
  UseMethod("feature_stats")
}

#' @rdname feature_stats
#' @export
feature_stats.default <- function(x,
                                  feature_name   = "Feature",
                                  labels         = NULL,
                                  run_anova      = TRUE,
                                  plot           = TRUE,
                                  plot_type      = c("both", "boxplot", "violin"),
                                  palette        = "Set1",
                                  jitter         = TRUE,
                                  point_alpha    = 0.35,
                                  facet_segments = FALSE,
                                  ncol           = NULL,
                                  ...) {

  stats_df <- x

  # Validate inputs
  if (!is.data.frame(stats_df)) {
    stop("'x' must be a data frame (e.g. sap$features$motif$stats_fund_freq).")
  }

  required_cols <- c("label", "segment_id", "mean", "sd", "n_samples")
  missing_cols  <- setdiff(required_cols, names(stats_df))
  if (length(missing_cols) > 0) {
    stop("'x' is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  plot_type <- match.arg(plot_type)

  ensure_pkgs("ggplot2", "dplyr")

  # Label filtering
  if (!is.null(labels)) {
    available <- unique(stats_df$label)
    missing_l <- setdiff(labels, available)
    if (length(missing_l) > 0) {
      stop("The following labels were not found in stats_df: ",
           paste(missing_l, collapse = ", "),
           "\nAvailable labels: ", paste(available, collapse = ", "))
    }
    stats_df <- stats_df[stats_df$label %in% labels, ]
  }

  all_labels    <- sort(unique(stats_df$label))
  n_labels      <- length(all_labels)
  all_segments  <- sort(unique(stats_df$segment_id))
  n_segments    <- length(all_segments)

  # Summary statistics
  cat("\n")
  cat(rep("=", 60), sep = "")
  cat(sprintf("\n  Feature: %s\n", feature_name))
  cat(sprintf("  Labels  : %s\n", paste(all_labels, collapse = ", ")))
  cat(sprintf("  Segments: %d | Renditions: %d\n",
              n_segments, nrow(stats_df)))
  cat(rep("=", 60), sep = "")
  cat("\n\n")

  # Overall per-label summary (pooling all segments)
  cat("## Per-label summary (pooled across segments)\n\n")

  label_summary <- stats_df |>
    dplyr::group_by(label) |>
    dplyr::summarise(
      n_renditions = dplyr::n(),
      n_segments   = dplyr::n_distinct(segment_id),
      grand_mean   = mean(mean,       na.rm = TRUE),
      grand_sd     = sd(mean,         na.rm = TRUE),
      grand_sem    = grand_sd / sqrt(dplyr::n()),
      grand_median = median(mean,     na.rm = TRUE),
      grand_min    = min(min_val,     na.rm = TRUE),
      grand_max    = max(max_val,     na.rm = TRUE),
      .groups      = "drop"
    ) |>
    dplyr::arrange(label)

  print(label_summary, n = Inf)
  cat("\n")

  # Per-segment x label breakdown
  if (n_segments > 1) {
    cat("## Per-segment breakdown\n\n")
    seg_summary <- stats_df |>
      dplyr::group_by(segment_id, label) |>
      dplyr::summarise(
        n      = dplyr::n(),
        mean   = mean(mean, na.rm = TRUE),
        sd     = sd(mean,   na.rm = TRUE),
        median = median(mean, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(segment_id, label)
    print(seg_summary, n = Inf)
    cat("\n")
  }

  # One-way ANOVA
  anova_result <- NULL

  if (run_anova && n_labels > 1) {
    cat(rep("-", 60), sep = "")
    cat("\n## One-way ANOVA (per segment)\n\n")
    anova_result <- tryCatch(
      anova_analysis(stats_df, plot = FALSE),
      error = function(e) {
        warning("anova_analysis() failed: ", conditionMessage(e))
        NULL
      }
    )
    cat("\n")
  } else if (n_labels == 1) {
    cat("## ANOVA skipped: only one label present ('",
        all_labels, "').\n\n", sep = "")
  } else if (!run_anova) {
    cat("## ANOVA skipped (run_anova = FALSE).\n\n")
  }

  # Plot
  p <- NULL

  if (plot) {
    # Build a tidy data frame: one row per rendition, using 'mean' as the
    # per-rendition summary value.
    plot_df <- stats_df |>
      dplyr::mutate(
        label      = factor(label, levels = all_labels),
        segment_id = factor(segment_id)
      )

    # Colour palette
    n_colors  <- n_labels
    pal_colors <- if (requireNamespace("RColorBrewer", quietly = TRUE) &&
                       n_colors <= RColorBrewer::brewer.pal.info[palette, "maxcolors"]) {
      RColorBrewer::brewer.pal(max(3L, n_colors), palette)[seq_len(n_colors)]
    } else {
      scales::hue_pal()(n_colors)
    }
    names(pal_colors) <- all_labels

    # Base ggplot
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = label, y = mean, fill = label, colour = label)
    )

    # Layer(s) based on plot_type
    if (plot_type %in% c("violin", "both")) {
      p <- p +
        ggplot2::geom_violin(
          alpha    = 0.35,
          colour   = NA,
          trim     = TRUE,
          scale    = "width"
        )
    }

    if (plot_type %in% c("boxplot", "both")) {
      box_fill <- if (plot_type == "both") "white" else NULL
      p <- p +
        ggplot2::geom_boxplot(
          width    = if (plot_type == "both") 0.25 else 0.55,
          outlier.shape = if (jitter && plot_type == "boxplot") NA else 19,
          fill     = box_fill,
          colour   = "grey30",
          alpha    = 0.8
        )
    }

    # Jittered individual points
    if (jitter) {
      p <- p +
        ggplot2::geom_jitter(
          width   = 0.12,
          size    = 1.2,
          alpha   = point_alpha,
          colour  = "grey20"
        )
    }

    # Mean +- SEM cross-bar
    p <- p +
      ggplot2::stat_summary(
        fun      = mean,
        fun.min  = function(x) mean(x) - sd(x) / sqrt(length(x)),
        fun.max  = function(x) mean(x) + sd(x) / sqrt(length(x)),
        geom     = "pointrange",
        shape    = 18,
        size     = 0.7,
        colour   = "grey10"
      )

    # Scales and theme
    p <- p +
      ggplot2::scale_fill_manual(values   = pal_colors) +
      ggplot2::scale_colour_manual(values = pal_colors) +
      ggplot2::labs(
        title    = feature_name,
        subtitle = sprintf(
          "%s | n labels = %d | n segments = %d",
          ifelse(n_labels > 1 && run_anova, "One-way ANOVA applied", ""),
          n_labels, n_segments
        ),
        x        = "Label",
        y        = feature_name,
        fill     = "Label",
        colour   = "Label"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey40"),
        axis.title       = ggplot2::element_text(face = "bold"),
        legend.position  = "right",
        strip.text       = ggplot2::element_text(face = "bold"),
        panel.grid.major.x = ggplot2::element_blank()
      )

    # Optional faceting by segment
    if (facet_segments && n_segments > 1) {
      p <- p +
        ggplot2::facet_wrap(
          ~ segment_id,
          ncol   = ncol,
          labeller = ggplot2::label_both
        )
    }

    print(p)
  }

  # Return
  invisible(list(
    summary = label_summary,
    anova   = anova_result,
    plot    = p
  ))
}

#' @rdname feature_stats
#' @export
feature_stats.Sap <- function(x,
                              feature      = c("fund_freq", "wiener_entropy"),
                              segment_type = c("motifs", "syllables", "segments"),
                              ...) {

  feature      <- match.arg(feature)
  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)   # "motifs" -> "motif"

  # Build the slot name: stats_fund_freq  |  stats_wiener_entropy
  slot_name <- paste0("stats_", feature)

  stats_df <- x$features[[feature_type]][[slot_name]]

  if (is.null(stats_df)) {
    stop(sprintf(
      paste0("Slot 'sap$features$%s$%s' is NULL or does not exist.\n",
             "Run refine_FF() or refine_sh() with stats = TRUE first."),
      feature_type, slot_name
    ))
  }

  # Human-readable feature name for plot labels
  feature_name <- switch(
    feature,
    "fund_freq"      = "Fundamental Frequency (Hz)",
    "wiener_entropy" = "Wiener Entropy"
  )

  feature_stats(stats_df, feature_name = feature_name, ...)
}

