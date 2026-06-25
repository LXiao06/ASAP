# Trajectory Analysis Utilities
# Update date : Jun. 23, 2026

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
#' Internal helper used by trajectory plotting and similarity functions to
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
