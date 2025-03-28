# Plot Traces -------------------------------------------------------------
# Update date: Mar. 27, 2025

#' Plot Traces of Amplitude Envelope or Fundamental Frequency
#'
#' @description
#' Creates line plot visualizations of amplitude envelopes or fundamental frequency
#' traces from audio segments, supporting multiple visualization options.
#'
#' @param x An object to visualize (matrix or SAP object)
#' @param labels Optional vector of labels to include (default: NULL, uses all labels)
#' @param plot_type Type of plot: "individual", "average", or "combined" (default: "combined")
#' @param feature Type of feature to plot: "env" (amplitude envelope) or "pitch" (fundamental frequency)
#' @param alpha Transparency for individual traces (default: 0.2)
#' @param ncol Number of columns for facet wrapping (default: 1)
#' @param palette Color palette for plotting (default: "Set1")
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables', 'segments')
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For matrices:
#' \itemize{
#'   \item Requires matrix with 'time_window' attribute
#'   \item Column names used as labels
#'   \item Supports direct visualization of pre-computed traces
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Supports multiple segment types
#'   \item Can plot amplitude envelope or fundamental frequency
#'   \item Flexible visualization options
#' }
#'
#' Plot Types:
#' \describe{
#'   \item{individual}{Shows each rendition's trace, faceted by label}
#'   \item{average}{Displays mean trace with standard error}
#'   \item{combined}{Shows both individual traces and mean trace}
#' }
#'
#' @return A ggplot object with the specified trace visualization
#'
#' @examples
#' \dontrun{
#' # Plot amplitude envelope traces from a matrix
#' plot_traces(sap_obj$features$motif$amp_env,
#'             feature = "env",
#'             plot_type = "combined")
#'
#' # Plot fundamental frequency traces from a SAP object
#' plot_traces(sap_obj,
#'             segment_type = "motifs",
#'             feature = "pitch",
#'             plot_type = "individual")
#'
#' # Customize visualization
#' plot_traces(sap_obj$features$motif$amp_env,
#'             labels = c("BL", "Post"),
#'             feature = "env",
#'             plot_type = "average",
#'             alpha = 0.1,
#'             ncol = 2)
#' }
#'
#' @importFrom dplyr mutate select filter arrange left_join bind_rows %>%
#'
#' @export
plot_traces <- function(x, ...) {
  UseMethod("plot_traces")
}


#' @rdname plot_traces
#' @export
plot_traces.default  <- function(x,
                              labels = NULL,
                              plot_type = c("combined", "individual", "average"),
                              feature = c("env", "pitch"),
                              alpha = 0.2,
                              ncol = 1,
                              palette = "Set1"
                              ) {

  # Input validation
  plot_type <- match.arg(plot_type)
  feature <- match.arg(feature)

  # Check if input is matrix
  if (!is.matrix(x)) {
    stop("Input must be a matrix")
  }

  # Check for required attributes
  if (is.null(attr(x, "time_window"))) {
    stop("Matrix must have 'time_window' attribute")
  }

  if (is.null(colnames(x))) {
    stop("Matrix must have column names as labels")
  }

  # Create data frame from matrix and make column names unique
  dat <- as.data.frame(x)
  original_names <- colnames(dat)
  colnames(dat) <- make.unique(colnames(dat))

  # Calculate time points
  time_window <- attr(x, "time_window")
  time_step <- time_window / nrow(x)
  time_points <- seq(0, time_window, by = time_step)[1:nrow(x)]

  # Create mapping dataframe with original and unique names
  mapping_df <- data.frame(
    col_name = colnames(dat),
    label = original_names,
    rendition_no = sequence(rle(original_names)$lengths)
  )

  # Create long format data
  res_long <- dat |>
    dplyr::mutate(time = time_points) |>
    tidyr::pivot_longer(cols = -time, names_to = "col_name", values_to = "value") |>
    dplyr::left_join(mapping_df, by = "col_name") |>
    dplyr::select(-col_name) |>
    dplyr::arrange(rendition_no, time)

  # Filter labels if specified
  if (!is.null(labels)) {
    if (!all(labels %in% unique(res_long$label))) {
      stop("Some provided labels not found in matrix")
    }
    res_long <- res_long |>
      dplyr::filter(label %in% labels)
  }

  # Load mean_se function from ggplot2
  mean_se <- function(x) {
    x <- stats::na.omit(x)
    se <- stats::sd(x) / sqrt(length(x))
    mean <- mean(x)
    data.frame(y = mean,
               ymin = mean - se,
               ymax = mean + se)
  }

  # Set y-axis label based on trace type
  y_label <- switch(feature,
                    "env" = "Amplitude Envelope",
                    "pitch" = "Fundamental Frequency (kHz)")

  # Create plot based on type
  if (plot_type == "individual") {
    p <- ggplot2::ggplot(res_long,
                         ggplot2::aes(x = time, y = value,
                                      group = rendition_no, color = label)) +
      ggplot2::geom_line(alpha = alpha) +
      ggplot2::facet_wrap(~label, ncol = ncol) +
      ggplot2::labs(x = "Time (s)",
                    y = y_label,
                    color = "Label") +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_brewer(palette = palette) +
      ggplot2::ggtitle(paste(y_label, "by Label")) +
      ggplot2::theme(legend.position = "right",
                     strip.text = ggplot2::element_text(size = 12, face = "bold"))

  } else if (plot_type == "average") {
    p <- ggplot2::ggplot(res_long,
                         ggplot2::aes(x = time, y = value,
                                      color = label, fill = label)) +
      ggplot2::stat_summary(fun.data = mean_se, geom = "ribbon",
                            alpha = 0.3, color = NA) +
      ggplot2::stat_summary(fun = mean, geom = "line", size = 1) +
      ggplot2::labs(x = "Time (s)", y = y_label) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_brewer(palette = palette) +
      ggplot2::scale_fill_brewer(palette = palette) +
      ggplot2::ggtitle(paste("Mean", y_label, "with Standard Error"))

  } else if (plot_type == "combined") {
    # Create combined data for both plots
    res_long <- res_long |>
      dplyr::mutate(plot_type = "Individual Traces")

    res_long_mean <- res_long |>
      dplyr::mutate(plot_type = "Mean \u00B1 SE")

    res_combined <- dplyr::bind_rows(res_long, res_long_mean)

    p <- ggplot2::ggplot(res_combined,
                         ggplot2::aes(x = time, y = value,
                                      color = label, fill = label)) +
      ggplot2::geom_line(data = . %>% dplyr::filter(plot_type == "Individual Traces"),
                         ggplot2::aes(group = interaction(label, rendition_no)),
                         alpha = alpha) +
      ggplot2::stat_summary(data = . %>% dplyr::filter(plot_type == "Mean \u00B1 SE"),
                            fun.data = mean_se, geom = "ribbon",
                            alpha = 0.3, color = NA) +
      ggplot2::stat_summary(data = . %>% dplyr::filter(plot_type == "Mean \u00B1 SE"),
                            fun = mean, geom = "line", size = 1) +
      ggplot2::facet_wrap(~plot_type, ncol = ncol) +
      ggplot2::labs(x = "Time (s)",
                    y = y_label,
                    color = "Label",
                    fill = "Label") +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_brewer(palette = palette) +
      ggplot2::scale_fill_brewer(palette = palette) +
      ggplot2::ggtitle(y_label) +
      ggplot2::theme(legend.position = "right",
                     strip.text = ggplot2::element_text(size = 12, face = "bold"))
  } else {
    stop("Invalid plot_type. Choose 'individual', 'average', or 'combined'")
  }

  return(p)
}

#' @rdname plot_traces
#' @export
plot_traces.Sap <- function(x,
                            segment_type = c("motifs", "syllables", "segments"),
                            feature = c("env", "pitch"),
                            labels = NULL,
                            plot_type = c("combined", "individual", "average"),
                            alpha = 0.2,
                            ncol = 1,
                            palette = "Set1",
                            ...) {

  # Input validation
  segment_type <- match.arg(segment_type)
  feature <- match.arg(feature)
  plot_type <- match.arg(plot_type)

  # Update feature type by removing 's' from segment_type
  feature_type <- sub("s$", "", segment_type)

  # Determine the matrix based on segment and feature type
  if (feature == "env") {
    matrix_to_plot <- x$features[[feature_type]]$amp_env
  } else if (feature == "pitch") {
    matrix_to_plot <- x$features[[feature_type]]$fund_freq
  } else {
    stop("Invalid feature type")
  }

  # Call the default method with the extracted matrix
  p <- plot_traces(matrix_to_plot,
              labels = labels,
              plot_type = plot_type,
              feature = feature,
              alpha = alpha,
              ncol = ncol,
              palette = palette,
              ...)

  print(p)

  invisible(x)
}
