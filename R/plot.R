# Plot Traces -------------------------------------------------------------
# Update date: Mar. 27, 2025

#' Plot Traces of Song Features
#'
#' @description
#' Creates line plot visualizations of various acoustic features from audio segments,
#' including amplitude envelope, fundamental frequency, pitch goodness, and
#' Wiener entropy. The function supports multiple visualization options.
#'
#' @param x An object to visualize (matrix or SAP object)
#' @param labels Optional vector of labels to include (default: NULL, uses all labels)
#' @param plot_type Type of plot: "individual", "average", or "combined" (default: "combined")
#' @param feature Type of feature to plot: "env" (amplitude envelope), "pitch" (fundamental frequency),
#'        "goodness" (pitch goodness), or "entropy" (Wiener entropy)
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
#'   \item Supports multiple segment types(motifs, syllables, segments)
#'   \item Can plot multiple acoustic features (amplitude envelope, fundamental frequency,
#'         pitch goodness, and Wiener entropy)
#'   \item Flexible visualization options
#' }
#'
#' Features:
#' \describe{
#'   \item{env}{Amplitude envelope - represents the overall loudness/intensity of sound over time}
#'   \item{pitch}{Fundamental frequency (in kHz) - represents the pitch contour}
#'   \item{goodness}{Pitch goodness - represents the periodicity/quality of pitch estimation}
#'   \item{entropy}{Wiener entropy - measures the width and uniformity of the power spectrum}
#' }
#'
#' Plot Types:
#' \describe{
#'   \item{individual}{Shows each rendition's trace, faceted by label}
#'   \item{average}{Displays mean trace with standard error}
#'   \item{combined}{Shows both individual traces and mean trace}
#' }
#'
#'
#' @return A ggplot object with the specified trace visualization. For SAP objects,
#'         the function returns the input object invisibly after printing the plot.
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
#' # Plot Wiener entropy traces
#' plot_traces(sap_obj,
#'             segment_type = "syllables",
#'             feature = "entropy",
#'             plot_type = "average")
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
#' @export
plot_traces <- function(x, ...) {
  UseMethod("plot_traces")
}


#' @rdname plot_traces
#' @export
plot_traces.default  <- function(x,
                              labels = NULL,
                              plot_type = c("combined", "individual", "average"),
                              feature = c("env", "pitch", "goodness", "entropy"),
                              alpha = 0.2,
                              ncol = 1,
                              palette = "Set1",
                              ...
                              ) {
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' needed for this function to work. Please install it.",
         call. = FALSE)
  }

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
    dplyr::select(-.data$col_name) |>
    dplyr::arrange(.data$rendition_no, .data$time)

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
                    "pitch" = "Fundamental Frequency (kHz)",
                    "goodness" = "Pitch Goodness",
                    "entropy" = "Wiener Entropy")

  # Create plot based on type
  if (plot_type == "individual") {
    p <- ggplot2::ggplot(res_long,
                         ggplot2::aes(x = .data$time, y = .data$value,
                                      group = .data$rendition_no, color = label)) +
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
                         ggplot2::aes(x = .data$time, y = .data$value,
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
                         ggplot2::aes(x = .data$time, y = .data$value,
                                      color = label, fill = label)) +
      ggplot2::geom_line(data = function(x) dplyr::filter(plot_type == "Individual Traces"),
                         ggplot2::aes(group = interaction(label, .data$rendition_no)),
                         alpha = alpha) +
      ggplot2::stat_summary(data = function(x) dplyr::filter(plot_type == "Mean \u00B1 SE"),
                            fun.data = mean_se, geom = "ribbon",
                            alpha = 0.3, color = NA) +
      ggplot2::stat_summary(data = function(x) dplyr::filter(plot_type == "Mean \u00B1 SE"),
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
                            feature = c("env", "pitch", "goodness","entropy"),
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

  # Check if the feature type exists in x$features
  if (is.null(x$features) || is.null(x$features[[feature_type]])) {
    stop(sprintf("Feature type '%s' not found in the Sap object", feature_type))
  }

  # Determine the matrix based on segment and feature type
  if (feature == "env") {
    matrix_name <- "amp_env"
  } else if (feature == "pitch") {
    matrix_name <- "fund_freq"
  } else if (feature == "goodness") {
    matrix_name <- "pitch_goodness"
  } else if (feature == "entropy") {
    matrix_name <- "weiner_entropy"
  } else {
    stop("Invalid feature type")
  }

  # Check if the specific feature matrix exists
  if (is.null(x$features[[feature_type]][[matrix_name]])) {
    stop(sprintf("Feature '%s' not available for '%s'. Check if this feature has been calculated.",
                 feature, segment_type))
  }

  # Get the matrix to plot
  matrix_to_plot <- x$features[[feature_type]][[matrix_name]]

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

# Plot clusters -------------------------------------------------------------
# Update date: Apr. 15, 2025

#' Plot Clusters in Time Series Data
#'
#' @description
#' Generic function for plotting clusters in time series data, primarily designed for
#' visualizing bird song syllable clusters and their manual labels. Supports both
#' numeric cluster IDs (pre manual labeling) and alphabetic syllable labels (post manual labeling).
#'
#' @param x Object to plot clusters from. Can be either:
#'   * A matrix with cluster/syllable data
#'   * A Sap object containing bird song analysis data at motif and syllable levels
#' @param data_type Type of data to analyze (for Sap objects):
#'   * "segment": Uses segment-level features and clusters
#'   * "syllable": Uses syllable-level features and clusters
#' @param label_type Type of labels to visualize (for Sap objects):
#'   * "pre": Shows numerical cluster IDs from automatic clustering
#'   * "post": Shows alphabetic syllable labels assigned via manual_label()
#' @param time_resolution Number of time points for plotting (default: 1000)
#' @param cluster_colors Optional named vector of colors for clusters
#' @param sample_percent Optional percentage of data to sample
#' @param balanced Whether to balance samples across labels
#' @param labels Optional vector of experimental condition labels to subset (e.g., "BL", "Post", "Rec")
#' @param motif_clusters Optional vector of specific clusters to include
#' @param ordered Whether to order motifs by UMAP coordinates
#' @param descending Direction of UMAP-based ordering
#' @param seed Random seed for reproducibility
#' @param cores Number of cores for parallel processing
#' @param verbose Whether to print progress messages
#' @param main Plot title (for matrix method)
#' @param ... Additional arguments passed between methods
#'
#' @details
#' The function provides two main visualization approaches:
#'
#' 1. Pre-labeling visualization (label_type = "pre"):
#'    * Shows numerical cluster IDs from automatic clustering
#'    * Uses data from x$features$segment/syllable$feat.embeds
#'    * Clusters are represented by numbers (1, 2, 3, etc.)
#'    * Useful for evaluating automatic clustering results
#'
#' 2. Post-labeling visualization (label_type = "post"):
#'    * Shows alphabetic syllable labels from manual annotation
#'    * Uses data from x$syllables (created by manual_label())
#'    * Syllables are represented by letters (a, b, c, etc.)
#'    * Useful for viewing manual syllable classifications
#'
#' For matrix input:
#' * Must have column names as labels
#' * Requires a 'time_window' attribute
#' * Can contain either numeric clusters or character syllable labels
#'
#' For Sap objects:
#' * data_type determines the feature level for analysis:
#'   - "segment": Fine-grained analysis of song segments
#'   - "syllable": Analysis at the syllable level
#' * Supports experimental condition labels (BL, Post, Rec)
#' * Can order motifs by UMAP coordinates for pattern visualization
#' * Allows balanced sampling across conditions
#'
#' @return Invisibly returns x
#'
#' @examples
#' \dontrun{
#' # Matrix method example
#' mat <- matrix(sample(1:5, 1000, replace = TRUE), ncol = 10)
#' attr(mat, "time_window") <- 1.2
#' colnames(mat) <- rep(c("BL", "Post"), each = 5)
#' plot_clusters(mat)
#'
#' # Sap object - Pre-labeling (numerical clusters)
#' # Basic cluster visualization
#' plot_clusters(sap,
#'             data_type = "syllable",
#'             label_type = "pre")
#'
#' # Ordered by UMAP with specific conditions
#' plot_clusters(sap,
#'             data_type = "syllable",
#'             label_type = "pre",
#'             ordered = TRUE,
#'             labels = c("BL", "Post"))
#'
#' # Sap object - Post-labeling (syllable letters)
#' # After running manual_label()
#' plot_clusters(sap,
#'             data_type = "syllable",
#'             label_type = "post")
#'
#' # Balanced sampling across conditions
#' plot_clusters(sap,
#'             data_type = "syllable",
#'             label_type = "post",
#'             balanced = TRUE,
#'             labels = c("BL", "Post", "Rec"))
#' }
#'
#' @seealso
#' * manual_label() for creating syllable labels
#' * auto_label() for automatic cluster generation
#'
#' @export
plot_clusters <- function(x, ...) {
  UseMethod("plot_clusters")
}

#' @rdname plot_clusters
#' @export
plot_clusters.matrix <- function(x,
                                labels = NULL,
                                cluster_colors = NULL,
                                main = NULL,
                                ...) {
  # Validations
  if (!is.matrix(x)) stop("Input must be a matrix")
  if (is.null(attr(x, "time_window"))) stop("Matrix must have 'time_window' attribute")
  if (is.null(colnames(x))) stop("Matrix must have column names as labels")

  # Store time_window before subsetting
  time_window <- attr(x, "time_window")
  if (!is.numeric(time_window) || length(time_window) != 1) {
    stop("time_window attribute must be a single numeric value")
  }

  # Handle label subsetting if provided
  if (!is.null(labels)) {
    if (!all(labels %in% colnames(x))) {
      stop("Some provided labels not found in matrix")
    }
    x <- x[, colnames(x) %in% labels]
    attr(x, "time_window") <- time_window  # Preserve time_window attribute
  }

  # Create reversed matrix for plotting
  reversed_cluster_matrix <- x[, ncol(x):1]

  # Calculate positions using reversed order
  current_labels <- colnames(x)
  reversed_labels <- rev(unique(current_labels))

  ordered_labels <- factor(current_labels, levels = unique(current_labels))
  samples_per_label <- rev(table(ordered_labels))  # rev() for bottom-to-top plotting

  cumulative_positions <- cumsum(c(0, head(samples_per_label, -1)))
  label_positions <- cumulative_positions + samples_per_label/2
  hline_positions <- cumsum(samples_per_label)[-length(samples_per_label)]

  # Set up colors and handle unique values
  unique_values <- sort(unique(na.omit(as.vector(x))))

  if (is.null(cluster_colors)) {
    n_colors <- length(unique_values)
    palette_colors <- Polychrome::createPalette(
      n_colors + 2,
      seedcolors = c("#ffffff", "#000000"),
      range = c(10, 90)
    )[-(1:2)]
    cluster_colors <- setNames(palette_colors, as.character(unique_values))
  }

  # Convert character values to numeric for plotting
  if (!is.numeric(x)) {
    value_mapping <- setNames(seq_along(unique_values), unique_values)
    numeric_matrix <- matrix(value_mapping[reversed_cluster_matrix],
                             nrow = nrow(reversed_cluster_matrix))
    at_values <- seq(0.5, length(unique_values) + 0.5, by = 1)
    label_values <- unique_values
  } else {
    numeric_matrix <- reversed_cluster_matrix
    at_values <- seq(min(unique_values) - 0.5,
                     max(unique_values) + 0.5,
                     by = 1)
    label_values <- paste("Cluster", unique_values)
  }

  # Create heatmap
  heatmap <- lattice::levelplot(numeric_matrix,
                                col.regions = cluster_colors,
                                at = at_values,
                                border = "transparent",
                                ylab = "Labels",
                                xlab = "Time (s)",
                                main = main,
                                scales = list(
                                  x = list(
                                    at = seq(0, nrow(numeric_matrix), length.out = 5),
                                    labels = sprintf("%.1f", seq(0, time_window, length.out = 5))
                                  ),
                                  y = list(
                                    at = label_positions,
                                    labels = reversed_labels
                                  )
                                ),
                                panel = function(x, y, z, ...) {
                                  lattice::panel.levelplot(x, y, z, ...)
                                  lattice::panel.abline(h = hline_positions,
                                                        col = "white",
                                                        lwd = 2)
                                },
                                aspect = "fill",
                                colorkey = list(
                                  space = "right",
                                  width = 1,
                                  height = 0.75,
                                  labels = list(
                                    at = seq_along(unique_values),
                                    labels = label_values
                                  )
                                ))

  print(heatmap)

  invisible(NULL)
}


#' @rdname plot_clusters
#' @export
plot_clusters.Sap <- function(x,
                             data_type = c("segment", "syllable"),
                             label_type = c("pre", "post"),
                             time_resolution = 1000,  # Number of time points
                             cluster_colors = NULL, # Optional custom color palette
                             sample_percent = NULL,
                             balanced = FALSE,
                             labels = NULL,
                             motif_clusters = NULL,
                             ordered = FALSE,
                             descending = TRUE,
                             seed = 222,
                             cores = NULL,
                             verbose = TRUE,
                             ...) {
  # Match and validate arguments
  data_type <- match.arg(data_type)
  label_type <- match.arg(label_type)

  if(verbose) {
    message(sprintf("\n=== Starting Plotting %s in Motif ===\n",
                    tools::toTitleCase(data_type)))
  }

  # Select appropriate data and plot column based on label_type
  if(label_type == "pre") {
    feat_embeds <- switch(data_type,
                          "segment" = x$features$segment$feat.embeds,
                          "syllable" = x$features$syllable$feat.embeds)
    plot_col <- "cluster"

  } else {  # post-labeling
    if(is.null(x$syllables)) {
      stop("No syllable labels found. Run manual_label() first.")
    }
    feat_embeds <- x$syllables
    plot_col <- "syllable"
  }

  # Validation for feature embeddings
  if(is.null(feat_embeds)) {
    stop(sprintf("Feature embeddings for %s required", data_type))
  }

  if(!plot_col %in% names(feat_embeds)) {
    stop(sprintf("No %s column found in the data", plot_col))
  }

  # Use original motifs if no ordering/cluster filtering needed
  if (!ordered && is.null(motif_clusters)) {
    segments_df <- x[["motifs"]]
  } else {
    if (is.null(x$features$motif$feat.embeds)) {
      stop("Feature embeddings required for ordered/clustered motif plots")
    }

    segments_df <- x$features$motif$feat.embeds |>
      as_segment()

    # Apply UMAP-based ordering if requested
    if (ordered) {
      segments_df <- segments_df |>
        dplyr::arrange(
          day_post_hatch,
          if (descending) dplyr::desc(UMAP2) else UMAP2,
          if (descending) dplyr::desc(UMAP1) else UMAP1
        )
    }
  }

  # Validation
  if (!inherits(segments_df, "segment") || nrow(segments_df) == 0) {
    stop("No segments found in the specified segment type")
  }

  # Select and balance segments
  segments_df <- select_segments(segments_df,
                                 labels = labels,
                                 clusters =  motif_clusters,
                                 balanced = balanced,
                                 sample_percent = sample_percent,
                                 seed = seed)

  # Check if segments_df is empty after subsetting
  if (nrow(segments_df) == 0) {
    stop("No segments remaining after subsetting. Check labels/clusters.")
  }

  # Remove UMAP and cluster columns after ordering if they exist
  cols_to_remove <- intersect(c("UMAP1", "UMAP2", "cluster"), names(segments_df))
  if (length(cols_to_remove) > 0) {
    segments_df <- segments_df |>
      dplyr::select(-dplyr::all_of(cols_to_remove))
  }

  # Prepare df for cluster plotting
  motifs <- segments_df |>
    dplyr::mutate(original_order = dplyr::row_number()) |>
    dplyr::rename(motif_start = start_time, motif_end = end_time)|>
    dplyr::group_by(filename) |>
    dplyr::mutate(motif_index = dplyr::row_number()) |>
    dplyr::ungroup()

  # Join with appropriate feature embeddings
  df <- feat_embeds |>
    dplyr::inner_join(motifs, by = c("filename", "day_post_hatch", "label"),
                      relationship = "many-to-many") |>
    dplyr::filter(start_time >= .data$motif_start & start_time <= .data$motif_end) |>
    dplyr::mutate(
      #segment_duration = end_time - start_time,
      .start = start_time - .data$motif_start,
      .end = end_time - .data$motif_start,
      motif_duration = .data$motif_end - .data$motif_start,
      motif_label = paste(filename, .data$motif_index))  |>
    dplyr::arrange(.data$original_order)   # Sort by original order
  # |>dplyr::select(filename, day_post_hatch, label, motif_index, .start, .end,
  #                 cluster, motif_duration, original_order, motif_label)

  # Get ordered unique motifs based on original UMAP order
  motif_order <- df |>
    dplyr::distinct(.data$motif_label, .data$original_order) |>
    dplyr::arrange(.data$original_order) |>
    dplyr::pull(.data$motif_label)

  # Check if motif_duration is unique
  unique_durations <- unique(df$motif_duration)
  time_window <- if (length(unique_durations) == 1) unique_durations else max(unique_durations)

  # Function to process each motif
  process_motif <- function(motif_id) {
    motif_rows <- dplyr::filter(df, .data$motif_label == motif_id)
    motif_vector <- rep(NA, time_resolution)

    for (j in seq_len(nrow(motif_rows))) {
      start_idx <- max(1, round(motif_rows$.start[j] / motif_rows$motif_duration[j] * time_resolution))
      end_idx <- min(time_resolution, round(motif_rows$.end[j] / motif_rows$motif_duration[j] * time_resolution))
      motif_vector[start_idx:end_idx] <- motif_rows[[plot_col]][j]
    }
    return(motif_vector)
  }

  # Use parallel processing
  cluster_list <- parallel_apply(motif_order, process_motif, cores = cores)

  # Convert list to matrix
  cluster_matrix <- do.call(cbind, cluster_list)

  # Store attributes in amp_matrix
  attr(cluster_matrix, "time_window") <- time_window

  # Set column names as labels, preserving order
  colnames(cluster_matrix) <- sapply(motif_order,
                                     function(m) df$label[df$motif_label == m][1])

  # Plot using matrix method with appropriate title
  title <- ifelse(label_type == "pre",
                  sprintf("%s cluster heatmap", data_type),
                  "syllable heatmap")

  result <- plot_clusters.matrix(cluster_matrix,
                                cluster_colors = cluster_colors,
                                main =  title)
  invisible(x)
}
