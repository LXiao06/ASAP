
# Run UMAP ----------------------------------------------------------------
# Update date : Feb. 7, 2025

#' Run UMAP Dimensionality Reduction
#'
#' @description
#' A generic function to perform UMAP dimensionality reduction on feature data.
#'
#' @param x An object to analyze, either a data frame or SAP object
#' @param metadata_cols Column indices for metadata (for default method)
#' @param scale Whether to scale features before UMAP
#' @param n_neighbors Number of neighbors (default: 15)
#' @param n_components Number of output dimensions (default: 2)
#' @param min_dist Minimum distance parameter (default: 0.1)
#' @param seed Random seed for reproducibility
#' @param n_threads Number of computation threads
#' @param verbose Whether to print progress messages
#' @param segment_type For SAP objects: Type of segments to analyze ('motifs', 'syllables', 'bouts', 'segments')
#' @param data_type For SAP objects: Type of feature data ('spectral_feature','spectrogram', 'traj_mat')
#' @param label For SAP objects: Specific label to filter data
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' This generic function supports UMAP analysis through two methods:
#' \itemize{
#'   \item Default method for feature data frames
#'   \item SAP object method for organized song features
#' }
#'
#' @return
#' For default method: Matrix of UMAP coordinates
#' For SAP objects: Updated SAP object with UMAP coordinates stored in features slot
#'
#' @examples
#' \dontrun{
#' # Run UMAP on feature data frame
#' coords <- run_umap(features,
#'                    metadata_cols = c(1:5),
#'                    n_neighbors = 15)
#'
#' # Run UMAP on SAP object
#' sap_obj <- run_umap(sap_object,
#'                     segment_type = "motifs",
#'                     data_type = "spectral_feature")
#'
#' # UMAP with specific parameters
#' coords <- run_umap(features,
#'                    metadata_cols = 1:3,
#'                    scale = TRUE,
#'                    n_neighbors = 20,
#'                    seed = 123)
#'
#' # UMAP with label filtering
#' sap_obj <- run_umap(sap_obj,
#'                     segment_type = "syllables",
#'                     data_type = "spectral_feature",
#'                     label = "a")
#' }
#'
#' @rdname run_umap
#' @export
run_umap <- function(x, ...) {
  UseMethod("run_umap")
}

#' @rdname run_umap
#' @export
run_umap.default <- function(x,
                             metadata_cols = NULL,  # Make it optional
                             scale = TRUE,
                             n_neighbors = 15,
                             n_components = 2,
                             min_dist = 0.1,
                             seed = 222,
                             n_threads = NULL,
                             verbose = TRUE,
                             ...) {

  # Validate input
  if(!is.data.frame(x)) stop("Input data must be a data frame")

  # Handle case when metadata_cols is provided
  if(!is.null(metadata_cols)) {
    if(any(metadata_cols < 1) || any(metadata_cols > ncol(x))) {
      stop("metadata_cols must be between 1 and the number of columns in the data frame")
    }

    # Print column information if metadata exists
    if(verbose) {
      message("Metadata columns:")
      print(colnames(x)[metadata_cols])

      feature_cols <- colnames(x)[-metadata_cols]
      total_feature_cols <- length(feature_cols)

      message("\nFeature columns (showing first 30 of ", total_feature_cols, "):")

      if(total_feature_cols > 30) {
        print(head(feature_cols, 30))
        message("... and ", total_feature_cols - 30, " more columns not shown")
      } else {
        print(feature_cols)
      }
    }

    # Extract features excluding metadata columns
    features <- x[,-metadata_cols]
  } else {
    # If no metadata columns, use all columns as features
    if(verbose) {
      message("No metadata columns specified. Using all columns as features:")
      print(colnames(x))
    }
    features <- x
  }

  if(is.null(n_threads)){
    n_threads <- parallel::detectCores() - 1
  }

  # Convert to matrix if needed
  if (is.data.frame(features)) {
    feature_matrix <- as.matrix(features)
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Run UMAP
  result <- uwot::umap2(
    X = feature_matrix,
    scale = scale,
    n_neighbors = n_neighbors,
    n_components = n_components,
    min_dist = min_dist,
    n_threads = n_threads,
    verbose = verbose
  )

  return(result)
}

#' @rdname run_umap
#' @export
run_umap.Sap <- function(x,
                         segment_type = c("motifs", "syllables", "bouts", "segments"),
                         data_type = c("spectral_feature", "spectrogram","traj_mat"),
                         label = NULL,
                         scale = TRUE,
                         n_neighbors = 20,
                         n_components = 2,
                         min_dist = 0.1,
                         seed = 222,
                         n_threads = NULL,
                         verbose = TRUE,
                         ...) {
  if(verbose) message(sprintf("\n=== Starting Run UMAP for %s using %s ===\n",
                              segment_type[1],
                              data_type[1]))

  # Validate input
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  # Match arguments
  segment_type <- match.arg(segment_type)
  data_type <- match.arg(data_type)

  # Get feature type
  feature_type <- sub("s$", "", segment_type)

  # Get appropriate data based on data_type
  if(data_type %in% c("spectral_feature", "spectrogram")) {
    feature_data <- x$features[[feature_type]][[data_type]]
    if (is.null(feature_data)) {
      stop(sprintf("%s features not found in SAP object", segment_type))
    }

    # Handle label filtering and metadata columns as before
    if (!is.null(label)) {
      if (!"label" %in% names(feature_data)) {
        stop("Label column not found in feature data")
      }
      available_labels <- unique(feature_data$label)
      if (!label %in% available_labels) {
        stop(sprintf("Label '%s' not found. Available labels: %s",
                     label, paste(available_labels, collapse = ", ")))
      }
      feature_data <- feature_data[feature_data$label == label, ]
    }

    # Identify metadata columns
    standard_meta <- c("filename", "day_post_hatch", "label", "start_time", "end_time")
    metadata_cols <- which(names(feature_data) %in% standard_meta)
    if ("duration" %in% names(feature_data) && length(unique(feature_data$duration)) == 1) {
      metadata_cols <- c(metadata_cols, which(names(feature_data) == "duration"))
    }

  } else if(data_type == "traj_mat") {
    feature_data <- x$features[[feature_type]][[data_type]]
    if (is.null(feature_data)) {
      stop(sprintf("%s trajectory matrix not found in SAP object", segment_type))
    }
    metadata_cols <- NULL
  }

  # Get UMAP coordinates
  umap_coords <- run_umap.default(
    x = feature_data,
    metadata_cols = metadata_cols,
    scale = scale,
    n_neighbors = n_neighbors,
    n_components = n_components,
    min_dist = min_dist,
    seed = seed,
    n_threads = n_threads,
    verbose = verbose,
    ...
  )

  # Store UMAP parameters
  params_list <- list(
    data_type = data_type,
    scale = scale,
    n_neighbors = n_neighbors,
    n_components = n_components,
    min_dist = min_dist,
    seed = seed,
    label = label
  )

  # Handle results based on data_type
  if(data_type %in% c("spectral_feature", "spectrogram")) {
    # Create UMAP result with metadata
    umap_result <- data.frame(
      feature_data[, metadata_cols, drop = FALSE],
      UMAP1 = umap_coords[,1],
      UMAP2 = umap_coords[,2]
    )

    # Merge with clusters if they exist
    if (!is.null(x$features[[feature_type]][["clusters"]])) {
      feat.embeds <- merge(
        x$features[[feature_type]][["clusters"]],
        umap_result,
        by = names(feature_data)[metadata_cols],
        all = TRUE
      )
      x$features[[feature_type]][["feat.embeds"]] <- feat.embeds
      x$features[[feature_type]][["clusters"]] <- NULL
    } else {
      x$features[[feature_type]][["feat.embeds"]] <- umap_result
    }

    # Add parameters as attributes to feat.embeds
    attr(x$features[[feature_type]][["feat.embeds"]], "umap_params") <- params_list

  } else {
    # For trajectory matrix, add UMAP coordinates to existing traj.embeds
    traj_embeds <- x$features[[feature_type]][["traj.embeds"]]
    traj_embeds$PC1 <-feature_data[,1]
    traj_embeds$PC2 <- feature_data[,2]
    traj_embeds$UMAP1 <- umap_coords[,1]
    traj_embeds$UMAP2 <- umap_coords[,2]

    # Add parameters as attributes to traj.embeds
    attrs <- attributes(traj_embeds)  # Preserve existing attributes
    attrs$umap_params <- params_list
    x$features[[feature_type]][["traj.embeds"]] <- traj_embeds
    attributes(x$features[[feature_type]][["traj.embeds"]]) <- attrs
  }

  # Print access information
  if(verbose) {
    cat("\nUMAP completed successfully!")
    if(data_type == "spectral_feature") {
      cat(sprintf("\nAccess results via: sap$features$%s$feat.embeds", feature_type))
      cat("\nAccess parameters via: attributes(sap$features$", feature_type, "$feat.embeds)$umap_params", sep="")
    } else {
      cat(sprintf("\nAccess results via: sap$features$%s$traj.embeds", feature_type))
      cat("\nAccess parameters via: attributes(sap$features$", feature_type, "$traj.embeds)$umap_params", sep="")
    }
  }

  invisible(x)
}


# Plot UMAP ---------------------------------------------------------------
# Update date : Feb. 7, 2025

#' Plot UMAP Visualization
#'
#' @description
#' Creates customizable UMAP visualizations with options for grouping, highlighting, and customization.
#'
#' @param x An object to visualize, either a data frame with UMAP coordinates or a SAP object
#' @param dims UMAP dimensions to plot (default: c("UMAP1", "UMAP2"))
#' @param group.by Column name for grouping points
#' @param split.by Column name for faceting plots
#' @param subset.by Column name for subsetting data
#' @param subset.value Values to subset by
#' @param cols Custom colors for groups
#' @param pt.size Point size (default: 0.5)
#' @param stroke Point stroke width (default: 0.5)
#' @param alpha Point transparency (default: 0.3)
#' @param highlight.alpha Transparency for highlighted points
#' @param label Whether to add labels (default: FALSE)
#' @param label.size Size of labels (default: 4)
#' @param repel Whether to use repelling labels (default: FALSE)
#' @param highlight.by Column name for highlighting
#' @param highlight.value Values to highlight
#' @param cols.highlight Colors for highlighted points (default: '#DE2D26')
#' @param sizes.highlight Size for highlighted points (default: 1)
#' @param background.value Background group value
#' @param na.value Color for NA values (default: 'grey80')
#' @param ncol Number of columns in multi-plot layout
#' @param combine Whether to combine multiple plots (default: TRUE)
#' @param segment_type For SAP objects: Type of segments to visualize ('motifs', 'syllables', 'bouts', 'segments')
#' @param verbose For SAP objects: Whether to print progress messages
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' This function creates UMAP visualizations with the following features:
#' \itemize{
#'   \item Flexible grouping and highlighting
#'   \item Customizable point appearance
#'   \item Optional labels and faceting
#'   \item Multiple plot combinations
#' }
#'
#' @return
#' For default method: A ggplot object or list of plots
#' For SAP objects: Updated SAP object with plot as side effect
#'
#' @examples
#' \dontrun{
#' # Basic UMAP plot from data frame
#' plot_umap(umap_df, group.by = "cluster")
#'
#' # Plot with highlighting
#' plot_umap(umap_df,
#'           group.by = "cluster",
#'           highlight.by = "label",
#'           highlight.value = "a")
#'
#' # Plot with faceting
#' plot_umap(umap_df,
#'           group.by = "cluster",
#'           split.by = "day_post_hatch")
#'
#' # Plot from SAP object
#' plot_umap(sap_obj,
#'           segment_type = "motifs",
#'           group.by = "label")
#'
#' # SAP object plot with custom grouping and highlighting
#' plot_umap(sap_obj,
#'           segment_type = "syllables",
#'           group.by = "label",
#'           highlight.by = "cluster",
#'           highlight.value = c(1, 2))
#' }
#'
#' @rdname plot_umap
#' @export
plot_umap <- function(x, ...) {
  UseMethod("plot_umap")
}

#' @rdname plot_umap
#' @export
plot_umap.default <- function(x,
                              dims = c("UMAP1", "UMAP2"),
                              subset.by = NULL,
                              subset.value = NULL,
                              group.by = NULL,
                              split.by = NULL,
                              cols = NULL,
                              pt.size = 0.5,
                              stroke = 0.5,
                              alpha = 0.3,
                              highlight.alpha = NULL,
                              label = FALSE,
                              label.size = 4,
                              repel = FALSE,
                              highlight.by = NULL,
                              highlight.value = NULL,
                              cols.highlight = '#DE2D26',
                              sizes.highlight = 1,
                              background.value = NULL,
                              na.value = 'grey80',
                              ncol = NULL,
                              combine = TRUE,
                              ...) {
  # Validate input
  if(!is.data.frame(x)) stop("Input must be a data frame")

  # Create a copy of the input data frame
  data <- x

  # Validate and handle subsetting
  if (!is.null(subset.by) && !is.null(subset.value)) {
    if (!subset.by %in% colnames(data)) {
      stop(paste("Cannot find", subset.by, "in data"))
    }

    # Check if all subset values exist in the data
    existing_values <- unique(data[[subset.by]])
    missing_values <- setdiff(subset.value, existing_values)
    if (length(missing_values) > 0) {
      stop(sprintf("The following values were not found in column '%s': %s",
                   subset.by,
                   paste(missing_values, collapse = ", ")))
    }

    # Apply subsetting
    data <- data[data[[subset.by]] %in% subset.value, ]
    if (nrow(data) == 0) {
      stop("No data remaining after subsetting")
    }
  }

  # Validate highlight values before proceeding
  if (!is.null(highlight.by) && !is.null(highlight.value)) {
    if (!highlight.by %in% colnames(data)) {
      stop(paste("Cannot find", highlight.by, "in data"))
    }
    # Check if all highlight values exist in the data
    existing_values <- unique(data[[highlight.by]])
    missing_values <- setdiff(highlight.value, existing_values)
    if (length(missing_values) > 0) {
      stop(sprintf("The following highlight values were not found in column '%s': %s",
                   highlight.by,
                   paste(missing_values, collapse = ", ")))
    }
  }

  # Set default grouping if not specified
  group.by <- group.by %||% "cluster"
  if (!is.null(group.by) && !all(group.by %in% colnames(data))) {
    stop("Cannot find grouping variable(s) in data")
  }

  # Generate colors if not provided
  if (is.null(cols) && !is.null(group.by)) {
    unique_levels <- unique(as.character(data[[group.by[1]]]))
    unique_levels <- sort(as.numeric(unique_levels))
    n_colors <- length(unique_levels)
    palette_colors <- Polychrome::createPalette(
      n_colors + 2,
      seedcolors = c("#ffffff", "#000000"),
      range = c(10, 90)
    )[-(1:2)]
    cols <- setNames(palette_colors, as.character(unique_levels))
  }

  # Ensure cluster levels are sorted numerically
  data$cluster <- factor(as.character(data$cluster), levels = as.character(sort(as.numeric(unique(as.character(data$cluster))))))

  # Create plots
  plots <- lapply(group.by, function(x) {
    if (!is.factor(data[[x]])) {
      data[[x]] <- factor(data[[x]], levels = sort(as.numeric(unique(as.character(data[[x]])))))
    }

    p <- plot_single_umap(
      data = data,
      dims = dims,
      cols = cols,
      pt.size = pt.size,
      stroke = stroke,
      alpha = alpha,
      highlight.alpha = highlight.alpha,
      label = label,
      repel = repel,
      label.size = label.size,
      highlight.by = highlight.by,
      highlight.value = highlight.value,
      cols.highlight = cols.highlight,
      sizes.highlight = sizes.highlight,
      background.value = background.value,
      na.value = na.value
    )

    if (!is.null(split.by)) {
      if (!split.by %in% colnames(data)) {
        stop(paste("Cannot find", split.by, "in data"))
      }
      if (!is.factor(data[[split.by]])) {
        data[[split.by]] <- factor(data[[split.by]])
      }
      p <- p + ggplot2::facet_wrap(as.formula(paste("~", split.by)),
                          ncol = if (length(group.by) > 1) NULL else ncol)
    }

    return(p)
  })

  # Return results
  if (length(plots) == 1) {
    return(plots[[1]])
  } else if (combine) {
    ensure_pkgs("patchwork")
    return(patchwork::wrap_plots(plots, ncol = ncol))
  } else {
    return(plots)
  }
}

#' @rdname plot_umap
#' @export
plot_umap.Sap <- function(x,
                          segment_type = c("motifs", "syllables", "bouts", "segments"),
                          dims = c("UMAP1", "UMAP2"),
                          group.by = NULL,
                          split.by = NULL,
                          subset.by = NULL,
                          subset.value = NULL,
                          cols = NULL,
                          pt.size = 0.5,
                          stroke = 0.5,
                          alpha = 0.3,
                          highlight.alpha = NULL,
                          label = FALSE,
                          label.size = 4,
                          repel = FALSE,
                          highlight.by = NULL,
                          highlight.value = NULL,
                          cols.highlight = '#DE2D26',
                          sizes.highlight = 1,
                          background.value = NULL,
                          na.value = 'grey80',
                          ncol = NULL,
                          combine = TRUE,
                          verbose = TRUE,
                          ...) {
  if(verbose) message(sprintf("\n=== Starting UMAP Plotting for %s ===\n",
                              segment_type[1]))

  # Validate input
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  # Match segment_type argument
  segment_type <- match.arg(segment_type)

  # Get feature type
  feature_type <- sub("s$", "", segment_type)

  # Check if feat.embeds exists
  if (is.null(x$features[[feature_type]][["feat.embeds"]])) {
    stop(sprintf("No embeddings found for %s. Run run_umap() first.", segment_type))
  }

  # Get the data
  data <- x$features[[feature_type]][["feat.embeds"]]

  # Check if UMAP coordinates exist
  if (!all(dims %in% colnames(data))) {
    stop("UMAP coordinates not found. Run run_umap() first.")
  }

  # Call default method with all parameters
  p <- plot_umap.default(
    x = data,
    dims = dims,
    group.by = group.by,
    split.by = split.by,
    subset.by = subset.by,
    subset.value = subset.value,
    cols = cols,
    pt.size = pt.size,
    stroke = stroke,
    alpha = alpha,
    highlight.alpha = highlight.alpha,
    label = label,
    label.size = label.size,
    repel = repel,
    highlight.by = highlight.by,
    highlight.value = highlight.value,
    cols.highlight = cols.highlight,
    sizes.highlight = sizes.highlight,
    background.value = background.value,
    na.value = na.value,
    ncol = ncol,
    combine = combine,
    ...
  )

  print(p)

  invisible(x)
}

#' Internal UMAP Plotting Function
#'
#' @description
#' Internal function for creating individual UMAP plots.
#'
#' @keywords internal
plot_single_umap <- function(
    data,
    dims = c("X1", "X2"),
    cols = NULL,               # Cluster colors
    pt.size = 0.5,
    stroke = 0.5,
    alpha = 1,
    highlight.alpha = NULL,
    label = FALSE,
    repel = FALSE,
    label.size = 4,
    highlight.by = NULL,
    highlight.value = NULL,
    cols.highlight = '#DE2D26',
    sizes.highlight = 1,
    background.value = NULL,
    na.value = 'grey80'
) {
  # Input validation
  if (!all(dims %in% colnames(data))) {
    stop("Cannot find dimensions to plot in data")
  }

  # Ensure 'cluster' is a factor with levels sorted numerically
  data$cluster <- factor(as.character(data$cluster), levels = sort(as.numeric(unique(as.character(data$cluster)))))

  # Initialize 'color_group' variable in data
  data$color_group <- NA_character_

  # Validate highlighting parameters
  if (!is.null(highlight.by) && !is.null(highlight.value)) {
    if (!highlight.by %in% colnames(data)) {
      stop(paste("Cannot find", highlight.by, "in data"))
    }

    # Assign 'color_group' based on highlight and background values
    if (is.null(background.value)) {
      # If background.value is NULL, non-highlighted points are assigned to "Other"
      data$color_group <- ifelse(
        data[[highlight.by]] %in% highlight.value,
        as.character(data[[highlight.by]]),
        "Other"
      )
    } else {
      # Assign 'color_group' with specified background.value
      data$color_group <- ifelse(
        data[[highlight.by]] %in% highlight.value,
        as.character(data[[highlight.by]]),
        ifelse(
          data[[highlight.by]] == background.value,
          as.character(data$cluster),
          "Other"
        )
      )
    }

    # Prepare color mapping
    # Combine cluster colors and highlight colors
    if (length(cols.highlight) == 1) {
      highlight_colors <- setNames(rep(cols.highlight, length(highlight.value)), highlight.value)
    } else {
      highlight_colors <- setNames(cols.highlight, highlight.value)
    }

    if (is.null(background.value)) {
      # Only include highlight colors and "Other"
      color_mapping <- c(highlight_colors, "Other" = na.value)
    } else {
      # Include cluster colors, highlight colors, and "Other"
      color_mapping <- c(cols, highlight_colors, "Other" = na.value)
    }
    # Ensure no duplicated names
    color_mapping <- color_mapping[!duplicated(names(color_mapping))]

  } else {
    # If no highlighting, use cluster colors
    data$color_group <- data$cluster
    color_mapping <- cols
  }

  # Create base plot
  ensure_pkgs("ggplot2")
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[dims[1]]], y = .data[[dims[2]]], color = .data$color_group))

  # Adjust point sizes and alpha
  if (!is.null(highlight.by) && !is.null(highlight.value)) {
    data$is_highlight <- data[[highlight.by]] %in% highlight.value
    # Set alpha for highlighted and non-highlighted points
    point_alpha <- ifelse(data$is_highlight,
                          ifelse(is.null(highlight.alpha), alpha, highlight.alpha),
                          alpha)
    # Set point sizes
    point_size <- ifelse(data$is_highlight, sizes.highlight * pt.size, pt.size)
  } else {
    # No highlighting
    point_alpha <- alpha
    point_size <- pt.size
  }

  # Add points without stroke
  plot <- plot + ggplot2::geom_point(
    size = point_size,
    alpha = point_alpha,
    stroke = stroke
  )

  # Determine the legend name
  legend_name <- if (!is.null(highlight.by)) {
    highlight.by
  } else {
    "Cluster"
  }

  # Adjust the factor levels of 'color_group' to ensure correct legend order
  if (!is.null(highlight.by) && !is.null(highlight.value)) {
    if (is.null(background.value)) {
      # When background.value is NULL, legend includes highlight values and "Other"
      levels_order <- c(highlight.value, "Other")
      data$color_group <- factor(data$color_group, levels = levels_order)
    } else {
      # When background.value is provided, include clusters
      levels_order <- c(as.character(sort(as.numeric(levels(data$cluster)))), highlight.value, "Other")
      data$color_group <- factor(data$color_group, levels = levels_order)
    }
  } else {
    # No highlighting, cluster levels
    data$color_group <- factor(data$color_group, levels = as.character(sort(as.numeric(levels(data$cluster)))))
  }

  # Add color scale
  plot <- plot + ggplot2::scale_color_manual(
    values = color_mapping,
    na.value = na.value,
    name = legend_name
  )

  # Add labels for clusters if needed
  if (label) {
    # Only label clusters
    label_data <- aggregate(
      data[, dims],
      by = list(cluster = data$cluster),
      FUN = mean
    )
    colnames(label_data)[-1] <- dims

    if (repel) {
      plot <- plot + ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes_string(x = dims[1], y = dims[2], label = "cluster"),
        size = label.size,
        fontface = "bold",
        inherit.aes = FALSE
      )
    } else {
      plot <- plot + ggplot2::geom_text(
        data = label_data,
        ggplot2::aes_string(x = dims[1], y = dims[2], label = "cluster"),
        size = label.size,
        fontface = "bold",
        inherit.aes = FALSE
      )
    }
  }

  # Add theme and styling
  plot <- plot +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_blank()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))

  return(plot)
}

#' Plot UMAP Visualization for Trajectory Analysis
#'
#' @description
#' Creates UMAP visualizations optimized for trajectory analysis, with support
#' for continuous color mapping and overlay comparisons.
#'
#' @param x An object to visualize, either a data frame or SAP object
#' @param dims UMAP dimensions to plot (default: c("UMAP1", "UMAP2"))
#' @param color.by Column for continuous color mapping (default: ".time")
#' @param split.by Column for faceting (default: "label")
#' @param order.by Column for ordering facets (default: "day_post_hatch")
#' @param pt.size Point size (default: 1.2)
#' @param alpha_range Range for alpha transparency (default: c(0.1, 0.5))
#' @param ncol Number of columns in layout
#' @param title Plot title
#' @param overlay_mode Whether to create overlay comparisons (default: FALSE)
#' @param base_label Base label for overlay comparison
#' @param compare_labels Labels to compare against base
#' @param base_color Color for base label (default: "steelblue")
#' @param compare_color Color for comparison labels (default: "orangered")
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables', 'bouts', 'segments')
#' @param data_type For SAP objects: Type of embedding data ('feat.embeds', 'traj.embeds')
#' @param verbose For SAP objects: Whether to print progress messages
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' Supports two visualization modes:
#' \itemize{
#'   \item Standard mode: Continuous color mapping for trajectory visualization
#'   \item Overlay mode: Direct comparison between trajectory patterns
#' }
#'
#' @return
#' For default method: A ggplot object
#' For SAP objects: Updated SAP object with plot as side effect
#'
#' @examples
#' \dontrun{
#' # Basic trajectory plot
#' plot_umap2(traj_df, color.by = ".time")
#'
#' # Overlay comparison plot
#' plot_umap2(traj_df,
#'            overlay_mode = TRUE,
#'            base_label = "a",
#'            compare_labels = c("b", "c"))
#'
#' # Plot motif trajectories from SAP object
#' plot_umap2(sap_obj,
#'            segment_type = "motifs",
#'            data_type = "traj.embeds",
#'            color.by = ".time")
#'
#' # Compare trajectories between labels
#' plot_umap2(sap_obj,
#'            segment_type = "motifs",
#'            data_type = "traj.embeds",
#'            overlay_mode = TRUE,
#'            base_label = "pre")
#' }
#'
#' @rdname plot_umap2
#' @export
plot_umap2 <- function(x, ...) {
  UseMethod("plot_umap2")
}

#' @rdname plot_umap2
#' @export
plot_umap2.default <- function(x,
                               dims = c("UMAP1", "UMAP2"),
                               color.by = ".time",
                               split.by = "label",
                               order.by = "day_post_hatch",
                               pt.size = 1.2,
                               alpha_range = c(0.1, 0.5),
                               ncol = NULL,
                               title = NULL,
                               overlay_mode = FALSE,
                               base_label = NULL,
                               compare_labels = NULL,
                               base_color = "steelblue",
                               compare_color = "orangered",
                               ...) {

  # Validate input
  if(!is.data.frame(x)) stop("Input must be a data frame")

  if(overlay_mode) {
    # Validate overlay parameters
    if(is.null(base_label)) stop("base_label must be provided in overlay mode")
    if(!base_label %in% x[[split.by]]) stop("base_label not found in data")

    # If compare_labels is NULL, use all labels except base_label
    if(is.null(compare_labels)) {
      compare_labels <- unique(x[[split.by]])[unique(x[[split.by]]) != base_label]
    }

    # Auto-determine number of columns if not specified
    if(is.null(ncol)) {
      n_plots <- length(compare_labels)
      # Square root heuristic for layout
      ncol <- ceiling(sqrt(n_plots))
    }

    # Create overlay plots
    ensure_pkgs("ggplot2")
    all_plots <- lapply(compare_labels, function(compare_label) {
      filtered_data <- x |>
        dplyr::filter(!!rlang::sym(split.by) %in% c(base_label, compare_label))

      # Create color mapping
      color_values <- c(base_color, compare_color)
      names(color_values) <- c(base_label, compare_label)

      ggplot2::ggplot(filtered_data,
                      ggplot2::aes(!!rlang::sym(dims[1]),
                                   !!rlang::sym(dims[2]))) +
        ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(split.by),
                                         alpha = !!rlang::sym(color.by)),
                            size = pt.size,
                            stroke = 0) +
        ggplot2::scale_color_manual(values = color_values) +
        ggplot2::scale_alpha_continuous(range = alpha_range) +
        ggplot2::ggtitle(paste(title %||% "UMAP Overlay of",
                               base_label, "vs", compare_label)) +
        ggplot2::theme_minimal()
    })
    ensure_pkgs("patchwork")
    return(patchwork::wrap_plots(all_plots, ncol = ncol))

  } else {
    # Original functionality for continuous color mapping
    if(!is.null(split.by) && !is.null(order.by)) {
      if(!all(c(split.by, order.by) %in% colnames(x))) {
        stop(sprintf("Cannot find %s or %s in data", split.by, order.by))
      }

      # Calculate order based on mean values
      label_order <- x |>
        dplyr::group_by(!!rlang::sym(split.by)) |>
        dplyr::summarise(mean_val = mean(!!rlang::sym(order.by), na.rm = TRUE)) |>
        dplyr::arrange(.data$mean_val) |>
        dplyr::pull(!!rlang::sym(split.by))

      x[[split.by]] <- factor(x[[split.by]], levels = label_order)
    }

    # Create regular plot
    p <- ggplot2::ggplot(x, ggplot2::aes(!!rlang::sym(dims[1]),
                                         !!rlang::sym(dims[2]))) +
      ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(color.by)),
                          size = pt.size,
                          alpha = alpha_range[2],
                          stroke = 0) +
      ggplot2::scale_color_viridis_c(option = "inferno") +
      ggplot2::theme_minimal()

    if (!is.null(split.by)) {
      p <- p + ggplot2::facet_wrap(as.formula(paste("~", split.by)),
                                   ncol = ncol)
    }

    if (!is.null(title)) {
      p <- p + ggplot2::ggtitle(title)
    }

    return(p)
  }
}

#' @rdname plot_umap2
#' @export
plot_umap2.Sap <- function(x,
                           segment_type = c("motifs", "syllables", "bouts", "segments"),
                           data_type = c("feat.embeds", "traj.embeds"),
                           dims = c("UMAP1", "UMAP2"),
                           color.by = ".time",
                           split.by = "label",
                           order.by = "day_post_hatch",
                           pt.size = 1.2,
                           alpha_range = c(0.1, 0.5),
                           ncol = NULL,
                           title = NULL,
                           overlay_mode = FALSE,
                           base_label = NULL,
                           compare_labels = NULL,
                           base_color = "steelblue",
                           compare_color = "orangered",
                           verbose = TRUE,
                           ...) {

  if(verbose) message(sprintf("\n=== Starting UMAP Plotting for %s using %s ===\n",
                              segment_type[1],
                              data_type[1]))

  # Validate input
  if (!inherits(x, "Sap")) stop("Input must be a SAP object")

  # Match arguments
  segment_type <- match.arg(segment_type)
  data_type <- match.arg(data_type)

  # Get feature type
  feature_type <- sub("s$", "", segment_type)

  # Get appropriate data frame
  plot_data <- x$features[[feature_type]][[data_type]]
  if(is.null(plot_data)) {
    stop(sprintf("No %s data found in %s features", data_type, feature_type))
  }

  # Check if UMAP coordinates exist
  if(!all(dims %in% colnames(plot_data))) {
    stop(sprintf("UMAP coordinates (%s) not found. Run UMAP first.",
                 paste(dims, collapse=", ")))
  }

  # Create title if not provided
  if(is.null(title)) {
    title <- sprintf("UMAP visualization of %s %s", segment_type, data_type)
  }


  # Create plot using default method
  p <- plot_umap2.default(
    x = plot_data,
    dims = dims,
    color.by = color.by,
    split.by = split.by,
    order.by = order.by,
    pt.size = pt.size,
    alpha_range = alpha_range,
    ncol = ncol,
    title = title,
    overlay_mode = overlay_mode,
    base_label = base_label,
    compare_labels = compare_labels,
    base_color = base_color,
    compare_color = compare_color,
    ...
  )

  # Print the plot
  print(p)

  invisible(x)
}
