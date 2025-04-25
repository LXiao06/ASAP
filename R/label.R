
# Labeling syllables in song motif ----------------------------------------
# Update date: Apr. 16, 2025

#' Automatic Syllable Labeling in Song Motif
#'
#' @description
#' Performs two-stage clustering of bird song segments to identify syllables:
#' 1. Temporal clustering using weighted time and duration features
#' 2. UMAP-based refinement of temporal clusters
#' Finally merges similar clusters based on UMAP proximity.
#'
#' @param x A Sap object containing bird song analysis data
#' @param eps_time Epsilon parameter for temporal DBSCAN clustering (default: 0.1)
#' @param eps_umap Epsilon parameter for UMAP-based DBSCAN clustering (default: 0.6)
#' @param min_pts Minimum points parameter for DBSCAN clustering (default: 5)
#' @param outlier_threshold Threshold for removing small clusters (default: 0.01)
#' @param umap_threshold Distance threshold for merging similar clusters (default: 1)
#' @param weight_time Weight factor for temporal features (default: 4)
#' @param ... Additional arguments (not currently used)
#'
#' @details
#' The function performs clustering in multiple stages:
#'
#' 1. Temporal Clustering:
#'    * Combines time position and duration information
#'    * Uses weighted features (controlled by weight_time):
#'      - weight_time > 1: Emphasizes temporal position over duration
#'      - weight_time = 1: Equal weighting
#'      - weight_time < 1: Emphasizes duration over temporal position
#'
#' 2. UMAP-based Refinement:
#'    * Further splits temporal clusters based on UMAP coordinates
#'    * Helps distinguish syllables with similar timing but different acoustic features
#'
#' 3. Cluster Cleaning:
#'    * Removes small clusters (controlled by outlier_threshold)
#'    * Merges similar clusters based on UMAP proximity (controlled by umap_threshold)
#'
#' Parameter Tuning Guidelines:
#' * eps_time: Controls temporal separation sensitivity
#'   - Smaller values create more temporal splits
#'   - Typical range: 0.05-0.2
#' * eps_umap: Controls acoustic feature sensitivity
#'   - Smaller values create more acoustic splits
#'   - Typical range: 0.4-0.8
#' * weight_time: Controls temporal vs duration importance
#'   - Default (4) weights time 4x more than duration
#'   - Increase for more temporal separation
#'   - Decrease for more duration-based separation
#'
#' @return
#' Returns the input Sap object with updated syllable clusters in:
#' x$features$syllable$feat.embeds
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' sap <- auto_label(sap)
#'
#' # Emphasize temporal separation
#' sap <- auto_label(sap,
#'                  eps_time = 0.08,
#'                  weight_time = 6)
#'
#' # More acoustic feature sensitivity
#' sap <- auto_label(sap,
#'                  eps_umap = 0.4,
#'                  min_pts = 3)
#'
#' # Stricter cluster cleaning
#' sap <- auto_label(sap,
#'                  outlier_threshold = 0.02,
#'                  umap_threshold = 0.8)
#' }
#'
#' @seealso
#' * manual_label() for manual refinement of automatic clusters
#' * plot_cluster() for visualizing clustering results
#'
#' @export
auto_label <- function(x,
                       eps_time = 0.1,
                       eps_umap = 0.6,
                       min_pts = 5,
                       outlier_threshold = 0.01,
                       umap_threshold = 1,
                       weight_time = 4,
                       ...) {
  # Check if dbscan is installed and available
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required for auto_label function.\n",
         "Please install it using: install.packages('dbscan')")
  }

  # Check if input is a Sap object
  if (!inherits(x, "Sap")) {
    stop("Input must be a Sap object")
  }

  # Validation
  if (is.null(x$motifs)) stop("motifs data required for auto_label")
  if (is.null(x$features$segment$feat.embeds)) stop("Feature embeddings for segments required for auto_label")

  # Prepare df for clustering
  motifs <- x$motifs |>
    dplyr::rename(motif_start = start_time, motif_end = end_time)|>
    dplyr::group_by(filename) |>
    dplyr::mutate(motif_index = dplyr::row_number()) |>
    dplyr::ungroup()

  data <- x$features$segment$feat.embeds |>
    dplyr::inner_join(motifs, by = c("filename", "day_post_hatch", "label"),
                      relationship = "many-to-many") |>
    dplyr::filter(start_time >= .data$motif_start & start_time <= .data$motif_end) |>
    dplyr::mutate(
      segment_duration = end_time - start_time,
      # motif_duration = motif_end - motif_start,
      # .start = start_time - motif_start,
      # .end = end_time - motif_start,
      .time =(start_time + end_time)/2 - .data$motif_start
      )


  # Initialize list to store temporary results
  results <- list()

  # Process individual label
  for(label in unique(data$label)) {
    # Subset data for current label
    label_data <- data[data$label == label, ]

    # Calculate total points for this label
    total_points <- nrow(label_data)

    # Scale features within this label
    scaled_features <- data.frame(
      UMAP1 = scale(label_data$UMAP1),
      UMAP2 = scale(label_data$UMAP2),
      time = scale(label_data$.time),
      duration = scale(label_data$segment_duration)
    )

    # Calculate normalized weights
    time_weight <- weight_time / (weight_time + 1)
    duration_weight <- 1 / (weight_time + 1)

    # Apply weights to temporal features
    weighted_features <- cbind(
      scaled_features$time * time_weight,
      scaled_features$duration * duration_weight
    )

    # First stage: Time-based clustering
    time_clusters <- dbscan::dbscan(
      weighted_features,
      eps = eps_time,
      minPts = min_pts
    )$cluster

    # Process each temporal cluster
    refined_clusters <- numeric(nrow(label_data))
    cluster_counter <- 0

    for(t_cluster in unique(time_clusters)) {
      if(t_cluster == -1) next

      temp_idx <- which(time_clusters == t_cluster)
      temp_points <- scaled_features[temp_idx, c("UMAP1", "UMAP2")]

      # Second stage: UMAP-based clustering
      umap_clusters <- dbscan::dbscan(
        temp_points,
        eps = eps_umap,
        minPts = min_pts
      )$cluster

      valid_clusters <- umap_clusters != -1
      refined_clusters[temp_idx[valid_clusters]] <-
        umap_clusters[valid_clusters] + cluster_counter

      cluster_counter <- max(refined_clusters, na.rm = TRUE) + 1
    }

    # Calculate cluster statistics
    label_results <- data.frame(
      label_data,
      cleaned_cluster = refined_clusters
    )

    # Remove small clusters
    cluster_sizes <- table(refined_clusters)
    small_clusters <- names(cluster_sizes)[cluster_sizes/total_points < outlier_threshold]

    # Store motif IDs to remove
    motifs_to_remove <- unique(label_results$filename[label_results$cleaned_cluster %in% small_clusters])

    # Remove entire motifs that contain small clusters
    label_results <- label_results[!label_results$filename %in% motifs_to_remove, ]

    # Store results
    results[[label]] <- label_results
  }

  # Merge all labels into single dataframe
  cleaned_data <- do.call(rbind, results)

  # Merge similar clusters based on UMAP embeddings
  merged_clusters <- merge_similar_clusters(cleaned_data, umap_threshold = umap_threshold) |>
    dplyr::select(filename, day_post_hatch, label, start_time, end_time, .data$merged_cluster, UMAP1, UMAP2) |>
    dplyr::rename(cluster = .data$merged_cluster)

  # Upadte syllable data
  x[["features"]][["syllable"]][["feat.embeds"]] <- merged_clusters

  invisible(x)
}

#' Merge Similar Clusters Based on UMAP Proximity
#'
#' @description
#' Internal function to merge clusters that are close in UMAP space.
#'
#' @param data Data frame containing cluster information
#' @param umap_threshold Distance threshold for merging clusters
#'
#' @return Data frame with merged cluster assignments
#'
#' @keywords internal
merge_similar_clusters <- function(data, umap_threshold = NULL) {
  # Calculate cluster centers in UMAP space
  cluster_centers <- data |>
    dplyr::group_by(.data$label, .data$cleaned_cluster) |>
    dplyr::reframe(
      mean_UMAP1 = mean(UMAP1),
      mean_UMAP2 = mean(UMAP2),
      n_points = dplyr::n(),
      cluster_id = paste(.data$label, .data$cleaned_cluster, sep = "_")
    ) |>
    dplyr::distinct(.data$cluster_id, .keep_all = TRUE)

  # Create initial cluster ID mapping
  data$cluster_id <- paste(data$label, data$cleaned_cluster, sep = "_")

  # Calculate pairwise distances using matrix operations
  umap_coords <- as.matrix(cluster_centers[, c("mean_UMAP1", "mean_UMAP2")])
  dist_matrix <- as.matrix(dist(umap_coords))

  # Create adjacency matrix based on threshold
  adj_matrix <- dist_matrix < umap_threshold
  diag(adj_matrix) <- FALSE  # Remove self-connections

  # Convert to edge list format
  edges <- which(adj_matrix, arr.ind = TRUE) |>
    as.data.frame() |>
    dplyr::filter(row > col) |>  # Keep only unique pairs
    dplyr::mutate(
      from_id = cluster_centers$cluster_id[row],
      to_id = cluster_centers$cluster_id[col]
    )

  if(nrow(edges) > 0) {
    # Create connected components
    n_clusters <- nrow(cluster_centers)
    merged_clusters <- seq_len(n_clusters)
    names(merged_clusters) <- cluster_centers$cluster_id

    # Merge connected components
    for(i in seq_len(nrow(edges))) {
      from_cluster <- merged_clusters[edges$from_id[i]]
      to_cluster <- merged_clusters[edges$to_id[i]]
      if(!is.na(from_cluster) && !is.na(to_cluster)) {
        min_cluster <- min(from_cluster, to_cluster)
        merged_clusters[merged_clusters == from_cluster |
                          merged_clusters == to_cluster] <- min_cluster
      }
    }

    # Make cluster numbers continuous
    unique_merged <- sort(unique(merged_clusters))
    cluster_mapping <- data.frame(
      old_cluster = unique_merged,
      new_cluster = seq_along(unique_merged)
    )

    # Create final mapping
    result <- data |>
      dplyr::left_join(
        data.frame(
          cluster_id = names(merged_clusters),
          old_cluster = merged_clusters
        ),
        by = "cluster_id"
      ) |>
      dplyr::left_join(cluster_mapping, by = "old_cluster") |>
      dplyr::mutate(merged_cluster = .data$new_cluster) |>
      dplyr::select(all_of(c("old_cluster", "new_cluster")))

  } else {
    # If no edges, create continuous cluster numbers from 1 to n
    result <- data |>
      dplyr::mutate(merged_cluster = .data$cleaned_cluster)
  }

  return(result)
}

#' Manual Syllable Labeling for Bird Song Analysis
#'
#' @description
#' Provides interactive or map-based assignment of alphabetic labels to automatically
#' clustered syllables. Supports both interactive console input and predefined label
#' maps, with persistent storage of label assignments.
#'
#' @param x A Sap object containing bird song analysis data
#' @param data_type Type of data to label: "segment" or "syllable"
#' @param label_map Optional data frame for non-interactive labeling with columns:
#'   * cluster: numeric cluster IDs
#'   * syllable: corresponding letter labels (a-z)
#' @param interactive Logical; force interactive mode even if label_map exists
#' @param verbose Logical; print progress messages and summaries
#'
#' @details
#' The function provides two main methods for syllable labeling:
#'
#' 1. Interactive Console Mode:
#'    * Prompts for letter input for each cluster
#'    * Validates input (single letters a-z only)
#'    * Allows quitting mid-process ('q')
#'    * Automatically stores label map for future use
#'
#' 2. Label Map Mode:
#'    * Uses predefined mapping of clusters to letters
#'    * Validates completeness and format of mapping
#'    * Updates stored label map
#'
#' Label Storage and Retrieval:
#'    * Label maps are stored in the Sap object attributes
#'    * Automatically retrieved in subsequent runs
#'    * Can be overridden with new map or interactive session
#'    * Accessible via attr(x$features$<data_type>$feat.embeds, "label_map")
#'
#' @return
#' Returns modified Sap object with:
#'   * Updated syllable labels in x$syllables
#'   * Stored label map in feature embeddings attributes
#'   * All labels converted to lowercase
#'
#' @examples
#' \dontrun{
#' # Interactive labeling (first time)
#' sap <- manual_label(sap, data_type = "syllable")
#'
#' # Using stored labels (subsequent runs)
#' sap <- manual_label(sap, data_type = "syllable")
#'
#' # Force new interactive session
#' sap <- manual_label(sap,
#'                    data_type = "syllable",
#'                    interactive = TRUE)
#'
#' # Using predefined label map
#' label_map <- data.frame(
#'   cluster = 1:5,
#'   syllable = c("a", "b", "b", "c", "a")
#' )
#' sap <- manual_label(sap,
#'                    data_type = "syllable",
#'                    label_map = label_map)
#'
#' # Access stored label map
#' stored_map <- attr(sap$features$syllable$feat.embeds, "label_map")
#' print(stored_map)
#'
#' # Typical workflow:
#' # 1. Run automatic clustering
#' sap <- auto_label(sap)
#'
#' # 2. View clusters
#' plot_cluster(sap,
#'             data_type = "syllable",
#'             label_type = "pre")
#'
#' # 3. Manually label clusters
#' sap <- manual_label(sap, data_type = "syllable")
#'
#' # 4. View syllable labels
#' plot_cluster(sap,
#'             data_type = "syllable",
#'             label_type = "post")
#' }
#'
#' @section Label Map Format:
#' The label_map data frame must have:
#' ```r
#' data.frame(
#'   cluster = numeric_vector,  # Cluster IDs from auto_label
#'   syllable = character_vector  # Single letters a-z
#' )
#' ```
#'
#' @section Tips for Interactive Labeling:
#' * View clusters first using plot_cluster(..., label_type = "pre")
#' * Consider similar acoustic features when assigning same letter
#' * Use 'q' to exit without saving if mistakes are made
#' * Review label summary after completion
#' * Store label_map externally for reproducibility
#'
#' @section Common Workflows:
#' 1. First-time labeling:
#'    ```r
#'    sap <- auto_label(sap)  # Generate clusters
#'    sap <- manual_label(sap, data_type = "syllable")  # Label interactively
#'    ```
#'
#' 2. Using stored labels:
#'    ```r
#'    sap <- manual_label(sap, data_type = "syllable")  # Uses stored map
#'    ```
#'
#' 3. Updating labels:
#'    ```r
#'    # Extract and modify existing map
#'    map <- attr(sap$features$syllable$feat.embeds, "label_map")
#'    map$syllable[map$cluster == 3] <- "d"  # Change cluster 3 to 'd'
#'    sap <- manual_label(sap, data_type = "syllable", label_map = map)
#'    ```
#'
#' @seealso
#' * auto_label() for automatic cluster generation
#' * plot_cluster() for visualizing clusters and syllables
#'
#' @export
manual_label <- function(x,
                         data_type = c("segment", "syllable"),
                         label_map = NULL,
                         interactive = FALSE,
                         verbose = TRUE) {
  # Input validation
  data_type <- match.arg(data_type)
  if (!inherits(x, "Sap")) stop("Input must be a Sap object")

  # Select appropriate data based on data_type
  feat_data <- switch(data_type,
                      "segment" = x[["features"]][["segment"]][["feat.embeds"]],
                      "syllable" = x[["features"]][["syllable"]][["feat.embeds"]])

  if(is.null(feat_data)) {
    stop(sprintf("No feature embeddings found for %s", data_type))
  }

  if(!"cluster" %in% names(feat_data)) {
    stop("No cluster column found in the data")
  }

  # Get unique clusters
  clusters <- sort(unique(feat_data$cluster))
  if(verbose) {
    message("Found clusters: ", paste(clusters, collapse = ", "))
  }

  # Check for existing label map in attributes
  stored_map <- attr(feat_data, "label_map")
  if(!is.null(stored_map) && is.null(label_map) && !interactive) {
    if(verbose) {
      message("\nFound stored label map in: x$features$", data_type, "$feat.embeds")
      message("Access it using: attr(x$features$", data_type, "$feat.embeds, 'label_map')")
      message("Using stored label map:")
      print(stored_map)
      message("\nTo override, provide new label_map or set interactive = TRUE")
    }
    label_map <- stored_map
  }

  # Determine whether to use interactive mode
  use_interactive <- interactive || (is.null(label_map) && is.null(stored_map))

  if(use_interactive) {
    if(verbose) {
      message("\nEntering interactive labeling mode...")
      message("Please enter a letter (a-z) for each cluster.")
      message("Press Enter after each input, or 'q' to quit.")
    }

    # Interactive labeling
    labels <- character(length(clusters))
    names(labels) <- clusters

    for(i in seq_along(clusters)) {
      repeat {
        if(verbose) message(sprintf("\nCluster %d:", clusters[i]))
        input <- readline(prompt = "Enter letter (a-z): ")

        if(input == "q") {
          return(x)  # Exit without saving
        }

        if(nchar(input) == 1 && grepl("^[a-zA-Z]$", input)) {
          labels[i] <- tolower(input)
          break
        } else {
          message("Invalid input. Please enter a single letter.")
        }
      }
    }

    # Create label_map from interactive input
    label_map <- data.frame(
      cluster = clusters,
      syllable = labels
    )

    # Store label_map in attributes
    if(verbose) {
      message("\nStoring label map in: x$features$", data_type, "$feat.embeds")
      message("Access it using: attr(x$features$", data_type, "$feat.embeds, 'label_map')")
    }

  } else {
    # Validate label_map
    if(!all(c("cluster", "syllable") %in% names(label_map))) {
      stop("label_map must contain 'cluster' and 'syllable' columns")
    }

    if(!all(clusters %in% label_map$cluster)) {
      stop("label_map must contain all clusters: ",
           paste(setdiff(clusters, label_map$cluster), collapse = ", "))
    }

    # Validate syllable labels
    invalid_labels <- !grepl("^[a-zA-Z]$", label_map$syllable)
    if(any(invalid_labels)) {
      stop("Invalid syllable labels found: ",
           paste(unique(label_map$syllable[invalid_labels]), collapse = ", "))
    }
  }

  # Update data with new labels
  label_map$syllable <- tolower(label_map$syllable)
  feat_data$syllable <- label_map$syllable[match(feat_data$cluster, label_map$cluster)]

  # Store the label_map in attributes
  attr(feat_data, "label_map") <- label_map

  # Update Sap object with new feat_data (including attributes)
  x[["features"]][[data_type]][["feat.embeds"]] <- feat_data
  x[["syllables"]] <- feat_data

  if(verbose) {
    message("\nLabeling complete. Summary:")
    print(table(feat_data$syllable))
  }

  invisible(x)
}
