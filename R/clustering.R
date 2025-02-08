# Clustering
# Update date : Feb. 7, 2025

# Find Clusters -----------------------------------------------------------


#' Find Clusters in Feature Data
#'
#' @description
#' Performs cluster analysis on feature data using shared nearest neighbor (SNN) clustering.
#'
#' @param x An object to analyze, either a data frame or SAP object
#' @param metadata_cols For default method: Column indices for metadata
#' @param k.param Number of nearest neighbors (default: 20)
#' @param prune.SNN Pruning threshold for SNN graph (default: 1/15)
#' @param n.pcs Number of principal components to use (default: 20)
#' @param resolution Resolution parameter for clustering (default: 0.2)
#' @param n.start Number of random starts (default: 10)
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables', 'bouts', 'segments')
#' @param data_type For SAP objects: Type of feature data ('spectral_feature')
#' @param label For SAP objects: Specific label to filter data
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' For feature data frames:
#' \itemize{
#'   \item Separates metadata and feature columns
#'   \item Finds nearest neighbors using PCA
#'   \item Constructs SNN graph
#'   \item Performs community detection
#' }
#'
#' For SAP objects:
#' \itemize{
#'   \item Supports multiple segment types
#'   \item Optional label filtering
#'   \item Stores results in features slot
#'   \item Updates feature embeddings
#' }
#'
#' The clustering approach is similar to that used in Seurat V3, implementing:
#' \itemize{
#'   \item PCA-based neighbor finding
#'   \item SNN graph construction
#'   \item Louvain community detection
#' }
#'
#' @return
#' For default method: Data frame with metadata columns and cluster assignments
#' For SAP objects: Updated object with clustering results in features slot
#'
#' @examples
#' \dontrun{
#' # Basic clustering of feature data
#' clusters <- find_clusters(features,
#'                          metadata_cols = c(1:5))
#'
#' # Clustering with custom parameters
#' clusters <- find_clusters(features,
#'                          metadata_cols = c(1:5),
#'                          k.param = 30,
#'                          resolution = 0.3)
#'
#' # Cluster SAP object features
#' sap_obj <- find_clusters(sap_object,
#'                         segment_type = "motifs",
#'                         data_type = "spectral_feature")
#'
#' # Label-specific clustering
#' sap_obj <- find_clusters(sap_object,
#'                         segment_type = "syllables",
#'                         label = "a",
#'                         resolution = 0.4)
#' }
#'
#' @seealso \code{\link{run_umap}} for visualization of clusters
#'
#' @rdname find_clusters
#' @export
find_clusters <- function(x, ...) {
  UseMethod("find_clusters")
}

#' @rdname find_clusters
#' @export
find_clusters.default <- function(x,
                                  metadata_cols,
                                  k.param = 20,
                                  prune.SNN = 1/15,
                                  n.pcs = 20,
                                  resolution = 0.2,
                                  n.start = 10,
                                  verbose = TRUE,
                                  ...) {

  # Validate input
  if(!is.data.frame(x)) stop("Input must be a data frame")
  if(!all(metadata_cols %in% seq_len(ncol(x)))) {
    stop("Invalid metadata_cols indices")
  }

  # Print column information
  if(verbose) {
    message("Metadata columns:")
    print(colnames(x)[metadata_cols])
    message("\nFeature columns:")
    print(colnames(x)[-metadata_cols])
  }

  # Extract features
  features <- x[,-metadata_cols]

  # Find neighbors
  message("\nFinding neighbors...")
  neighbors <- find_segment_neighbors(
    segments_data = features,
    k.param = k.param,
    prune.SNN = prune.SNN,
    compute.SNN = TRUE,
    n.pcs = n.pcs
  )

  # Print neighbor results
  message("\nNeighbor graph dimensions:")
  print(dim(neighbors$snn))


  # Find clusters
  message("\nFinding clusters...")
  clusters <- find_segment_clusters(
    neighbor_graphs = neighbors,
    resolution = resolution,
    n.start = n.start
  )

  # Print cluster results
  message("\nCluster assignments:")
  print(table(clusters[,1]))

  # Combine results
  seg_cluster <- cbind(
    x[,metadata_cols, drop = FALSE],
    cluster = clusters[,1]
  )

  return(seg_cluster)
}

#' @rdname find_clusters
#' @export
find_clusters.Sap <- function(x,
                              segment_type = c("motifs", "syllables", "bouts", "segments"),
                              data_type = c("spectral_feature"),
                              label = NULL,
                              k.param = 20,
                              prune.SNN = 1/15,
                              n.pcs = 20,
                              resolution = 0.2,
                              n.start = 10,
                              verbose = TRUE,
                              ...) {

  if(verbose) message(sprintf("\n=== Starting cluster analysis for %s using %s ===\n",
                              segment_type[1],
                              data_type[1]))

  # Validate input
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }

  # Match segment_type argument
  segment_type <- match.arg(segment_type)
  #data_type <- match.arg(data_type)

  # Get feature data
  feature_type <- sub("s$", "", segment_type)  # Remove 's' from end
  feature_data <- x$features[[feature_type]][["spectral_feature"]]

  if (is.null(feature_data)) {
    stop(sprintf("%s features not found in SAP object", segment_type))
  }

  # Check if label filtering is requested
  if (!is.null(label)) {
    if (!"label" %in% names(feature_data)) {
      stop("Label column not found in feature data")
    }

    # Verify label exists in the data
    available_labels <- unique(feature_data$label)
    if (!label %in% available_labels) {
      stop(sprintf("Label '%s' not found. Available labels: %s",
                   label,
                   paste(available_labels, collapse = ", ")))
    }

    # Filter data by label
    feature_data <- feature_data[feature_data$label == label, ]

    if (verbose) {
      message(sprintf("\nFiltered data for label '%s': %d samples",
                      label, nrow(feature_data)))
    }
  }

  # Identify metadata columns
  standard_meta <- c("filename", "day_post_hatch", "label", "start_time", "end_time")
  if ("duration" %in% names(feature_data) && length(unique(feature_data$duration)) == 1) {
    metadata_cols <- which(names(feature_data) %in% c(standard_meta, "duration"))
  } else {
    metadata_cols <- which(names(feature_data) %in% standard_meta)
  }

  # Call default method
  result <- find_clusters.default(
    x = feature_data,
    metadata_cols = metadata_cols,
    k.param = k.param,
    prune.SNN = prune.SNN,
    n.pcs = n.pcs,
    resolution = resolution,
    n.start = n.start,
    verbose = verbose,
    ...
  )

  x$features[[feature_type]][["clusters"]] <- result
  if (verbose) {
    message("\nCluster assignments stored in clusters data frame")
    message("Access via: sap$features$", feature_type, "$clusters")
  }

  # Check if feat.embeds exists
  if (!is.null(x$features[[feature_type]][["feat.embeds"]])) {
    x$features[[feature_type]][["feat.embeds"]] <- NULL

    if (verbose) {
      message("\nPrevious feat.embeds removed")
      message("Run run_umap() again to create new feat.embeds with current clusters")
    }
  }

  # Store clustering parameters
  x$features[[feature_type]][["cluster_params"]] <- list(
    k.param = k.param,
    prune.SNN = prune.SNN,
    n.pcs = n.pcs,
    resolution = resolution,
    n.start = n.start
  )

  invisible(x)
}

#' Internal Functions from Seurat Package
#'
#' @name find_clusters
#' @description
#' Internal functions adapted from the Seurat package for neighbor finding
#' and cluster analysis. These functions are modified versions of the
#' original Seurat code.
#'
#' @details
#' Original source: Seurat package
#' (\url{https://github.com/satijalab/seurat/blob/main/inst/CITATION})
#'
#' These functions are based on methods developed in the following publications:
#' \itemize{
#'   \item Hao Y, et al. (2023). Dictionary learning for integrative, multimodal
#'         and scalable single-cell analysis. \emph{Nature Biotechnology}.
#'         doi:10.1038/s41587-023-01767-y
#'   \item Hao Y, Hao S, et al. (2021). Integrated analysis of multimodal
#'         single-cell data. \emph{Cell}.
#'         doi:10.1016/j.cell.2021.04.048
#'   \item Stuart T, Butler A, et al. (2019). Comprehensive Integration of
#'         Single-Cell Data. \emph{Cell}.
#'         doi:10.1016/j.cell.2019.05.031
#'   \item Butler A, et al. (2018). Integrating single-cell transcriptomic data
#'         across different conditions, technologies, and species.
#'         \emph{Nature Biotechnology}.
#'         doi:10.1038/nbt.4096
#'   \item Satija R, et al. (2015). Spatial reconstruction of single-cell gene
#'         expression data. \emph{Nature Biotechnology}.
#'         doi:10.1038/nbt.3192
#' }
#'
#' @keywords internal
NULL


#' Find Nearest Neighbors for Segments
#'
#' @description
#' Internal function to compute nearest neighbors and SNN graph.
#'
#' @param segments_data Matrix or data frame of segment features
#' @param k.param Number of nearest neighbors
#' @param prune.SNN SNN pruning threshold
#' @param compute.SNN Whether to compute SNN graph
#' @param n.pcs Number of PCs to use
#'
#' @return List containing neighbor graphs
#'
#' @keywords internal
find_segment_neighbors <- function(
    segments_data,
    k.param = 20,
    prune.SNN = 1/15,
    compute.SNN = TRUE,
    n.pcs = 10
) {

  # # Source C++ files
  # cpp_files <- c("snn.cpp", "RModularityOptimizer.cpp")
  # for(file in cpp_files) {
  #   cpp_path <- file.path("src", file)
  #   if(!file.exists(cpp_path)) {
  #     stop(paste("C++ file not found:", cpp_path))
  #   }
  #   Rcpp::sourceCpp(cpp_path)
  # }

  # Input validation
  if (!is.matrix(segments_data) && !is.data.frame(segments_data)) {
    stop("segments_data must be a matrix or data frame")
  }

  # Convert data frame to matrix if needed
  segments_data <- as.matrix(segments_data)

  # Ensure we have rownames
  if (is.null(rownames(segments_data))) {
    rownames(segments_data) <- paste0("segment_", 1:nrow(segments_data))
  }

  # Check number of segments vs k.param
  n.segments <- nrow(segments_data)
  if (n.segments < k.param) {
    warning("k.param set larger than number of segments. Setting k.param to number of segments - 1.")
    k.param <- n.segments - 1
  }

  # Apply PCA if number of features > 10
  if (ncol(segments_data) > 10) {
    message("Number of features > 10, performing PCA...")
    pca_result <- prcomp(segments_data, scale. = TRUE)
    segments_data <- pca_result$x[, 1:min(n.pcs, ncol(pca_result$x))]
  }

  # Find k-nearest neighbors using RANN
  message("Computing nearest neighbors...")
  nn.results <- RANN::nn2(
    data = segments_data,
    k = k.param
  )

  # Convert to graph
  message("Converting to neighbor graph...")
  nn.idx <- nn.results$nn.idx
  j <- as.numeric(t(nn.idx))
  i <- ((1:length(j)) - 1) %/% k.param + 1

  # Create nearest neighbors matrix
  nn.matrix <-  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(n.segments, n.segments)
  )
  rownames(nn.matrix) <- rownames(segments_data)
  colnames(nn.matrix) <- rownames(segments_data)

  neighbor.graphs <- list(nn = nn.matrix)

  # Compute SNN using C++ function
  if (compute.SNN) {
    message("Computing SNN graph...")
    snn.matrix <- ComputeSNN(nn.idx, prune.SNN)
    rownames(snn.matrix) <- rownames(segments_data)
    colnames(snn.matrix) <- rownames(segments_data)
    neighbor.graphs$snn <- snn.matrix
  }

  return(neighbor.graphs)
}


#' Find Clusters Using Modularity Optimization
#'
#' @description
#' Internal function to perform cluster detection.
#'
#' @param neighbor_graphs List of neighbor graphs
#' @param modularity.fxn Modularity function type
#' @param resolution Resolution parameter
#' @param algorithm Clustering algorithm
#' @param n.start Number of random starts
#' @param n.iter Maximum iterations
#' @param random.seed Random seed
#' @param group.singletons Whether to group singleton clusters
#' @param verbose Print progress messages
#'
#' @return Data frame with cluster assignments
#'
#' @keywords internal
find_segment_clusters <- function(
    neighbor_graphs,
    modularity.fxn = 1,
    resolution = 0.8,
    algorithm = 1,  # 1 = original Louvain
    n.start = 10,
    n.iter = 10,
    random.seed = 0,
    group.singletons = TRUE,
    verbose = TRUE
) {
  if (!"snn" %in% names(neighbor_graphs)) {
    stop("SNN graph required for clustering")
  }

  # Run modularity clustering using Seurat's implementation
  cluster_ids <- RunModularityClustering(
    SNN = neighbor_graphs$snn,
    modularity = modularity.fxn,
    resolution = resolution,
    algorithm = algorithm,
    n.start = n.start,
    n.iter = n.iter,
    random.seed = random.seed,
    print.output = verbose
  )

  # Name the clusters
  names(cluster_ids) <- rownames(neighbor_graphs$snn)

  # Handle singletons using Seurat's function
  if (group.singletons) {
    cluster_ids <- GroupSingletons(
      ids = cluster_ids,
      SNN = neighbor_graphs$snn,
      group.singletons = TRUE,
      verbose = verbose
    )
  }

  # Create results data frame
  results <- data.frame(
    row.names = names(cluster_ids)
  )
  results[paste0("res.", resolution)] <- factor(cluster_ids)

  return(results)
}

#' Run Modularity Clustering
#'
#' @description
#' Internal function to run modularity optimization.
#'
#' @keywords internal
RunModularityClustering <- function(
    SNN = matrix(),
    modularity = 1,
    resolution = 0.8,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = 0,
    print.output = TRUE,
    temp.file.location = NULL,
    edge.file.name = NULL
) {
  edge_file <- edge.file.name %||% ''
  clusters <- RunModularityClusteringCpp(
    SNN,
    modularity,
    resolution,
    algorithm,
    n.start,
    n.iter,
    random.seed,
    print.output,
    edge_file
  )
  return(clusters)
}

#' Group Singleton Clusters
#'
#' @description
#' Internal function to handle singleton clusters.
#'
#' @keywords internal
GroupSingletons <- function(ids, SNN, group.singletons = TRUE, verbose = TRUE) {
  # identify singletons
  singletons <- c()
  singletons <- names(x = which(x = table(ids) == 1))
  singletons <- intersect(x = unique(x = ids), singletons)
  if (!group.singletons) {
    ids[which(ids %in% singletons)] <- "singleton"
    return(ids)
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(x = unique(x = ids))
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  new.ids <- ids
  for (i in singletons) {
    i.cells <- names(which(ids == i))
    for (j in cluster_names) {
      j.cells <- names(which(ids == j))
      subSNN <- SNN[i.cells, j.cells]
      set.seed(1) # to match previous behavior, random seed being set in WhichCells
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }
  return(ids)
}
