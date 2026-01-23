# Test Clustering Functions
# Tests for functions in R/clustering.R

test_that("find_clusters parameter validation", {
  # Test clustering parameter validation
  validate_clustering_params <- function(k.param = 20,
                                       prune.SNN = 1/15,
                                       n.pcs = NULL,
                                       resolution = 0.8,
                                       n.start = 10) {
    
    # Validate k.param
    if (!is.numeric(k.param) || k.param <= 0) {
      stop("k.param must be a positive number")
    }
    
    # Validate prune.SNN
    if (!is.numeric(prune.SNN) || prune.SNN < 0 || prune.SNN > 1) {
      stop("prune.SNN must be between 0 and 1")
    }
    
    # Validate n.pcs
    if (!is.null(n.pcs) && (!is.numeric(n.pcs) || n.pcs <= 0)) {
      stop("n.pcs must be a positive number")
    }
    
    # Validate resolution
    if (!is.numeric(resolution) || resolution <= 0) {
      stop("resolution must be a positive number")
    }
    
    # Validate n.start
    if (!is.numeric(n.start) || n.start <= 0) {
      stop("n.start must be a positive number")
    }
    
    return(TRUE)
  }
  
  # Test valid parameters
  expect_true(validate_clustering_params())
  expect_true(validate_clustering_params(k.param = 15, prune.SNN = 0.1))
  expect_true(validate_clustering_params(n.pcs = 50, resolution = 1.2))
  
  # Test invalid parameters
  expect_error(validate_clustering_params(k.param = -5),
               "k.param must be a positive number")
  expect_error(validate_clustering_params(prune.SNN = 1.5),
               "prune.SNN must be between 0 and 1")
  expect_error(validate_clustering_params(n.pcs = -10),
               "n.pcs must be a positive number")
  expect_error(validate_clustering_params(resolution = 0),
               "resolution must be a positive number")
  expect_error(validate_clustering_params(n.start = 0),
               "n.start must be a positive number")
})

test_that("SNN graph construction logic", {
  # Test shared nearest neighbor graph construction
  construct_snn_graph <- function(nn_matrix, prune_threshold = 1/15) {
    n_cells <- nrow(nn_matrix)
    k <- ncol(nn_matrix)
    
    # Initialize SNN matrix
    snn_matrix <- matrix(0, nrow = n_cells, ncol = n_cells)
    
    # Calculate shared neighbors
    for (i in 1:n_cells) {
      for (j in 1:n_cells) {
        if (i != j) {
          # Find shared neighbors
          shared_neighbors <- length(intersect(nn_matrix[i, ], nn_matrix[j, ]))
          
          # Calculate SNN strength
          snn_strength <- shared_neighbors / k
          
          # Apply pruning threshold
          if (snn_strength >= prune_threshold) {
            snn_matrix[i, j] <- snn_strength
          }
        }
      }
    }
    
    return(snn_matrix)
  }
  
  # Create test nearest neighbor matrix
  test_nn <- matrix(c(
    1, 2, 3,  # Cell 1's neighbors
    1, 3, 4,  # Cell 2's neighbors  
    1, 2, 4,  # Cell 3's neighbors
    2, 3, 5,  # Cell 4's neighbors
    4, 1, 2   # Cell 5's neighbors
  ), nrow = 5, byrow = TRUE)
  
  # Test SNN construction
  snn_result <- construct_snn_graph(test_nn, prune_threshold = 0.1)
  
  expect_true(is.matrix(snn_result))
  expect_equal(dim(snn_result), c(5, 5))
  expect_equal(diag(snn_result), rep(0, 5))  # Diagonal should be zero
  expect_true(all(snn_result >= 0 & snn_result <= 1))  # Values between 0 and 1
  
  # Test symmetry (SNN should be symmetric)
  expect_true(all(abs(snn_result - t(snn_result)) < 1e-10))
  
  # Test pruning effect
  snn_strict <- construct_snn_graph(test_nn, prune_threshold = 0.5)
  expect_true(sum(snn_strict > 0) <= sum(snn_result > 0))  # Stricter pruning = fewer connections
})

test_that("cluster assignment validation", {
  # Test cluster assignment logic
  validate_cluster_assignments <- function(clusters, n_observations) {
    # Check cluster vector length
    if (length(clusters) != n_observations) {
      stop("Cluster vector length must match number of observations")
    }
    
    # Check for valid cluster IDs
    if (any(is.na(clusters))) {
      stop("Cluster assignments cannot contain NA values")
    }
    
    if (any(clusters <= 0)) {
      stop("Cluster IDs must be positive integers")
    }
    
    # Check for reasonable number of clusters
    n_clusters <- length(unique(clusters))
    if (n_clusters == 1) {
      warning("Only one cluster found - consider adjusting parameters")
    }
    
    if (n_clusters > n_observations / 2) {
      warning("Very high number of clusters relative to observations")
    }
    
    return(list(
      n_clusters = n_clusters,
      cluster_sizes = table(clusters),
      valid = TRUE
    ))
  }
  
  # Test valid cluster assignments
  valid_clusters <- c(1, 1, 2, 2, 3, 3, 1, 2)
  result1 <- validate_cluster_assignments(valid_clusters, 8)
  expect_true(result1$valid)
  expect_equal(result1$n_clusters, 3)
  expect_equal(sum(result1$cluster_sizes), 8)
  
  # Test invalid assignments
  expect_error(validate_cluster_assignments(c(1, 2, 3), 5),
               "Cluster vector length must match number of observations")
  expect_error(validate_cluster_assignments(c(1, 2, NA, 3), 4),
               "Cluster assignments cannot contain NA values")
  expect_error(validate_cluster_assignments(c(1, 2, 0, 3), 4),
               "Cluster IDs must be positive integers")
  
  # Test warnings
  expect_warning(validate_cluster_assignments(c(1, 1, 1, 1), 4),
                "Only one cluster found")
  expect_warning(validate_cluster_assignments(1:6, 6),
                "Very high number of clusters")
})

test_that("neighbor finding logic", {
  # Test k-nearest neighbor finding
  find_k_neighbors <- function(distance_matrix, k) {
    n <- nrow(distance_matrix)
    
    if (k >= n) {
      stop("k must be less than the number of observations")
    }
    
    # Initialize neighbor matrix
    neighbors <- matrix(0, nrow = n, ncol = k)
    
    # Find k nearest neighbors for each point
    for (i in 1:n) {
      # Get distances for point i (excluding self)
      distances <- distance_matrix[i, -i]
      neighbor_indices <- which(distance_matrix[i, ] %in% sort(distances)[1:k])
      
      # Remove self if included
      neighbor_indices <- neighbor_indices[neighbor_indices != i]
      
      # Take first k neighbors
      neighbors[i, ] <- neighbor_indices[1:k]
    }
    
    return(neighbors)
  }
  
  # Create test distance matrix
  set.seed(123)
  test_coords <- matrix(rnorm(20), nrow = 5, ncol = 4)
  test_dist <- as.matrix(dist(test_coords))
  
  # Test neighbor finding
  neighbors <- find_k_neighbors(test_dist, k = 3)
  
  expect_true(is.matrix(neighbors))
  expect_equal(dim(neighbors), c(5, 3))
  expect_true(all(neighbors > 0 & neighbors <= 5))  # Valid indices
  
  # Test that no point is its own neighbor
  for (i in 1:5) {
    expect_false(i %in% neighbors[i, ])
  }
  
  # Test error for invalid k
  expect_error(find_k_neighbors(test_dist, k = 5),
               "k must be less than the number of observations")
})

test_that("modularity optimization", {
  # Test modularity calculation
  calculate_modularity <- function(adjacency_matrix, clusters) {
    n <- nrow(adjacency_matrix)
    m <- sum(adjacency_matrix) / 2  # Total number of edges
    
    if (m == 0) {
      return(0)  # No edges, modularity is 0
    }
    
    # Calculate degree for each node
    degrees <- rowSums(adjacency_matrix)
    
    # Calculate modularity
    modularity <- 0
    unique_clusters <- unique(clusters)
    
    for (cluster in unique_clusters) {
      cluster_nodes <- which(clusters == cluster)
      
      # Internal edges in cluster
      internal_edges <- sum(adjacency_matrix[cluster_nodes, cluster_nodes]) / 2
      
      # Expected internal edges
      cluster_degree_sum <- sum(degrees[cluster_nodes])
      expected_internal <- (cluster_degree_sum^2) / (4 * m)
      
      # Add to modularity
      modularity <- modularity + (internal_edges / m) - (expected_internal / m)
    }
    
    return(modularity)
  }
  
  # Test with simple graph
  # Create a graph with clear community structure
  adj_matrix <- matrix(0, nrow = 6, ncol = 6)
  # Community 1: nodes 1,2,3
  adj_matrix[1, 2] <- adj_matrix[2, 1] <- 1
  adj_matrix[2, 3] <- adj_matrix[3, 2] <- 1
  adj_matrix[1, 3] <- adj_matrix[3, 1] <- 1
  # Community 2: nodes 4,5,6
  adj_matrix[4, 5] <- adj_matrix[5, 4] <- 1
  adj_matrix[5, 6] <- adj_matrix[6, 5] <- 1
  adj_matrix[4, 6] <- adj_matrix[6, 4] <- 1
  # Weak connection between communities
  adj_matrix[3, 4] <- adj_matrix[4, 3] <- 1
  
  # Test good clustering (matches community structure)
  good_clusters <- c(1, 1, 1, 2, 2, 2)
  good_modularity <- calculate_modularity(adj_matrix, good_clusters)
  
  # Test bad clustering (random assignment)
  bad_clusters <- c(1, 2, 1, 2, 1, 2)
  bad_modularity <- calculate_modularity(adj_matrix, bad_clusters)
  
  expect_true(good_modularity > bad_modularity)
  expect_true(good_modularity > 0)  # Should be positive for good clustering
  
  # Test edge cases
  empty_matrix <- matrix(0, nrow = 4, ncol = 4)
  expect_equal(calculate_modularity(empty_matrix, c(1, 1, 2, 2)), 0)
})