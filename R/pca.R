# Principal Component Analysis
# Update date : Feb. 4, 2025

#' Run Principal Component Analysis
#'
#' @description
#' A generic function to perform PCA with support for multiple methods
#' and large-scale data processing.
#'
#' @param x An object to analyze, either a matrix/data frame
#'          or a SAP object
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' This generic function supports PCA through two methods:
#' \itemize{
#'   \item Default method for matrices with multiple PCA implementations
#'   \item SAP object method for organized trajectory data
#' }
#'
#' @return
#' PCA results or updated SAP object
#'
#' @examples
#' \dontrun{
#' # Run PCA on matrix
#' pca <- run_pca(matrix, method = "irlba", n_components = 50)
#'
#' # Run PCA on SAP object
#' sap_obj <- run_pca(sap_object,
#'                    segment_type = "motifs",
#'                    data_type = "traj_mat")
#' }
#'
#' @export
run_pca <- function(x, ...) {
  UseMethod("run_pca")
}


#' Run PCA on Matrix Data
#'
#' @description
#' Performs PCA using various methods optimized for different data scales.
#'
#' @param x Matrix or data frame to analyze
#' @param method PCA method ('irlba', 'base', or 'parallel')
#' @param n_components Number of principal components
#' @param n_cores Number of cores for parallel processing
#' @param diagnostic_plots Whether to create diagnostic plots
#' @param output_scores Whether to return PC scores
#' @param scale_scores Whether to scale PC scores
#' @param ... Additional arguments
#'
#' @details
#' Supports three PCA methods:
#' \itemize{
#'   \item irlba: Fast truncated SVD for large matrices
#'   \item base: Standard R prcomp implementation
#'   \item parallel: Distributed computation for very large matrices
#' }
#'
#' @return
#' PCA results or PC scores matrix
#'
#' @export
run_pca.default <- function(x,
                   method = c("irlba", "base", "parallel"),
                   n_components = 50,
                   n_cores = NULL,
                   diagnostic_plots = TRUE,
                   output_scores = TRUE,
                   scale_scores = FALSE) {


  method <- match.arg(method)

  # Set default cores if NULL
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  # Validate inputs
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  # Check if n_components is greater than or equal to number of columns
  if (n_components >= ncol(x)) {
    stop(sprintf("n_components (%d) must be less than number of columns (%d) in input matrix",
                 n_components, ncol(x)))
  }

  if (n_components > min(dim(x))) {
    warning("n_components larger than minimum dimension, reducing to ", min(dim(x)))
    n_components <- min(dim(x))
  }

  # # Start timing
  # tic(paste("Time elapsed for", method, "PCA:"))

  result <- tryCatch({
    switch(
      method,
      "base" = {
        message("Performing PCA using r base prcomp method...")
        prcomp(x, scale. = TRUE)
      },
      "irlba" = {
        message("Performing PCA using irlba method...")
        fast_pca(x, n_components, scale = TRUE)
      },
      "parallel" = {
        message("Performing PCA using parallel processing method...")
        parallel_pca(x, n_components, n_cores = n_cores, scale = TRUE)
      }
    )
  }, error = function(e) {
    message("Error in PCA computation: ", e$message)
    return(NULL)
  })

  # # Print timing
  # toc()

  # Add diagnostic plots and analysis if requested
  if (diagnostic_plots && !is.null(result)) {
    # Calculate variance explained based on method
    var_explained <- switch(
      method,
      "base" = {
        (result$sdev^2) / sum(result$sdev^2)
      },
      "irlba" = {
        (result$sdev^2) / sum(result$sdev^2)
      },
      "parallel" = {
        (result$d^2) / sum(result$d^2)
      }
    )

    # Create diagnostic plots
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par))  # Restore original par settings

    # Scree plot
    plot(var_explained[1:n_components], type = "b",
         xlab = "Principal Component",
         ylab = "Proportion of Variance Explained",
         main = "Scree Plot")

    # Cumulative variance plot
    cumulative_var <- cumsum(var_explained)
    plot(cumulative_var[1:n_components], type = "b",
         xlab = "Number of Components",
         ylab = "Cumulative Proportion of Variance Explained",
         main = "Cumulative Variance")

    # Print summary statistics
    cat("\nPCA Diagnostic Summary:\n")
    cat("------------------------\n")
    cat("Variance explained by first PC:", round(var_explained[1]*100, 2), "%\n")
    cat("Variance explained by last PC:", round(var_explained[n_components]*100, 2), "%\n")
    cat("Ratio between first and last PC:",
        round(var_explained[1]/var_explained[n_components], 2), "\n")
    cat("Number of PCs needed for 80% variance:",
        which(cumulative_var >= 0.8)[1], "\n")

    # Decision making for subsequent analysis
    if(var_explained[1] / var_explained[n_components] > 10) {
      cat("\nRECOMMENDATION: Consider scaling before further dimension reduction because:\n")
      cat("- First PC explains much more variance than last PC\n")
      cat("- Scaling will help consider patterns in lower PCs\n")
    } else {
      cat("\nRECOMMENDATION: Scaling might not be necessary because:\n")
      cat("- Variance is relatively evenly distributed across PCs\n")
    }
  }

  # Handle output based on output_scores parameter
  if (output_scores && !is.null(result)) {
    # Extract scores based on method
    scores <- switch(
      method,
      "base" = {
        result$x[, 1:n_components, drop = FALSE]
      },
      "irlba" = {
        result$x[, 1:n_components, drop = FALSE] %*% diag(result$sdev[1:n_components])
      },
      "parallel" = {
        result$u[, 1:n_components, drop = FALSE] %*% diag(result$d[1:n_components])
      }
    )

    # Scale scores if requested
    if (scale_scores) {
      scores <- scale(scores)
      cat("\nPCA scores are scaled.")
    }
    return(scores)
  } else {
    return(result)
  }
}

#' Run PCA on SAP Object
#'
#' @description
#' Performs PCA on trajectory matrices stored in a SAP object.
#'
#' @param x A SAP object
#' @param segment_type Type of segments to analyze
#' @param data_type Type of data to analyze
#' @param method PCA method to use
#' @param n_components Number of components
#' @param n_cores Number of cores
#' @param diagnostic_plots Whether to show diagnostics
#' @param scale_scores Whether to scale scores
#' @param verbose Whether to print progress
#' @param ... Additional arguments
#'
#' @details
#' Performs PCA with the following features:
#' \itemize{
#'   \item Support for different segment types
#'   \item Multiple PCA methods
#'   \item Diagnostic visualizations
#'   \item Results storage in SAP object
#' }
#'
#' @return
#' Updated SAP object with PCA results stored in features slot
#'
#' @export
run_pca.Sap <- function(x,
                        segment_type = c("motifs", "syllables", "bouts", "segments"),
                        data_type = "traj_mat",
                        method = c("irlba", "base", "parallel"),
                        n_components = 50,
                        n_cores = NULL,
                        diagnostic_plots = TRUE,
                        scale_scores = FALSE,
                        verbose = TRUE,
                        ...) {
  if(verbose) message(sprintf("\n=== Starting PCA analysis for %s using %s (method: %s) ===\n",
                              segment_type[1],
                              data_type[1],
                              method[1]))

  # Validate segment_type
  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)  # Remove 's' from end

  # Get trajectory matrix
  traj_mat <- x$features[[feature_type]][[data_type]]
  if(is.null(traj_mat)) {
    stop(sprintf("No %s data found in %s features", data_type, feature_type))
  }

  # Run PCA using default method
  pca_scores <- run_pca.default(
    x = traj_mat,
    method = method,
    n_components = n_components,
    n_cores = n_cores,
    diagnostic_plots = diagnostic_plots,
    output_scores = TRUE,
    scale_scores = scale_scores,
    ...
  )

  # Save PCA scores to traj.embeds
  if(!is.null(pca_scores)) {
    traj_mat <- as.data.frame(pca_scores)
    colnames(traj_mat) <- paste0("PC", 1:ncol(pca_scores))

    # Update traj_mat in Sap object
    x$features[[feature_type]][["traj_mat"]] <- traj_mat

    # Add PCA parameters as attributes
    attr(x$features[[feature_type]][["traj_mat"]], "pca_params") <- list(
      method = method[1],
      n_components = n_components,
      scale_scores = scale_scores
    )
  }

  if (verbose) {
    message("\nPCA score stored in traj_mat")
    message("Access via: sap$features$", feature_type, "$traj_mat")
    message(sprintf("\n- Access PCA parameters via: attributes(sap$features$%s$traj_mat)$pca_params", feature_type))
  }

  # Return modified Sap object
  invisible(x)
}


#' Fast PCA Using IRLBA
#'
#' @description
#' Internal function for fast PCA using truncated SVD.
#'
#' @param matrix Input matrix
#' @param n_components Number of components
#' @param center Whether to center data
#' @param scale Whether to scale data
#'
#' @importFrom irlba prcomp_irlba
#' @keywords internal
fast_pca <- function(matrix, n_components = 50, center = TRUE, scale = TRUE) {
  #require(irlba)
  pca_result <- prcomp_irlba(matrix,
                             n = n_components,
                             center = center,
                             scale. = scale)

  return(pca_result)
}

#' Parallel PCA Implementation
#'
#' @description
#' Internal function for parallel PCA computation.
#'
#' @param x Input matrix
#' @param n_components Number of components
#' @param n_cores Number of cores
#' @param scale Whether to scale data
#'
#' @keywords internal
parallel_pca <- function(x, n_components = 50, n_cores = NULL, scale = TRUE) {
  # Check if bigstatsr is installed
  if (!requireNamespace("bigstatsr", quietly = TRUE)) {
    warning("Package 'bigstatsr' is required for parallel PCA processing.
            Please install it with: install.packages('bigstatsr')")
    return(NULL)
  }

  # Check if input is a valid matrix
  if (!is.matrix(x)) {
    stop("Input 'x' must be a matrix")
  }

  # Safer core detection
  if (is.null(n_cores)) {
    n_cores <- min(bigstatsr::nb_cores(), 8)  # Setting a reasonable default maximum
  } else {
    n_cores <- min(n_cores, bigstatsr::nb_cores())
  }

  # Convert matrix to a Filebacked Big Matrix (FBM)
  message("Converting to big matrix format...")
  tryCatch({
    FBM <- bigstatsr::as_FBM(x)
  }, error = function(e) {
    stop("Failed to convert to FBM format: ", e$message)
  })

  # Set up scaling function
  scaling_fun <- if(scale) {
    message("Setting up scaling...")
    bigstatsr::big_scale(center = TRUE, scale = TRUE)
  } else {
    bigstatsr::big_scale(center = FALSE, scale = FALSE)
  }

  # Perform PCA
  message(sprintf("Performing PCA using %d cores...", n_cores))
  tryCatch({
    pca_result <- bigstatsr::big_randomSVD(
      X = FBM,
      fun.scaling = scaling_fun,
      k = min(n_components, ncol(FBM)-1),
      ncores = n_cores
    )
    return(pca_result)
  }, error = function(e) {
    message("Error with parallel processing, falling back to single core...")
    pca_result <- bigstatsr::big_randomSVD(
      X = FBM,
      fun.scaling = scaling_fun,
      k = min(n_components, ncol(FBM)-1),
      ncores = 1
    )
    return(pca_result)
  })
}
