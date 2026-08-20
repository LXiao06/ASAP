#' Construct file path for audio file
#'
#' Uses the \code{subfolder} column (the original directory name) if present,
#' otherwise falls back to \code{day_post_hatch}. This is important when
#' \code{day_post_hatch} was supplied manually by the user and does not match
#' the actual folder on disk.
#'
#' @param x A data frame row with file name
#' @param wav_dir Path to WAV files directory (default: NULL)
#'
#' @return Character string of full file path
#' @noRd
#' @keywords internal
construct_wav_path <- function(x,
                               wav_dir = NULL) {
  # Check required column
  if (!("filename" %in% names(x))) {
    stop("Input must contain 'filename' column")
  }

  # Handle wav_dir
  if (is.null(wav_dir) && is.null(attr(x, "wav_dir"))) {
    stop("wav_dir must be provided either as argument or attribute")
  }

  # Determine the subdirectory component: prefer 'subfolder' over 'day_post_hatch'
  subdir <- if ("subfolder" %in% names(x) && !is.na(x$subfolder)) {
    x$subfolder
  } else if ("day_post_hatch" %in% names(x) && !is.na(x$day_post_hatch)) {
    x$day_post_hatch
  } else {
    NULL
  }

  # Construct file path based on available information
  base <- if (!is.null(attr(x, "wav_dir"))) attr(x, "wav_dir") else wav_dir

  if (!is.null(subdir)) {
    sound_path <- file.path(base, subdir, x$filename)
  } else {
    sound_path <- file.path(base, x$filename)
  }

  # Verify file exists
  if (!file.exists(sound_path)) {
    stop(sprintf("Sound file not found: %s", sound_path))
  }

  return(sound_path)
}


#' Select and balance segments from a data frame
#'
#' @param segments_df Data frame containing segments
#' @param labels Optional vector of labels to select
#' @param balanced Logical. Whether to balance groups
#' @param sample_percent Numeric. Percentage of segments to sample from each group
#' @param max_samples_per_label Optional integer. Maximum number of segments to
#'   randomly sample from each label when \code{balanced = FALSE}.
#' @param seed Integer. Seed for random sampling (default: 222)
#'
#' @return List containing modified segments data frame and summary information
#' @noRd
#'
#' @keywords internal
select_segments <- function(segments_df,
                            labels = NULL,
                            clusters = NULL,
                            balanced = FALSE,
                            sample_percent = NULL,
                            max_samples_per_label = NULL,
                            seed = 222) {
  # Set seed
  set.seed(seed)

  # Handle cluster filtering first if specified
  if (!is.null(clusters)) {
    if (!"cluster" %in% names(segments_df)) {
      stop("No 'cluster' column found in the data")
    }

    available_clusters <- unique(segments_df$cluster)
    missing_clusters <- clusters[!clusters %in% available_clusters]

    if (length(missing_clusters) > 0) {
      stop(sprintf(
        "Clusters not found: %s\nAvailable clusters: %s",
        paste(missing_clusters, collapse = ", "),
        paste(available_clusters, collapse = ", ")
      ))
    }

    segments_df <- segments_df[segments_df$cluster %in% clusters, ]
    message(sprintf("\nFiltered for clusters: %s", paste(clusters, collapse = ", ")))
  }

  # Handle label selection
  if (!is.null(labels)) {
    available_labels <- unique(segments_df$label)
    missing_labels <- labels[!labels %in% available_labels]
    if (length(missing_labels) > 0) {
      stop(sprintf(
        "Labels not found: %s\nAvailable labels: %s",
        paste(missing_labels, collapse = ", "),
        paste(available_labels, collapse = ", ")
      ))
    }
    segments_df <- segments_df[segments_df$label %in% labels, ]
  }

  # Add order tracking index after initial filtering
  segments_df <- segments_df |>
    dplyr::mutate(.original_order = dplyr::row_number())

  message("Available labels: ", paste(unique(segments_df$label), collapse = ", "))

  # Print initial summary
  message("\nInitial data summary:")
  initial_summary <- segments_df |>
    dplyr::group_by(label) |>
    dplyr::summarise(
      dph = mean(day_post_hatch),
      n = dplyr::n(),
      time_window = mean(end_time - start_time)
    ) |>
    dplyr::arrange(dph)
  print(initial_summary)

  # Handle balancing first if requested
  if (balanced) {
    if (!is.null(max_samples_per_label)) {
      message("\nIgnoring max_samples_per_label because balanced = TRUE")
    }

    group_counts <- segments_df |>
      dplyr::count(label)
    n_sample <- floor(min(group_counts$n) / 10) * 10 # Round to nearest 10

    message(sprintf("\nBalancing groups to %d segments each", n_sample))

    segments_df <- segments_df |>
      dplyr::group_by(label) |>
      dplyr::slice_sample(n = n_sample) |>
      dplyr::arrange(.data$.original_order) |>
      dplyr::ungroup() |>
      dplyr::arrange(day_post_hatch, .data$.original_order)
  }

  # Handle sampling if specified
  if (!is.null(sample_percent)) {
    if (!is.numeric(sample_percent) || sample_percent <= 0 || sample_percent > 100) {
      stop("sample_percent must be between 0 and 100")
    }

    segments_df <- segments_df |>
      dplyr::group_by(label) |>
      dplyr::group_modify(~ {
        n_available <- nrow(.x)
        n_to_sample <- ceiling(n_available * sample_percent / 100)
        dplyr::slice_sample(.x, n = n_to_sample) |>
          dplyr::arrange(.data$.original_order)
      }) |>
      dplyr::ungroup() |>
      dplyr::arrange(day_post_hatch)

    message(sprintf("\nSampling %.1f%% from each label", sample_percent))
  }

  # Cap samples per label only for unbalanced workflows.
  if (!balanced && !is.null(max_samples_per_label)) {
    if (!is.numeric(max_samples_per_label) ||
        length(max_samples_per_label) != 1 ||
        is.na(max_samples_per_label) ||
        max_samples_per_label <= 0 ||
        max_samples_per_label != floor(max_samples_per_label)) {
      stop("max_samples_per_label must be a positive integer")
    }

    segments_df <- segments_df |>
      dplyr::group_by(label) |>
      dplyr::group_modify(~ {
        n_available <- nrow(.x)
        n_to_sample <- min(as.integer(max_samples_per_label), n_available)
        dplyr::slice_sample(.x, n = n_to_sample) |>
          dplyr::arrange(.data$.original_order)
      }) |>
      dplyr::ungroup() |>
      dplyr::arrange(day_post_hatch)

    message(sprintf(
      "\nSampling up to %d segments from each label",
      as.integer(max_samples_per_label)
    ))
  }

  # Final cleanup and ordering
  segments_df <- segments_df |>
    dplyr::arrange(.data$.original_order) |> # Global order preservation
    dplyr::select(-dplyr::all_of(".original_order"))

  # Print final summary if data was modified
  if (balanced || !is.null(sample_percent) || !is.null(max_samples_per_label)) {
    message("\nFinal data summary:")
    final_summary <- segments_df |>
      dplyr::group_by(label) |>
      dplyr::summarise(
        dph = mean(day_post_hatch),
        n = dplyr::n(),
        time_window = mean(end_time - start_time)
      ) |>
      dplyr::arrange(dph)
    print(final_summary)
  }

  # #  Store day ordering as an attribute
  # attr(segments_df, "label_order") <- unique(segments_df$label)

  return(segments_df)
}



#' parallel processing function
#'
#' @keywords internal
parallel_apply <- function(indices, FUN, cores, use_preschedule = FALSE, cl = NULL) {
  # Set number of cores
  if (is.null(cores)) {
    ensure_pkgs("parallel")
    detected <- suppressWarnings(as.integer(parallel::detectCores()))
    if (is.na(detected) || detected < 1) {
      cores <- 1L
    } else {
      cores <- max(1L, detected - 1L)
    }
  }

  # Determine which package is needed based on OS and environment
  if (cores > 1) {
    if (Sys.info()["sysname"] == "Darwin") {
      # macOS: use fork-based pbmclapply (fast and reliable on macOS)
      ensure_pkgs("pbmcapply")
      result <- pbmcapply::pbmclapply(
        indices,
        FUN,
        mc.cores = cores,
        mc.preschedule = use_preschedule
      )
    } else {
      # Linux/Windows: use PSOCK clusters.
      # If a pre-created cluster is passed in, reuse it (avoids repeated
      # makeCluster startup cost when called in a loop over many days).
      ensure_pkgs("pbapply", "parallel")
      if (is.null(cl)) {
        cl <- parallel::makeCluster(cores, type = "PSOCK")
        parallel::clusterEvalQ(cl, loadNamespace("ASAP"))
        on.exit(parallel::stopCluster(cl), add = TRUE)
      }
      result <- pbapply::pblapply(
        indices,
        FUN,
        cl = cl
      )
    }
  } else {
    # Single-core with progress still uses pbapply
    ensure_pkgs("pbapply")
    result <- pbapply::pblapply(
      indices,
      FUN
    )
  }

  return(result)
}


#' Check Python dependencies
#'
#' @keywords internal
check_python_dependencies <- function(verbose = FALSE) {
  # Check for reticulate package
  ensure_pkgs("reticulate")

  # Check for librosa
  tryCatch(
    {
      reticulate::py_module_available("librosa")
    },
    error = function(e) {
      if (verbose) {
        message("Librosa is not installed. Attempting to install...")
      }

      # Try to install librosa
      tryCatch(
        {
          reticulate::py_install("librosa", pip = TRUE)
        },
        error = function(e) {
          stop("Failed to install librosa. Please install it manually using pip: pip install librosa")
        }
      )
    }
  )

  # Check for numpy
  tryCatch(
    {
      reticulate::py_module_available("numpy")
    },
    error = function(e) {
      if (verbose) {
        message("Numpy is not installed. Attempting to install...")
      }

      # Try to install numpy
      tryCatch(
        {
          reticulate::py_install("numpy", pip = TRUE)
        },
        error = function(e) {
          stop("Failed to install numpy. Please install it manually using pip: pip install numpy")
        }
      )
    }
  )
}

#' @keywords internal
calculate_segment_stats <- function(feature_matrix,
                                    template,
                                    time_step,
                                    cores = NULL,
                                    ...) {
  # Validate inputs
  if (!is.matrix(feature_matrix)) stop("feature_matrix must be a matrix")
  if (!is.logical(template)) stop("template must be a logical vector")
  if (nrow(feature_matrix) != length(template)) {
    stop("feature_matrix rows must match template length")
  }

  # Get unique labels
  labels <- unique(colnames(feature_matrix))

  # Find valid segments
  runs <- rle(template)
  segment_starts <- cumsum(c(1, runs$lengths[-length(runs$lengths)]))
  segment_ends <- cumsum(runs$lengths)
  valid_segments <- which(runs$values)

  # Process all labels in parallel
  all_stats <- parallel_apply(
    labels,
    function(label) {
      process_label(
        label, feature_matrix, segment_starts, segment_ends,
        valid_segments, time_step
      )
    },
    cores = cores
  )

  # Combine results
  stats_df <- do.call(rbind, unlist(all_stats, recursive = FALSE))
  rownames(stats_df) <- NULL

  return(stats_df)
}

#' @keywords internal
process_label <- function(label, feature_matrix, segment_starts, segment_ends,
                          valid_segments, time_step) {
  # Get columns belonging to this label
  label_cols <- which(colnames(feature_matrix) == label)
  if (length(label_cols) == 0) {
    return(NULL)
  }

  # Create all combinations of label columns and segments
  combinations <- expand.grid(
    col = label_cols,
    segment = valid_segments,
    stringsAsFactors = FALSE
  )

  # Function to process a single combination
  process_combination <- function(combo) {
    col <- combo$col
    i <- combo$segment

    seg_range <- segment_starts[i]:segment_ends[i]
    seg_data <- feature_matrix[seg_range, col]
    valid_vals <- seg_data[!is.na(seg_data)]

    if (length(valid_vals) > 0) {
      # Filter out zeros
      non_zero_vals <- valid_vals[valid_vals != 0]

      # Only proceed if we have more then two non-zero values for statistics
      if (length(non_zero_vals) > 2) {
        data.frame(
          label = label,
          segment_id = i,
          duration = (segment_ends[i] - segment_starts[i] + 1) * time_step,
          mean = mean(non_zero_vals),
          median = median(non_zero_vals),
          sd = sd(non_zero_vals),
          sem = sd(non_zero_vals) / sqrt(length(non_zero_vals)),
          min_val = min(non_zero_vals),
          max_val = max(non_zero_vals),
          n_samples = length(non_zero_vals),
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  }

  # Use parallel_apply to process combinations
  parallel_apply(
    1:nrow(combinations),
    function(idx) process_combination(combinations[idx, ]),
    cores = NULL # Will use default core detection
  )
}


#' Perform ANOVA and Multiple Comparisons Analysis
#'
#' @description
#' Performs one-way ANOVA and Tukey's HSD test for multiple comparisons across different segments.
#' Provides both statistical results and optional visualization.
#'
#' @param stats_df A data frame containing columns:
#'   \itemize{
#'     \item segment_id: Numeric identifier for segments
#'     \item label: Factor or character indicating groups to compare
#'     \item mean: Numeric values for comparison
#'   }
#' @param plot Logical, whether to create a boxplot visualization (default: TRUE)
#'
#' @return A tibble containing ANOVA results with columns:
#'   \itemize{
#'     \item segment_id: Segment identifier
#'     \item term: Source of variation (label or Residuals)
#'     \item df: Degrees of freedom
#'     \item sumsq: Sum of squares
#'     \item meansq: Mean squares
#'     \item statistic: F-statistic
#'     \item p.value: P-value
#'     \item significant: Logical indicating if p.value < 0.05
#'   }
#'
#' @details
#' The function performs two main analyses:
#' \itemize{
#'   \item One-way ANOVA for each segment
#'   \item Tukey's HSD test for multiple comparisons with adjusted p-values
#' }
#'
#' The printed output includes:
#' \itemize{
#'   \item Tukey's HSD results with adjusted p-values
#'   \item Significance levels: *** (p<0.001), ** (p<0.01), * (p<0.05), ns (p>=0.05)
#'   \item Optional boxplot visualization
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- anova_results(stats_df)
#'
#' # Without plot
#' results <- anova_results(stats_df, plot = FALSE)
#' }
#'
#' @export
anova_analysis <- function(stats_df, plot = TRUE) {
  # Store ANOVA and Tukey results
  results_list <- stats_df |>
    dplyr::group_by(.data$segment_id) |>
    dplyr::filter(dplyr::n_distinct(.data$label) > 1) |> # Only compare segments with multiple labels
    dplyr::group_split() |>
    lapply(function(seg_data) {
      # Perform ANOVA
      aov_model <- aov(mean ~ label, data = seg_data)

      # Create anova result data frame
      aov_summary <- summary(aov_model)[[1]]
      anova_result <- data.frame(
        term = rownames(aov_summary),
        Df = aov_summary$Df,
        Sum_Sq = aov_summary$`Sum Sq`,
        Mean_Sq = aov_summary$`Mean Sq`,
        F_value = aov_summary$`F value`,
        p.value = aov_summary$`Pr(>F)`,
        stringsAsFactors = FALSE
      ) |>
        dplyr::mutate(
          significant = .data$p.value < 0.05,
          segment_id = unique(seg_data$segment_id)
        ) |>
        dplyr::select(.data$segment_id, dplyr::everything())

      # Perform Tukey's HSD test for multiple comparison adjustment
      tukey_result <- TukeyHSD(aov_model)
      tukey_df <- as.data.frame(tukey_result$label)
      tukey_df <- tukey_df |>
        dplyr::mutate(
          comparison = rownames(tukey_df),
          segment_id = unique(seg_data$segment_id),
          significance = dplyr::case_when(
            `p adj` < 0.001 ~ "***",
            `p adj` < 0.01 ~ "**",
            `p adj` < 0.05 ~ "*",
            TRUE ~ "ns"
          )
        )


      list(anova = anova_result, tukey = tukey_df)
    })

  # Combine ANOVA results
  anova_summary <- dplyr::bind_rows(lapply(results_list, function(x) x$anova))

  # Print pretty output for Tukey's HSD results
  cat("## Tukey multiple comparisons of means")
  cat("\n## 95% family-wise confidence level")
  cat("\n## Fit: aov(formula = mean ~ label, data = segment)")
  cat("\n## Note: p-values are adjusted for multiple comparisons using Tukey's method\n")
  cat("\nlabel       segment      diff         lwr          upr         p adj   sign")

  # Print Tukey results for each segment
  for (result in results_list) {
    tukey_df <- result$tukey
    for (i in 1:nrow(tukey_df)) {
      cat(sprintf(
        "\n%-10s %5.0f %12.4f %12.4f %12.4f %12.5f    %s",
        tukey_df$comparison[i],
        tukey_df$segment_id[i],
        tukey_df$diff[i],
        tukey_df$lwr[i],
        tukey_df$upr[i],
        tukey_df$`p adj`[i],
        tukey_df$significance[i]
      ))
    }
  }

  # Add significance code explanation
  cat("\n\nSignificance codes (adjusted p-values):")
  cat("\n'***' < 0.001")
  cat("\n'**'  < 0.01")
  cat("\n'*'   < 0.05")
  cat("\n'ns'  >= 0.05\n")

  # Create plot if requested
  if (plot) {
    ensure_pkgs("ggplot2")
    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = .data$label, y = .data$mean)) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~ .data$segment_id) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Mean Values by Label across Segments",
        x = "Label",
        y = "Mean Value"
      )
    print(p)
  }

  # Return ANOVA summary
  return(anova_summary)
}

#' Check if a package is available
#'
#' @param pkg_name Character name of the package
#' @return TRUE if package is available, FALSE otherwise
#' @keywords internal
check_pkg <- function(pkg_name) {
  # Get cached status or check if not in cache
  has_var <- paste0("has_", pkg_name)
  if (is.null(.pkg_env[[has_var]])) {
    .pkg_env[[has_var]] <- requireNamespace(pkg_name, quietly = TRUE)
  }

  return(.pkg_env[[has_var]])
}


#' Ensure required packages are available, auto-installing if needed
#'
#' @param ... Character names of required packages
#' @return Nothing, called for side effects
#' @keywords internal
ensure_pkgs <- function(...) {
  pkgs <- c(...)
  missing_pkgs <- pkgs[!sapply(pkgs, check_pkg)]

  # Auto-install missing packages
  if (length(missing_pkgs) > 0) {
    message("Installing required packages: ", paste(missing_pkgs, collapse = ", "))

    install_results <- lapply(missing_pkgs, function(pkg) {
      tryCatch(
        {
          utils::install.packages(pkg)
          pkg_available <- requireNamespace(pkg, quietly = TRUE)
          # Update cached status
          .pkg_env[[paste0("has_", pkg)]] <- pkg_available
          return(pkg_available)
        },
        error = function(e) {
          message("Failed to install ", pkg, ": ", conditionMessage(e))
          return(FALSE)
        }
      )
    })

    still_missing <- missing_pkgs[!unlist(install_results)]

    # If any packages couldn't be installed, throw an error
    if (length(still_missing) > 0) {
      stop("Failed to install required packages: ",
        paste(still_missing, collapse = ", "),
        ". Please install them manually with install.packages().",
        call. = FALSE
      )
    }
  }

  # At this point, all packages should be available
  return(invisible(TRUE))
}

#' Get Indices of WAV Files in SAP Object Metadata
#'
#' @description
#' Finds the row indices of specified WAV files within the metadata of a Sound Analysis Pro (SAP) object.
#'
#' @param sap A SAP object
#' @param wav_files A character vector of WAV file names
#'
#' @return A numeric vector representing the row indices of the specified WAV files
#'   in the SAP object's metadata.
#'
#' @export
get_wav_indices <- function(sap, wav_files) {
  # Validate inputs
  if (!inherits(sap, "Sap")) {
    stop("Input 'sap' must be a Sap object.")
  }

  if (!is.character(wav_files) || length(wav_files) == 0) {
    stop("Input 'wav_files' must be a non-empty character vector.")
  }

  # Find indices where the filename matches the requested files
  indices <- which(sap$metadata$filename %in% wav_files)

  # Check for missing files
  if (length(indices) == 0) {
    warning("None of the specified WAV files were found in the SAP metadata.")
  } else {
    found_files <- unique(sap$metadata$filename[indices])
    missing_files <- setdiff(wav_files, found_files)
    if (length(missing_files) > 0) {
      warning(
        "The following WAV files were not found in the SAP metadata:\n",
        paste(missing_files, collapse = ", ")
      )
    }
  }

  return(indices)
}

#' List Numeric Subdirectory Names
#'
#' @description
#' Detects all subfolders within a specified directory that have purely numeric names
#' (e.g., "1", "100", "042") and returns them as a character vector.
#'
#' @param directory Path to the parent directory to search within
#'
#' @return A character vector of the numeric subfolder names. Returns an empty character
#'   vector if none are found.
#'
#' @export
list_numeric_dirs <- function(directory) {
  # Input validation
  if (!is.character(directory) || length(directory) != 1) {
    stop("Input 'directory' must be a single character string.")
  }
  
  if (!dir.exists(directory)) {
    stop(sprintf("Directory does not exist: %s", directory))
  }

  # List all subdirectories, extracting only the folder names (not full paths)
  subdirs <- list.dirs(directory, full.names = FALSE, recursive = FALSE)

  # Filter for strings that consist exclusively of numbers
  # ^        Start of string
  # [0-9]+   One or more digits
  # $        End of string
  numeric_subdirs <- subdirs[grepl("^[0-9]+$", subdirs)]

  # Sort numerically (e.g., so "20" comes before "100")
  if (length(numeric_subdirs) > 0) {
    numeric_subdirs <- numeric_subdirs[order(as.numeric(numeric_subdirs))]
  }

  return(as.character(numeric_subdirs))
}


# Discriminant Analysis ---------------------------------------------------

#' Canonical Discriminant Analysis with Robustness Testing
#'
#' @description
#' Performs Canonical Discriminant Analysis (CDA/LDA) on multivariate feature data
#' (such as trajectory similarities, acoustic feature statistics, or PCA embeddings)
#' to quantify and visualize group separation.
#'
#' @param data A data frame containing group labels and numeric feature columns
#' @param group_col Character. Name of the grouping variable (default: \code{"label"})
#' @param feature_cols Character vector of feature column names. If \code{NULL} (default),
#'   all numeric columns (excluding common metadata columns) are automatically used
#' @param reference_group Optional character. Reference group label (e.g., \code{"tutor"})
#'   to compute multivariate Mahalanobis distances from all groups to the reference
#' @param scale Logical. Whether to standardize features to zero-mean and unit-variance
#'   prior to analysis (default: \code{TRUE})
#' @param n_perm Integer. Number of label permutations for empirical robustness
#'   testing (default: \code{1000}). Set to \code{0} to skip permutation testing
#' @param cv Logical. Whether to perform Leave-One-Out Cross-Validation (default: \code{TRUE})
#' @param seed Integer. Random seed for reproducible permutation testing (default: \code{222})
#' @param plot Logical. Whether to generate and display the canonical separation plot (default: \code{TRUE})
#' @param save_plot Logical. Whether to save the generated plot to disk (default: \code{FALSE})
#' @param plot_dir Optional character. Directory to save the plot when \code{save_plot = TRUE}
#' @param palette Character. RColorBrewer palette name for group colors (default: \code{"Set1"})
#' @param ellipse Logical. Whether to draw confidence ellipses around group clusters (default: \code{TRUE})
#' @param ellipse_level Numeric. Confidence level for ellipses (default: \code{0.95})
#' @param alpha Numeric. Transparency for individual data points (default: \code{0.4})
#' @param point_size Numeric. Size of individual data points (default: \code{2})
#' @param verbose Logical. Whether to print analytical summary to console (default: \code{TRUE})
#' @param ... Additional arguments passed to internal plot helper
#'
#' @details
#' Canonical Discriminant Analysis finds linear combinations of features (Canonical Variates:
#' \code{CV1}, \code{CV2}, ...) that maximize between-group variance relative to within-group
#' variance.
#'
#' **Robustness and Validation:**
#' \itemize{
#'   \item \strong{Leave-One-Out Cross-Validation (LOOCV):} Evaluates unbiased classification
#'     accuracy without sample-size inflation.
#'   \item \strong{Permutation Test:} Randomly permutes group labels \code{n_perm} times to construct
#'     an empirical null distribution and compute a non-parametric permutation p-value for group
#'     discriminability.
#'   \item \strong{Mahalanobis Distance:} Computes multivariate distances between group centroids
#'     accounting for feature covariance.
#' }
#'
#' @return A list of class \code{"DiscriminantAnalysis"} containing:
#' \describe{
#'   \item{\code{scores}}{Data frame of projected canonical variate scores per observation.}
#'   \item{\code{centroids}}{Data frame of group centroids in canonical space.}
#'   \item{\code{loadings}}{Standardized canonical coefficients for each feature.}
#'   \item{\code{prop_trace}}{Proportion of between-group variance explained by each canonical axis.}
#'   \item{\code{mahalanobis}}{Matrix of pairwise Mahalanobis distances between group centroids.}
#'   \item{\code{reference_distances}}{Distances from each group to \code{reference_group} (if specified).}
#'   \item{\code{loocv_accuracy}}{Overall cross-validated classification accuracy.}
#'   \item{\code{confusion_matrix}}{Confusion matrix of observed versus LOOCV-predicted classes.}
#'   \item{\code{permutation}}{List containing \code{n_perm}, empirical \code{p_value}, and null distribution.}
#'   \item{\code{manova}}{Multivariate ANOVA (Wilks' Lambda) summary.}
#'   \item{\code{plot}}{The generated \code{ggplot} object (or \code{NULL} if \code{plot = FALSE}).}
#'   \item{\code{features_used}}{Character vector of feature column names used in the analysis.}
#' }
#'
#' @examples
#' \dontrun{
#' # Using trajectory similarity results
#' sim_data <- sap$features$motif$trajectory_similarity$similarity
#' res <- discriminant_analysis(
#'   data = sim_data,
#'   group_col = "label",
#'   reference_group = "tutor",
#'   n_perm = 1000
#' )
#'
#' # Inspect classification accuracy and confusion matrix
#' res$loocv_accuracy
#' res$confusion_matrix
#'
#' # Access or modify the plot
#' res$plot + ggplot2::theme_minimal()
#' }
#'
#' @seealso \code{\link{trajectory_similarity}}, \code{\link{feature_stats}}
#'
#' @export
discriminant_analysis <- function(data,
                                  group_col = "label",
                                  feature_cols = NULL,
                                  reference_group = NULL,
                                  scale = TRUE,
                                  n_perm = 1000,
                                  cv = TRUE,
                                  seed = 222,
                                  plot = TRUE,
                                  save_plot = FALSE,
                                  plot_dir = NULL,
                                  palette = "Set1",
                                  ellipse = TRUE,
                                  ellipse_level = 0.95,
                                  alpha = 0.4,
                                  point_size = 2,
                                  verbose = TRUE,
                                  ...) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }
  if (!group_col %in% names(data)) {
    stop(sprintf("Group column '%s' not found in data.", group_col))
  }

  # Automatically select numeric feature columns if not specified
  if (is.null(feature_cols)) {
    ignore_cols <- c(
      group_col, "rendition", ".source_row", ".time", "time", "n_matched_time",
      "n_time", "day", "ref_day", "song_id", "bout_id", "motif_id", "cluster"
    )
    feature_cols <- setdiff(names(data)[vapply(data, is.numeric, logical(1))], ignore_cols)
  }

  if (length(feature_cols) < 1) {
    stop("At least 1 numeric feature column is required for discriminant analysis.")
  }

  # Clean complete cases
  analysis_df <- data[, c(group_col, feature_cols), drop = FALSE]
  analysis_df <- stats::na.omit(analysis_df)
  groups <- as.factor(analysis_df[[group_col]])
  features <- as.matrix(analysis_df[, feature_cols, drop = FALSE])

  # Remove zero-variance features if any
  feat_vars <- apply(features, 2, stats::var)
  if (any(feat_vars == 0 | is.na(feat_vars))) {
    keep_feats <- names(feat_vars)[feat_vars > 0 & !is.na(feat_vars)]
    if (length(keep_feats) < 1) stop("No features have non-zero variance.")
    features <- features[, keep_feats, drop = FALSE]
    feature_cols <- keep_feats
  }

  # Scale features if requested
  if (scale) {
    features <- scale(features)
  }

  n_groups <- length(levels(groups))
  if (n_groups < 2) {
    stop("At least 2 distinct groups are required for discriminant analysis.")
  }

  # Fit LDA / CDA
  lda_fit <- MASS::lda(x = features, grouping = groups)

  # Proportion of between-group variance explained per CV
  prop_trace <- lda_fit$svd^2 / sum(lda_fit$svd^2)
  names(prop_trace) <- paste0("CV", seq_along(prop_trace))

  # Project observations onto Canonical Variates (CV scores)
  cv_scores <- as.data.frame(features %*% lda_fit$scaling)
  colnames(cv_scores) <- paste0("CV", seq_len(ncol(cv_scores)))
  cv_scores[[group_col]] <- groups

  # Standardized canonical coefficients
  feature_sds <- apply(analysis_df[, feature_cols, drop = FALSE], 2, stats::sd)
  std_loadings <- lda_fit$scaling * feature_sds

  # Centroids in Canonical Space
  centroids <- stats::aggregate(
    cv_scores[, grep("^CV", names(cv_scores)), drop = FALSE],
    by = list(Group = groups),
    FUN = mean
  )
  rownames(centroids) <- centroids$Group
  centroids$Group <- NULL

  # Pairwise Mahalanobis distances between group centroids
  cov_pooled <- stats::cov(features)
  inv_cov <- tryCatch(solve(cov_pooled), error = function(e) MASS::ginv(cov_pooled))
  orig_centroids <- stats::aggregate(features, by = list(Group = groups), FUN = mean)
  rownames(orig_centroids) <- orig_centroids$Group
  orig_centroids$Group <- NULL

  group_names <- levels(groups)
  mahalanobis_mat <- matrix(
    0, nrow = n_groups, ncol = n_groups,
    dimnames = list(group_names, group_names)
  )
  for (i in seq_len(n_groups)) {
    for (j in seq_len(n_groups)) {
      if (i != j) {
        diff_vec <- as.numeric(orig_centroids[i, ]) - as.numeric(orig_centroids[j, ])
        mahalanobis_mat[i, j] <- sqrt(t(diff_vec) %*% inv_cov %*% diff_vec)
      }
    }
  }

  # Distances to reference group
  ref_distances <- NULL
  if (!is.null(reference_group) && reference_group %in% group_names) {
    ref_distances <- sort(mahalanobis_mat[reference_group, ])
  }

  # Leave-One-Out Cross-Validation (LOOCV)
  loocv_acc <- NA_real_
  conf_mat <- NULL
  if (cv) {
    lda_cv <- MASS::lda(x = features, grouping = groups, CV = TRUE)
    conf_mat <- table(Observed = groups, Predicted = lda_cv$class)
    loocv_acc <- mean(groups == lda_cv$class)
  }

  # Permutation test
  perm_p_val <- NA_real_
  perm_accuracies <- numeric(0)
  if (n_perm > 0 && cv) {
    set.seed(seed)
    perm_accuracies <- vapply(seq_len(n_perm), function(p) {
      shuffled_groups <- sample(groups)
      cv_perm <- MASS::lda(x = features, grouping = shuffled_groups, CV = TRUE)
      mean(shuffled_groups == cv_perm$class)
    }, numeric(1))

    perm_p_val <- (sum(perm_accuracies >= loocv_acc) + 1) / (n_perm + 1)
  }

  # Multivariate ANOVA (or univariate ANOVA if 1 feature)
  manova_summary <- if (ncol(features) > 1) {
    stats::manova(features ~ groups) |> summary(test = "Wilks")
  } else {
    stats::aov(features[, 1] ~ groups) |> summary()
  }

  # Generate plot if requested
  p <- NULL
  if (plot || save_plot) {
    p <- render_cda_plot(
      scores = cv_scores,
      centroids = centroids,
      group_col = group_col,
      prop_trace = prop_trace,
      loocv_accuracy = loocv_acc,
      perm_p_value = perm_p_val,
      palette = palette,
      ellipse = ellipse,
      ellipse_level = ellipse_level,
      alpha = alpha,
      point_size = point_size
    )

    if (plot) {
      print(p)
    }

    if (save_plot) {
      if (is.null(plot_dir)) plot_dir <- getwd()
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
      plot_file <- file.path(plot_dir, "discriminant_analysis.png")
      ggplot2::ggsave(plot_file, plot = p, width = 8, height = 6, dpi = 300)
    }
  }

  res <- structure(
    list(
      scores = cv_scores,
      centroids = centroids,
      loadings = std_loadings,
      raw_scalings = lda_fit$scaling,
      prop_trace = prop_trace,
      mahalanobis = mahalanobis_mat,
      reference_distances = ref_distances,
      reference_group = reference_group,
      loocv_accuracy = loocv_acc,
      confusion_matrix = conf_mat,
      permutation = list(
        n_perm = n_perm,
        p_value = perm_p_val,
        null_distribution = perm_accuracies
      ),
      manova = manova_summary,
      plot = p,
      features_used = feature_cols,
      group_col = group_col,
      n_obs = nrow(analysis_df)
    ),
    class = "DiscriminantAnalysis"
  )

  if (verbose) {
    cat("=== Canonical Discriminant Analysis ===\n")
    cat(sprintf(
      "Observations: %d | Groups (%d): %s\n",
      res$n_obs, n_groups, paste(group_names, collapse = ", ")
    ))
    cat(sprintf("Features (%d): %s\n", length(feature_cols), paste(feature_cols, collapse = ", ")))
    cat("\nVariance Explained by Canonical Variates:\n")
    print(round(prop_trace * 100, 2))
    if (cv) {
      cat(sprintf("\nLOOCV Classification Accuracy: %.2f%%\n", loocv_acc * 100))
      if (n_perm > 0) {
        cat(sprintf("Permutation Test: empirical p = %.4f (%d permutations)\n", perm_p_val, n_perm))
      }
      cat("\nConfusion Matrix:\n")
      print(conf_mat)
    }
    if (!is.null(ref_distances)) {
      cat(sprintf("\nMahalanobis Distances to Reference ('%s'):\n", reference_group))
      print(round(ref_distances, 3))
    }
  }

  invisible(res)
}


#' Render Canonical Discriminant Analysis Plot
#'
#' @description
#' Internal helper used by \code{discriminant_analysis()} to construct a 2D or 1D
#' canonical separation ggplot visualization.
#'
#' @param scores Data frame of canonical scores
#' @param centroids Data frame of group centroids
#' @param group_col Character. Group column name
#' @param prop_trace Named numeric vector of variance explained per CV
#' @param loocv_accuracy Numeric. Cross-validated accuracy
#' @param perm_p_value Numeric. Permutation test p-value
#' @param palette Character. Palette name
#' @param ellipse Logical. Draw ellipses
#' @param ellipse_level Numeric. Ellipse confidence level
#' @param alpha Numeric. Point alpha
#' @param point_size Numeric. Point size
#'
#' @keywords internal
#' @noRd
render_cda_plot <- function(scores,
                            centroids,
                            group_col,
                            prop_trace,
                            loocv_accuracy,
                            perm_p_value,
                            palette = "Set1",
                            ellipse = TRUE,
                            ellipse_level = 0.95,
                            alpha = 0.4,
                            point_size = 2) {
  ensure_pkgs("ggplot2")

  groups <- scores[[group_col]]
  n_cvs <- sum(grepl("^CV", names(scores)))
  pal_map <- make_pal(levels(groups), palette)

  sub_text <- sprintf("LOOCV Accuracy: %.1f%%", loocv_accuracy * 100)
  if (!is.na(perm_p_value)) {
    sub_text <- paste0(sub_text, sprintf(" | Permutation p = %.4f", perm_p_value))
  }

  if (n_cvs >= 2) {
    # 2D Canonical space plot (3 or more groups)
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = .data$CV1,
        y = .data$CV2,
        color = .data[[group_col]],
        fill = .data[[group_col]]
      )
    ) +
      ggplot2::geom_point(alpha = alpha, size = point_size) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_fill_manual(values = pal_map) +
      ggplot2::labs(
        title = "Canonical Discriminant Analysis",
        subtitle = sub_text,
        x = sprintf("CV1 (%.1f%% between-group var)", prop_trace["CV1"] * 100),
        y = sprintf("CV2 (%.1f%% between-group var)", prop_trace["CV2"] * 100),
        color = "Group", fill = "Group"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      )

    if (ellipse) {
      p <- p + ggplot2::stat_ellipse(
        ggplot2::aes(x = .data$CV1, y = .data$CV2, group = .data[[group_col]]),
        geom = "polygon", alpha = 0.15, level = ellipse_level, show.legend = FALSE
      )
    }

    cent_df <- as.data.frame(centroids)
    cent_df[[group_col]] <- rownames(cent_df)
    p <- p + ggplot2::geom_point(
      data = cent_df,
      ggplot2::aes(x = .data$CV1, y = .data$CV2, fill = .data[[group_col]]),
      shape = 21, size = 4.5, stroke = 1.5, color = "black", inherit.aes = FALSE
    )
  } else {
    # 1D Canonical space plot (2 groups)
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(x = .data[[group_col]], y = .data$CV1, fill = .data[[group_col]])
    ) +
      ggplot2::geom_violin(alpha = 0.6) +
      ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5) +
      ggplot2::scale_fill_manual(values = pal_map) +
      ggplot2::labs(
        title = "Canonical Discriminant Analysis (1D)",
        subtitle = sub_text,
        y = "Canonical Variate 1 (CV1)", x = NULL, fill = "Group"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      )
  }

  return(p)
}
