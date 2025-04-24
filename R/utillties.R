# utilities

#' @importClassesFrom tuneR Wave
#' @importFrom tuneR readWave
#' @importFrom methods getMethod
#' @importFrom rlang `%||%` .data
#' @importFrom tools file_path_sans_ext
#' @importFrom grDevices colorRampPalette dev.off png rgb
#' @importFrom utils head tail
#' @importFrom graphics abline axis box image layout legend lines mtext par
#'             plot.new plot.window points rect text title
#' @importFrom dplyr mutate select filter arrange left_join bind_rows %>%
#'             group_by n_distinct group_split case_when
#' @importFrom ggplot2 ggplot aes geom_line geom_boxplot facet_wrap labs
#'             theme_minimal scale_color_brewer scale_fill_brewer ggtitle
#'             theme element_text stat_summary
#' @importFrom stats TukeyHSD aggregate aov approx as.formula dist gaussian
#'             median na.omit prcomp quantile setNames sd time
#' @importFrom seewave spec inputw ftwindow sfm sh th meanspec afilter sspectro
#'
NULL

# Define global variables used in NSE contexts
utils::globalVariables(c(
  "day_post_hatch",
  "dph",
  "UMAP1",
  "UMAP2",
  "start_time",
  "end_time",
  "label",
  "filename"
))

#' Construct file path for audio file
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

  # Construct file path based on available information
  if (!is.null(attr(x, "wav_dir"))) {
    if ("day_post_hatch" %in% names(x)) {
      sound_path <- file.path(attr(x, "wav_dir"), x$day_post_hatch, x$filename)
    } else {
      sound_path <- file.path(attr(x, "wav_dir"), x$filename)
    }
  } else if (!is.null(wav_dir)) {
    if ("day_post_hatch" %in% names(x)) {
      sound_path <- file.path(wav_dir, x$day_post_hatch, x$filename)
    } else {
      sound_path <- file.path(wav_dir, x$filename)
    }
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
      stop(sprintf("Clusters not found: %s\nAvailable clusters: %s",
                   paste(missing_clusters, collapse = ", "),
                   paste(available_clusters, collapse = ", ")))
    }

    segments_df <- segments_df[segments_df$cluster %in% clusters, ]
    message(sprintf("\nFiltered for clusters: %s", paste(clusters, collapse = ", ")))
  }

  # Handle label selection
  if (!is.null(labels)) {
    available_labels <- unique(segments_df$label)
    missing_labels <- labels[!labels %in% available_labels]
    if (length(missing_labels) > 0) {
      stop(sprintf("Labels not found: %s\nAvailable labels: %s",
                   paste(missing_labels, collapse = ", "),
                   paste(available_labels, collapse = ", ")))
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
    group_counts <- segments_df |>
      dplyr::count(label)
    n_sample <- floor(min(group_counts$n) / 10) * 10  # Round to nearest 10

    message(sprintf("\nBalancing groups to %d segments each", n_sample))

    segments_df <- segments_df |>
      dplyr::group_by(label) |>
      dplyr::slice_sample(n = n_sample) |>
      dplyr::arrange(.original_order) |>
      dplyr::ungroup() |>
      dplyr::arrange(day_post_hatch, .original_order)
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
        dplyr::slice_sample(.x, n = n_to_sample)  |>
          dplyr::arrange(.original_order)
      }) |>
      dplyr::ungroup() |>
      dplyr::arrange(day_post_hatch)

    message(sprintf("\nSampling %.1f%% from each label", sample_percent))
  }

  # Final cleanup and ordering
  segments_df <- segments_df |>
    dplyr::arrange(.original_order) |>  # Global order preservation
    dplyr::select(-.original_order)

  # Print final summary if data was modified
  if (balanced || !is.null(sample_percent)) {
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
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#'
#' @keywords internal
parallel_apply <- function(indices, FUN, cores) {
  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
    cores <- max(1, cores)
  }

  if (cores > 1) {
    if (Sys.info()["sysname"] == "Darwin") {
      # macOS/Linux: Use pbmcapply with multicore
      result <- pbmcapply::pbmclapply(
        indices,
        FUN,
        mc.cores = cores,
        mc.preschedule = FALSE
      )
    } else {
      # Windows: Use pbapply with PSOCK cluster
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      result <- pbapply::pblapply(
        indices,
        FUN,
        cl = cl
      )
    }
  } else {
    # Single-core with progress
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
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required. Please install it using: install.packages('reticulate')")
  }

  # Check for librosa
  tryCatch({
    reticulate::py_module_available("librosa")
  }, error = function(e) {
    if (verbose) {
      message("Librosa is not installed. Attempting to install...")
    }

    # Try to install librosa
    tryCatch({
      reticulate::py_install("librosa", pip = TRUE)
    }, error = function(e) {
      stop("Failed to install librosa. Please install it manually using pip: pip install librosa")
    })
  })

  # Check for numpy
  tryCatch({
    reticulate::py_module_available("numpy")
  }, error = function(e) {
    if (verbose) {
      message("Numpy is not installed. Attempting to install...")
    }

    # Try to install numpy
    tryCatch({
      reticulate::py_install("numpy", pip = TRUE)
    }, error = function(e) {
      stop("Failed to install numpy. Please install it manually using pip: pip install numpy")
    })
  })
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
      process_label(label, feature_matrix, segment_starts, segment_ends,
                    valid_segments, time_step)
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
  if (length(label_cols) == 0) return(NULL)

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
      } else NULL
    } else NULL
  }

  # Use parallel_apply to process combinations
  parallel_apply(
    1:nrow(combinations),
    function(idx) process_combination(combinations[idx,]),
    cores = NULL  # Will use default core detection
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
    dplyr::filter(dplyr::n_distinct(.data$label) > 1) |>  # Only compare segments with multiple labels
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
  for(result in results_list) {
    tukey_df <- result$tukey
    for(i in 1:nrow(tukey_df)) {
      cat(sprintf("\n%-10s %5.0f %12.4f %12.4f %12.4f %12.5f    %s",
                  tukey_df$comparison[i],
                  tukey_df$segment_id[i],
                  tukey_df$diff[i],
                  tukey_df$lwr[i],
                  tukey_df$upr[i],
                  tukey_df$`p adj`[i],
                  tukey_df$significance[i]))
    }
  }

  # Add significance code explanation
  cat("\n\nSignificance codes (adjusted p-values):")
  cat("\n'***' < 0.001")
  cat("\n'**'  < 0.01")
  cat("\n'*'   < 0.05")
  cat("\n'ns'  >= 0.05\n")

  # Create plot if requested
  if(plot) {
    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = .data$label, y = .data$mean)) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~ .data$segment_id) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Mean Values by Label across Segments",
                    x = "Label",
                    y = "Mean Value")
    print(p)
  }

  # Return ANOVA summary
  return(anova_summary)
}
