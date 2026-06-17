# Compute Duration----------------------------------------------------------
# Update date : Mar. 20, 2025

#' Compute WAV File Durations
#'
#' @description
#' Calculates durations for all WAV files referenced in a SAP object's metadata
#' using parallel processing. Handles missing files gracefully by returning NA
#' durations and providing warnings.
#'
#' @param x A SAP object containing metadata and base path to audio files
#' @param cores Number of cores to use for parallel processing (NULL for auto-detection: total cores - 1)
#' @param verbose Logical flag to control progress messages and warnings (default: TRUE)
#'
#' @return Returns the modified SAP object with added duration column in metadata containing
#'         wave file durations in seconds. Missing files will have NA durations.
#'
#' @details
#' Key features:
#' - Parallel processing implementation using \code{parallel_apply}
#' - Automatic core detection with fallback to single-core processing
#' - Progress tracking and missing file warnings
#' - Preserves original object structure while adding duration information
#'
#' @export
compute_wav_durations <- function(x, cores = NULL, verbose = TRUE) {
  if(verbose) message(sprintf("\n=== Starting compute WAV file durations ===\n"))

  # Validate SAP structure
  if (!inherits(x, "Sap")) {
    stop("Input must be a SAP object")
  }
  if (is.null(x$base_path)) stop("sap$base_path must contain WAV directory path")
  if (is.null(x$metadata)) stop("sap$metadata is missing")

  # Set up parallel processing
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }

  metadata <- x$metadata
  base_path <- x$base_path

  # Main computation function
  # NOTE: construct_wav_path() internally calls stop() for missing files, so we
  # build the path directly here and wrap the entire read in tryCatch so that
  # any problem (missing file, corrupt WAV, etc.) returns NA_real_ instead of
  # crashing the parallel worker and producing an uncoercible list.
  process_row <- function(i) {
    row <- metadata[i, ]
    subdir <- if ("subfolder" %in% names(row)) row$subfolder else if ("day_post_hatch" %in% names(row)) row$day_post_hatch else NULL
    if (!is.null(subdir) && !is.na(subdir)) {
      wav_path <- file.path(base_path, subdir, row$filename)
    } else {
      wav_path <- file.path(base_path, row$filename)
    }

    if (!file.exists(wav_path)) {
      if (verbose) warning("File not found: ", wav_path)
      return(NA_real_)
    }

    tryCatch(
      as.numeric(seewave::duration(tuneR::readWave(wav_path))),
      error = function(e) {
        if (verbose) warning("Error reading ", wav_path, ": ", conditionMessage(e))
        NA_real_
      }
    )
  }

  # Execute in parallel
  durations <- parallel_apply(
    indices = seq_len(nrow(metadata)),
    FUN =  process_row,
    cores = cores
  )

  # Process results
  x$metadata$duration <- unlist(durations)

  # Reporting
  if (verbose) {
    valid_count <- sum(!is.na(x$metadata$duration))
    message("Computed durations for ", valid_count, "/", nrow(x$metadata), " files")
    if (any(is.na(x$metadata$duration))) {
      warning("Missing durations for ", sum(is.na(x$metadata$duration)), " files")
    }
  }

  invisible(x)
}


# Compute Motif Duration----------------------------------------------------------
# Update date : Mar. 20, 2025

#' Refine Motif Boundaries Using Segment Alignment
#'
#' @description
#' Adjusts motif boundaries based on underlying segments and optional label-specific
#' time adjustments. Integrates segment information to determine precise motif
#' onsets/offsets.
#'
#' @param x A SAP object containing 'segments' and 'motifs' data
#' @param adjustments_by_label Optional named list of time adjustments (in seconds)
#'        to apply to motif limits, where names correspond to labels. Each adjustment
#'        can be a single numeric value (adjusts only the end limit) or a numeric vector
#'        of length 2 (first adjusts the start limit, second adjusts the end limit).
#' @param sample_percent Percentage of motifs to refine from each label (0-100)
#' @param balanced Logical indicating whether to balance motif counts across labels
#' @param max_samples_per_label Optional integer. Maximum number of motifs to
#'   randomly refine from each label. If a label has fewer motifs, all available
#'   motifs are used. Unselected motifs are retained with missing boundary values.
#' @param labels Character vector of specific labels to refine
#' @param clusters Numeric vector of cluster IDs to refine. Requires motif feature
#'   embeddings from prior clustering.
#' @param seed Random seed for reproducible sampling (default: 222)
#' @param verbose Logical flag for printing progress messages (default: TRUE)
#'
#' @return Returns a modified SAP object with updated motifs containing:
#' \itemize{
#'   \item motif_onset - Precise start time based on segments
#'   \item motif_offset - Precise end time based on segments
#'   \item first_seg_index - First segment index in motif
#'   \item last_seg_index - Last segment index in motif
#'   \item motif_duration - Calculated motif duration
#' }
#'
#' @details
#' Key operations:
#' 1. Applies label-specific time adjustments to motif start and end limits
#' 2. Identifies segments contained within adjusted motif boundaries
#' 3. Calculates precise motif timing based on contained segments
#' 4. Preserves original motif structure while adding new timing columns
#'
#' @examples
#' \dontrun{
#' # Apply 0.1s extension to "BL" motif end limits
#' sap <- refine_motif_boundaries(
#'   sap,
#'   adjustments_by_label = list(BL = 0.1)
#' )
#'
#' # Adjust both start (-0.05s) and end (+0.1s) limits for "Rec" motifs
#' sap <- refine_motif_boundaries(
#'   sap,
#'   adjustments_by_label = list(Rec = c(-0.05, 0.1))
#' )
#' }
#'
#' @export
refine_motif_boundaries <- function(x,
                                    adjustments_by_label = NULL,
                                    sample_percent = NULL,
                                    balanced = FALSE,
                                    max_samples_per_label = NULL,
                                    labels = NULL,
                                    clusters = NULL,
                                    seed = 222,
                                    verbose = TRUE) {

  # Remove previous boundary calculations if they exist
  motifs <- x$motifs %>%
    dplyr::select(-dplyr::any_of(c("motif_onset", "motif_offset",
                                   "first_seg_index", "last_seg_index",
                                   "motif_duration")))

  # Extract segments and motifs
  segments <- x$segments

  motifs <-  motifs |>
    dplyr::group_by(filename) |>
    dplyr::mutate(motif_index = dplyr::row_number()) |>
    dplyr::ungroup()

  select_motifs <- balanced ||
    !is.null(sample_percent) ||
    !is.null(max_samples_per_label) ||
    !is.null(labels) ||
    !is.null(clusters)

  if (select_motifs) {
    motifs_for_selection <- motifs

    if (!is.null(clusters)) {
      feature_data <- x$features$motif$feat.embeds
      if (is.null(feature_data)) {
        stop("Feature embeddings required for cluster filtering")
      }
      if (!"cluster" %in% names(feature_data)) {
        stop("No 'cluster' column found in motif feature embeddings")
      }

      key_cols <- c("filename", "start_time", "end_time", "label", "day_post_hatch")
      key_cols <- intersect(key_cols, intersect(names(motifs), names(feature_data)))
      embed_cols <- unique(c(key_cols, "cluster"))
      embed_cols <- intersect(embed_cols, names(feature_data))

      motifs_for_selection <- motifs |>
        dplyr::inner_join(
          feature_data[, embed_cols],
          by = key_cols,
          relationship = "many-to-many"
        )
    }

    motifs_join_keys <- motifs_for_selection |>
      select_segments(
        labels = labels,
        clusters = clusters,
        balanced = balanced,
        sample_percent = sample_percent,
        max_samples_per_label = max_samples_per_label,
        seed = seed
      ) |>
      dplyr::select(dplyr::all_of(c("filename", "motif_index"))) |>
      dplyr::distinct()

    motifs_to_refine <- motifs |>
      dplyr::inner_join(motifs_join_keys, by = c("filename", "motif_index"))
  } else {
    motifs_to_refine <- motifs
  }

  # Create temporary motifs with limit columns
  motifs_join <- motifs_to_refine |>
    dplyr::rename(start_limit = start_time, end_limit = end_time)

  # Add label-based limit adjustments
  if(!is.null(adjustments_by_label)) {
    # Validate adjustments
    if(!is.list(adjustments_by_label) ||
       !all(names(adjustments_by_label) %in% unique(motifs$label))) {
      invalid_labels <- setdiff(names(adjustments_by_label), unique(motifs$label))
      stop("Invalid labels in adjustments: ", paste(invalid_labels, collapse = ", "))
    }

    # Parse and validate each adjustment value (numeric of length 1 or 2)
    parsed_adjs <- lapply(adjustments_by_label, function(adj) {
      if (!is.numeric(adj) || length(adj) < 1 || length(adj) > 2) {
        stop("Adjustments must be numeric vectors of length 1 or 2")
      }
      if (length(adj) == 1) {
        return(c(start = 0, end = adj))
      } else {
        if (!is.null(names(adj))) {
          start_val <- if ("start" %in% names(adj)) adj[["start"]] else adj[1]
          end_val <- if ("end" %in% names(adj)) adj[["end"]] else adj[2]
          return(c(start = start_val, end = end_val))
        }
        return(c(start = adj[1], end = adj[2]))
      }
    })

    # Apply adjustments
    ensure_pkgs("purrr")
    motifs_join <- motifs_join |>
      dplyr::mutate(
        start_limit = .data$start_limit + purrr::map_dbl(
          label,
          ~ if (.x %in% names(parsed_adjs)) parsed_adjs[[.x]]["start"] else 0
        ),
        end_limit = .data$end_limit + purrr::map_dbl(
          label,
          ~ if (.x %in% names(parsed_adjs)) parsed_adjs[[.x]]["end"] else 0
        )
      )

    if(verbose) {
      message("Applied limit adjustments by label:")
      print_df <- data.frame(
        Label = names(adjustments_by_label),
        Start_Adjustment = vapply(parsed_adjs, function(x) x["start"], numeric(1)),
        End_Adjustment = vapply(parsed_adjs, function(x) x["end"], numeric(1)),
        stringsAsFactors = FALSE
      )
      print(print_df, row.names = FALSE)
    }
  }

  # Find segments within motifs
  segments_in_motifs <- segments |>
    dplyr::inner_join(motifs_join, by = c("filename", "day_post_hatch", "label"),
                      relationship = "many-to-many") |>
    dplyr::filter(start_time >= .data$start_limit & start_time <= .data$end_limit)

  # Aggregate to find motif boundaries
  agg_results <- segments_in_motifs |>
    dplyr::group_by(filename, .data$motif_index) |>
    dplyr::summarize(
      motif_onset = min(start_time),
      motif_offset = max(end_time),
      first_seg_index = dplyr::first(.data$selec[start_time == .data$motif_onset]),
      last_seg_index = dplyr::last(.data$selec[end_time == .data$motif_offset]),
      .groups = "drop"
    ) |>
    dplyr::mutate(motif_duration = .data$motif_offset - .data$motif_onset)

  # Merge results back with original motifs
  updated_motifs <- motifs |>
    dplyr::left_join(agg_results, by = c("filename", "motif_index"))

  # Convert to segment object and update SAP
  x$motifs <- as_segment(updated_motifs)

  invisible(x)
}

#' Visualize Motif Boundaries in Audio Data
#'
#' @description
#' Generates a heatmap visualization of motif boundaries with optional cluster filtering
#' and UMAP-based ordering. Incorporates amplitude envelopes and boundary markers for
#' precise temporal analysis.
#'
#' @param x A SAP object containing motifs and metadata
#' @param sample_percent Percentage of data to sample from each group (0-100)
#' @param balanced Logical indicating whether to balance group sizes
#' @param max_samples_per_label Optional integer. Maximum number of samples to
#'   randomly analyze from each label when \code{balanced = FALSE}. If a label
#'   has fewer samples, all available samples are used. Ignored when
#'   \code{balanced = TRUE}.
#' @param labels Character vector of specific labels to include
#' @param clusters Numeric vector of cluster IDs to filter
#' @param ordered Logical indicating UMAP-based ordering
#' @param descending Logical for descending order in UMAP sorting
#' @param cores Number of processing cores (NULL for auto-detection)
#' @param seed Random seed for reproducibility
#' @param msmooth Vector of smoothing parameters for amplitude envelopes
#' @param color_palette Color palette function for heatmap
#' @param n_colors Number of color gradations in palette
#' @param contrast Contrast adjustment for color scaling
#' @param marginal_window Time window extension for motif margins (seconds)
#' @param verbose Logical flag for progress messages
#' @param ... Additional parameters passed to lattice::levelplot
#'
#' @return A lattice plot object showing motif boundaries heatmap
#'
#' @details
#' Key features:
#' - Integrates motif boundary data from refine_motif_boundaries()
#' - Supports cluster filtering using feature embeddings
#' - Enables UMAP-based ordering of motifs
#' - Visualizes boundaries with colored markers (green = onset, cyan = offset)
#' - Includes automatic duration calculation and marginal window extension
#'
#' Requires prior execution of refine_motif_boundaries() for boundary detection.
#' Uses parallel processing for efficient amplitude envelope calculation.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' plot_motif_boundaries(sap_object)
#'
#' # Cluster-specific plot with custom colors
#' plot_motif_boundaries(sap_object,
#'                       clusters = c(1,3),
#'                       ordered = TRUE)
#' }
#'
#' @export
plot_motif_boundaries <- function(x,
                                  sample_percent = NULL,
                                  balanced = FALSE,
                                  max_samples_per_label = NULL,
                                  labels = NULL,
                                  clusters = NULL,
                                  ordered = FALSE,
                                  descending = TRUE,
                                  cores = NULL,
                                  seed = 222,
                                  msmooth = c(256,50),
                                  color_palette = NULL,
                                  n_colors = 500,
                                  contrast = 3,
                                  marginal_window = 0.1,
                                  verbose = TRUE,
                                  ...) {
  if(verbose) message(sprintf("\n=== Starting Plotting Motif Boundaries===\n"))

  # Handle feature embeddings for clustering/ordering
  if (!is.null(clusters) | ordered) {
    if (is.null(x$features$motif$feat.embeds)) {
      stop("Feature embeddings required for cluster filtering or ordering")
    }

    # Define natural composite key columns
    key_cols <- c("filename", "start_time", "end_time", "label", "day_post_hatch")

    # Merge using natural keys
    segments_df <- x$motifs %>%
      dplyr::inner_join(
        x$features$motif$feat.embeds,
        by = key_cols,
        relationship = "one-to-one"  # Ensures 1:1 mapping between motifs and embeddings
      )
  } else {
    segments_df <- x$motifs
  }

  # Apply UMAP-based ordering if requested
  if(ordered) {
    if(!is.null(x$features$motif$feat.embeds) &&
       all(c("UMAP1", "UMAP2") %in% colnames(segments_df))) {
      segments_df <- segments_df |>
        dplyr::arrange(day_post_hatch,
                       if(descending) dplyr::desc(UMAP2) else UMAP2,
                       if(descending) dplyr::desc(UMAP1) else UMAP1)
    } else {
      warning("Ordering requested but UMAP coordinates not found - using original order")
    }
  }

  # Check for required columns
  if(!all(c("motif_onset", "motif_offset") %in% colnames(segments_df))) {
    stop("motif_onset and motif_offset columns must exist. Run refine_motif_boundaries() first.")
  }

  # Filter and report NA values
  original_count <- nrow(segments_df)
  segments_df <- segments_df %>%
    dplyr::filter(!is.na(.data$motif_onset), !is.na(.data$motif_offset))

  removed_count <- original_count - nrow(segments_df)

  if(verbose && removed_count > 0) {
    message(sprintf("Removed %d motifs without identified boundaries (%.1f%%)",
                    removed_count,
                    (removed_count/original_count)*100))
  }

  # Select and balance segments
  segments_df <- select_segments(segments_df,
                                  labels = labels,
                                  clusters = clusters,
                                  balanced = balanced,
                                  sample_percent = sample_percent,
                                  max_samples_per_label = max_samples_per_label,
                                  seed = seed)
  # Ensure WAV durations are computed
  if (!"duration" %in% names(x$metadata)) {
    x <- compute_wav_durations(x, cores = cores, verbose = verbose)
  }

  # add durations of WAV file
  segments_df <- segments_df |>
    dplyr::left_join(x$metadata |> dplyr::select(filename, wav_duration = duration),
              by = "filename")

  # Add marginal window
  segments_df <- add_marginal_window(segments_df, marginal_window)

  # Function to process a single row
  process_row <- function(i) {
    env <- amp_env(segments_df[i,],
                   wav_dir = x$base_path,
                   msmooth = msmooth)

    frame_indices <- compute_motif_frame_indices(segments_df[i,], length(env))

    list(
      envelope = env,
      onset = frame_indices$onset,
      offset = frame_indices$offset
    )
  }

  # Calculate time window based on the most frequent duration
  time_window <- names(table(segments_df$duration)[which.max(table(segments_df$duration))])
  time_window <- as.numeric(time_window)


  # Generate amplitude envelope matrix using parallel processing
  results <- parallel_apply(
    indices = seq_len(nrow(segments_df)),
    FUN = process_row,
    cores = cores
  )

  # Create amplitude matrix
  amp_matrix <- do.call(cbind, lapply(results, function(x) x$envelope))

  # Set column names as labels
  colnames(amp_matrix) <- segments_df$label

  # Create reversed matrix for plotting
  reversed_amp_matrix <- amp_matrix[, ncol(amp_matrix):1]

  # Get reversed labels and calculate positions
  reversed_labels <- rev(unique(colnames(amp_matrix)))
  samples_per_label <- rev(table(colnames(amp_matrix)))
  cumulative_positions <- cumsum(c(0, head(samples_per_label, -1)))
  label_positions <- cumulative_positions + samples_per_label/2
  hline_positions <- cumsum(samples_per_label)[-length(samples_per_label)]

  # Set color palette
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("black", "red", "yellow", "white"))
  }

  # Reverse onsets and offsets to match matrix orientation
  onsets <- rev(sapply(results, function(x) x$onset))
  offsets <- rev(sapply(results, function(x) x$offset))

  # Create heatmap
  heatmap <- lattice::levelplot(reversed_amp_matrix,
                                col.regions = color_palette(n_colors),
                                at = seq(min(reversed_amp_matrix),
                                         max(reversed_amp_matrix)/contrast,
                                         length.out = n_colors),
                                border = "transparent",
                                ylab = "Labels",
                                xlab = "Time (s)",
                                main = "Heatmap of motif boundaries",
                                scales = list(
                                  x = list(
                                    at = seq(0, nrow(reversed_amp_matrix), length.out = 5),
                                    labels = sprintf("%.1f", seq(0, time_window, length.out = 5))
                                  ),
                                  y = list(
                                    at = label_positions,
                                    labels = reversed_labels
                                  )
                                ),
                                panel = function(...) {
                                  lattice::panel.levelplot(...)

                                  # Add onset and offset lines
                                  for(i in 1:length(onsets)) {
                                    if(!is.na(onsets[i])) {
                                      lattice::panel.segments(y0 = i-0.5, y1 = i+0.5,
                                                              x0 = onsets[i], x1 = onsets[i],  # same y value for horizontal line
                                                              col = "green", lwd = 2)
                                    }
                                    if(!is.na(offsets[i])) {
                                      lattice::panel.segments(y0 = i-0.5, y1 = i+0.5,
                                                              x0 = offsets[i], x1 = offsets[i],  # same y value for horizontal line
                                                              col = "cyan2", lwd = 2)
                                    }
                                  }

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
                                    at = seq(min(amp_matrix),
                                             max(amp_matrix) / contrast,
                                             length.out = 3),
                                    cex = 1,
                                    col = "black",
                                    labels = sprintf("%.1f", seq(0, 1, length.out = 3))
                                  )
                                ))


  print(heatmap)

  invisible(x)
}

#' Adjust Segment Time Windows with Safety Margins
#'
#' @description
#' Internal function that expands segment time boundaries by a specified margin while
#' ensuring the adjusted times stay within audio file duration limits.
#'
#' @param segments_df Dataframe of audio segments with timing information
#' @param marginal_window Time margin (in seconds) to add to both ends of segments
#' @param verbose Logical flag for printing processing summary
#'
#' @return Returns modified dataframe with:
#' \itemize{
#'   \item Adjusted start/end times incorporating margins
#'   \item Removed invalid segments (if any)
#'   \item Updated duration calculations
#' }
#'
#' @details
#' Key operations:
#' 1. Adds specified margin to both ends of each segment
#' 2. Filters out segments that would:
#'    - Start before 0 seconds
#'    - End after file duration
#'    - Have missing duration data
#' 3. Provides detailed removal statistics when verbose=TRUE
#'
#' @examples
#' \dontrun{
#' segments <- add_marginal_window(segments, 0.1)
#' }
#'
#' @keywords internal
#' @noRd
add_marginal_window <- function(segments_df, marginal_window = 0.1, verbose = TRUE) {
  if (!is.null(marginal_window)) {
    # Input validation
    if (!is.numeric(marginal_window) || length(marginal_window) != 1) {
      stop("marginal_window must be a single numeric value")
    }
    if (!"wav_duration" %in% names(segments_df)) {
      stop("segments_df must contain wav_duration column")
    }

    original_count <- nrow(segments_df)

    # Create a copy to work with
    processed <- segments_df

    # Create adjusted time columns
    processed$adj_start <- processed$start_time - marginal_window
    processed$adj_end <- processed$end_time + marginal_window

    # Calculate removal reasons BEFORE filtering
    removal_stats <- list(
      negative_start = sum(processed$adj_start < 0, na.rm = TRUE),
      exceeds_duration = sum(processed$adj_end > processed$wav_duration, na.rm = TRUE),
      missing_duration = sum(is.na(processed$wav_duration))
    )

    # Now filter
    keep_idx <- processed$adj_start >= 0 &
      processed$adj_end <= processed$wav_duration &
      !is.na(processed$wav_duration)
    processed <- processed[keep_idx, ]

    # Update time columns
    processed$start_time <- processed$adj_start
    processed$end_time <- processed$adj_end
    processed$duration <- processed$end_time - processed$start_time
    processed$adj_start <- NULL
    processed$adj_end <- NULL

    removal_stats$total <- original_count - nrow(processed)
  }

  # Generate report
  if (verbose) {
    message(sprintf(
      paste0(
        "\nAdd marginal window (%.2fs) results:\n",
        "  Original segments: %d\n",
        "  Removed segments: %d (%.1f%%)\n",
        "  Remaining segments: %d\n",
        "Removal reasons:\n",
        "  - Negative start time: %d\n",
        "  - Exceeds file duration: %d\n",
        "  - Missing duration data: %d"
      ),
      marginal_window,
      original_count,
      removal_stats$total,
      (removal_stats$total/original_count)*100,
      nrow(processed),
      removal_stats$negative_start,
      removal_stats$exceeds_duration,
      removal_stats$missing_duration
    ))
  }

  return(processed)
}

#' Calculate Motif Boundary Frame Indices
#'
#' @description
#' Internal function that converts motif timing information to amplitude envelope
#' frame indices for precise visualization alignment.
#'
#' @param segment_row Single-row dataframe containing segment timing information
#' @param env_length Length of the amplitude envelope vector
#'
#' @return List with:
#' \itemize{
#'   \item onset - Frame index for motif start
#'   \item offset - Frame index for motif end
#' }
#'
#' @details
#' Key operations:
#' 1. Converts absolute timing to relative positions within segment
#' 2. Maps relative times to envelope vector indices
#' 3. Ensures indices stay within valid range
#'
#' @examples
#' \dontrun{
#' indices <- compute_motif_frame_indices(segment[1,], 500)
#' }
#'
#' @keywords internal
#' @noRd
compute_motif_frame_indices <- function(segment_row, env_length) {
  # Validate inputs
  if (!is.data.frame(segment_row) || nrow(segment_row) != 1) {
    stop("segment_row must be a single row from a data frame")
  }

  # Check for required timing columns
  required_cols <- c("start_time", "end_time", "motif_onset", "motif_offset")
  if (!all(required_cols %in% names(segment_row))) {
    stop("Input row must contain 'start_time', 'end_time', 'motif_onset', and 'motif_offset' columns")
  }

  # Calculate total segment duration
  duration <- segment_row$duration

  # Compute relative times for motif onset and offset
  relative_motif_onset <- segment_row$motif_onset - segment_row$start_time
  relative_motif_offset <- segment_row$motif_offset - segment_row$start_time

  # Convert relative times to frame indices
  onset_frame_index <- round(relative_motif_onset / duration * env_length)
  offset_frame_index <- round(relative_motif_offset / duration * env_length)

  # Ensure indices are within bounds
  onset_frame_index <- max(1, min(onset_frame_index, env_length))
  offset_frame_index <- max(1, min(offset_frame_index, env_length))

  return(list(
    onset = onset_frame_index,
    offset = offset_frame_index
  ))
}
