# Song Trajectory Analysis
# Update date : Feb. 7, 2025

#' Create Spectrogram Matrices for Song Trajectory Analysis
#'
#' @description
#' Creates spectrogram matrices from sliding windows for song trajectory analysis.
#'
#' @param x An object to process, either a data frame or SAP object
#' @param wav_dir Directory containing WAV files (for default method)
#' @param window_size Size of sliding window in seconds (default: 0.1)
#' @param step_size Step size between windows (default: 0.005)
#' @param wl Window length for spectrogram (default: 128)
#' @param ovlp Overlap percentage (default: 50)
#' @param fftw Logical, use FFTW or not (default: TRUE)
#' @param flim Frequency limits (default: c(1, 12))
#' @param cores Number of processing cores
#' @param segment_type For SAP objects: Type of segments ('motifs', 'syllables', 'bouts', 'segments')
#' @param data_type For SAP objects: Type of data to analyze
#' @param clusters For SAP objects: Specific clusters to include
#' @param sample_percent For SAP objects: Percentage to sample
#' @param balanced For SAP objects: Whether to balance across groups
#' @param labels For SAP objects: Specific labels to include
#' @param seed For SAP objects: Random seed
#' @param verbose For SAP objects: Whether to print progress
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' Creates trajectory matrix with the following steps:
#' \itemize{
#'   \item Generates sliding windows for each segment
#'   \item Computes spectrograms for each window
#'   \item Combines results into matrix form
#' }
#'
#' For SAP objects, additional features include:
#' \itemize{
#'   \item Support for different segment types
#'   \item Optional cluster/label filtering
#'   \item Balanced sampling options
#'   \item Results storage in features slot
#' }
#'
#' @return
#' For default method: A list containing:
#' \itemize{
#'   \item spectrogram_matrix: Matrix of spectrogram vectors
#'   \item sliding_windows: Data frame of window information
#' }
#'
#' For SAP objects: Updated SAP object with trajectory matrix stored in features slot
#'
#' @examples
#' \dontrun{
#' # Create trajectory matrix from segments
#' matrix <- create_trajectory_matrix(segments,
#'                                   wav_dir = "path/to/wavs",
#'                                   window_size = 0.1,
#'                                   step_size = 0.005)
#'
#' # Create trajectory matrix from SAP object
#' sap_obj <- create_trajectory_matrix(sap_object,
#'                                    segment_type = "motifs",
#'                                    balanced = TRUE,
#'                                    sample_percent = 80)
#'
#' # Create matrix with specific clusters
#' sap_obj <- create_trajectory_matrix(sap_object,
#'                                    segment_type = "syllables",
#'                                    clusters = c(1, 2),
#'                                    labels = c("a", "b"))
#' }
#'
#' @rdname create_trajectory_matrix
#' @export
create_trajectory_matrix <- function(x, ...) {
  UseMethod("create_trajectory_matrix")
}

#' @rdname create_trajectory_matrix
#' @export
create_trajectory_matrix.default <- function(
    x,
    wav_dir,
    window_size = 0.1,
    step_size = 0.005,
    wl = 128,
    ovlp = 50,
    fftw = TRUE,
    flim = c(1, 12),
    cores = NULL,
    ...
) {
  # Input validation
  if (missing(x) || missing(wav_dir)) {
    stop("Both data frame containing audio segment information and wave file directory must be provided")
  }
  if (!is.data.frame(x)) {
    stop("Input x must be a data frame")
  }
  required_cols <- c("filename", "day_post_hatch", "label", "start_time", "end_time")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  if (!dir.exists(wav_dir)) {
    stop("Wave file directory does not exist: ", wav_dir)
  }

  # Set number of cores
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }
  sys_type <- Sys.info()["sysname"]

  # Generate constrained sliding windows
  sliding_windows <- do.call(rbind, lapply(
    seq_len(nrow(x)),
    function(i) create_sliding_window(i, x, window_size, step_size)
  ))

  # Split windows by original segment
  segments <- split(sliding_windows, sliding_windows$rendition)

  # Define core processing function
  process_segment <- function(segment) {
    tryCatch({
      original_row <- x[segment$rendition[1], ]
      sound_path <- construct_wav_path(original_row, wav_dir = wav_dir)

      # Read audio with small buffer
      buffer <- max(window_size * 0.005, step_size)
      wave <- tuneR::readWave(
        sound_path,
        from = original_row$start_time,
        to = original_row$end_time + buffer,
        units = "seconds"
      )

      # Compute spectrogram
      sr <- wave@samp.rate
      if (fftw) ensure_pkgs("fftw")

      spg <- seewave::spectro(
        wave,
        f = sr,
        wl = wl,
        ovlp = ovlp,
        fftw = fftw,
        flim = flim,
        plot = FALSE,
        osc = FALSE,
        cont = FALSE
      )

      # Calculate time bin parameters
      time_vector <- spg$time
      time_per_bin <- time_vector[2] - time_vector[1]
      expected_bins <- floor(window_size / time_per_bin)

      # Process all windows in segment
      vectors <- lapply(seq_len(nrow(segment)), function(j) {
        sw <- segment[j, ]
        start_rel <- sw$.time
        end_rel <- start_rel + window_size

        # Find matching time bins
        idx <- which(time_vector >= start_rel & time_vector <= end_rel)

        # Extract exactly expected_bins bins
        idx <- idx[1:expected_bins]

        # Convert to vector and return
        as.vector(spg$amp[, idx, drop = FALSE])
      })

      # Combine vectors column-wise
      do.call(cbind, vectors)
    }, error = function(e) {
      warning(sprintf("Error processing segment %d: %s",
                      segment$rendition[1], e$message))
      return(NULL)
    })
  }

  # Parallel processing setup
  message(sprintf("Generating spectrograms using %d cores...", cores))
  if (fftw) ensure_pkgs("fftw")
  spc_list <- parallel_apply(segments, process_segment, cores)
  # if (cores > 1) {
  #   if (sys_type == "Linux") {
  #     requireNamespace("pbmcapply", quietly = TRUE)
  #     spc_list <- pbmcapply::pbmclapply(
  #       segments,
  #       process_segment,
  #       mc.cores = cores
  #     )
  #   } else {
  #     # cl <- parallel::makeCluster(cores)
  #     # on.exit(parallel::stopCluster(cl))
  #     # parallel::clusterExport(
  #     #   cl,
  #     #   c("construct_wav_path", "create_sliding_window"),
  #     #   envir = environment()
  #     # )
  #     spc_list <- pbapply::pblapply(
  #       segments,
  #       process_segment,
  #       cl = cores
  #     )
  #   }
  # } else {
  #   spc_list <- pbapply::pblapply(segments, process_segment)
  # }

  # Combine results
  spc_list <- spc_list[!sapply(spc_list, is.null)]
  if (length(spc_list) == 0) {
    stop("No valid spectrograms were generated")
  }
  spectrogram_matrix <- do.call(cbind, spc_list) |> t()

  cat("Spectrogram generation complete!\n")
  list(
    spectrogram_matrix = spectrogram_matrix,
    sliding_windows = sliding_windows
  )
}


#' @rdname create_trajectory_matrix
#' @export
create_trajectory_matrix.Sap <- function(
    x,
    segment_type = c("motifs", "syllables", "bouts", "segments"),
    data_type = NULL,
    clusters = NULL,
    sample_percent = NULL,
    balanced = FALSE,
    labels = NULL,
    seed = 222,
    window_size = 0.1,
    step_size = 0.005,
    wl = 128,
    ovlp = 50,
    fftw = TRUE,
    flim = c(1, 12),
    cores = NULL,
    verbose= TRUE,
    ...
) {

  if(verbose) message(sprintf("\n=== Starting trajectory analysis for %s using %s%s ===",
                              segment_type[1],
                              data_type[1],
                              ifelse(is.null(clusters), "", sprintf(" (cluster %s)", clusters))))

  # Validate segment_type
  segment_type <- match.arg(segment_type)
  feature_type <- sub("s$", "", segment_type)  # Remove 's' from end

  # Get appropriate data frame based on data_type
  if(is.null(data_type)) {
    segments_df <- x[[segment_type]]
    if(!inherits(segments_df, "segment") || nrow(segments_df) == 0) {
      stop("No segments found in the specified segment type")
    }
  } else if(data_type == "feat.embeds") {
    # Get feature data
    feature_data <- x$features[[feature_type]][[data_type]]

    if (is.null(feature_data)) {
      stop(sprintf("No %s data found in %s features", data_type, feature_type))
    }

    segments_df <- feature_data
  } else {
    stop("Invalid data_type. Must be NULL or 'feat.embeds'")
  }

  # Validate clusters
  if(!is.null(clusters)) {
    if(!"cluster" %in% colnames(segments_df)) {
      stop("No cluster column found in the data")
    }

    # Check which clusters are missing
    available_clusters <- unique(segments_df$cluster)
    missing_clusters <- setdiff(clusters, available_clusters)
    if(length(missing_clusters) > 0) {
      stop(sprintf("The following clusters were not found: %s",
                   paste(missing_clusters, collapse = ", ")))
    }
  }

  # Apply segment selection
  segments_df <- select_segments(segments_df,
                                 clusters = clusters,
                                 balanced = balanced,
                                 sample_percent = sample_percent,
                                 labels = labels,
                                 seed = seed)

  # Call default method to generate spectrograms
  if (fftw) ensure_pkgs("fftw")

  result <- create_trajectory_matrix.default(
    x = segments_df,
    wav_dir = x$base_path,
    window_size = window_size,
    step_size = step_size,
    wl = wl,
    ovlp = ovlp,
    fftw = fftw,
    flim = flim,
    cores = cores,
    ...
  )

  # Save results in the Sap object using feature_type
  x$features[[feature_type]][["traj_mat"]] <- result$spectrogram_matrix
  x$features[[feature_type]][["traj.embeds"]] <- result$sliding_windows

  # Add parameters as attributes
  attrs <- list(
    window_size = window_size,
    step_size = step_size,
    wl = wl,
    ovlp = ovlp,
    flim = flim,
    data_type = data_type,
    clusters = clusters,
    balanced = balanced,
    sample_percent = sample_percent,
    labels = labels,
    seed = seed
  )
  attributes(x$features[[feature_type]][["traj.embeds"]]) <- c(
    attributes(x$features[[feature_type]][["traj.embeds"]]),
    attrs
  )

  invisible(x)
}

#' Create Sliding Windows for Segment
#'
#' @description
#' Internal function to generate overlapping time windows for a segment.
#'
#' @param i Rendition number
#' @param x Data frame with segment information
#' @param window_size Window size in seconds
#' @param step_size Step size in seconds
#' @param ... Additional arguments
#'
#' @return
#' Data frame containing sliding window information
#'
#' @keywords internal
create_sliding_window <- function(
    i,
    x,
    window_size = 0.1,
    step_size = 0.005,
    ...
) {
  # Constrain window generation to ensure sufficient bins
  start <- x$start_time[i]
  end <- x$end_time[i]

  # Calculate max windows that fit in the segment
  max_windows <- floor((end - start - window_size) / step_size) + 1
  window_start_time <- start + (0:(max_windows-1)) * step_size

  # Create sliding window data frame
  sliding_windows <- data.frame(
    filename = x$filename[i],
    day_post_hatch = x$day_post_hatch[i],
    label = x$label[i],
    rendition = i,
    selec = seq_along(window_start_time),
    start_time = window_start_time,
    end_time = window_start_time  + window_size,
    .time = window_start_time - start
  )

  return(sliding_windows)
}
