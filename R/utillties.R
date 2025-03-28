# utilities

#' @importFrom tuneR readWave
#' @importFrom rlang `%||%`
#' @importFrom tools file_path_sans_ext
NULL

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
      mean_duration = mean(end_time - start_time)
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
        mean_duration = mean(end_time - start_time)
      ) |>
      dplyr::arrange(dph)
    print(final_summary)
  }

  #  Store day ordering as an attribute
  attr(segments_df, "label_order") <- unique(segments_df$label)

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

#' # internal use
#'
#' #' Helper function for SAP analysis
#' #'
#' #' @description Provides help information for SAP (Sound Analysis Pro) object operations
#' #' including object creation, clip creation, template creation, and motif detection.
#' #'
#' #' @param topic Character string specifying the help topic.
#' #'              Options are: "sap", "clip", "template", "motif"
#' #'
#' #' @return Prints help information to the console
#' #'
#' #' @examples
#' #' Help.sap()  # Shows available topics
#' #' Help.sap("sap")  # Shows help for SAP object creation
#' #'
#' #' @export
#' Help.sap <- function(topic = NULL) {
#'   if (is.null(topic)) {
#'     cat("Available topics: 'sap', 'clip', 'template', 'motif'\n")
#'     cat("Use Help.sap('topic') to get specific information\n")
#'     return()
#'   }
#'
#'   if (topic == "sap") {
#'     cat("SAP Object Creation and Basic Information:\n\n")
#'     cat("Example code:\n")
#'     cat('sap <- create_sap_object(\n  base_path = "/path/to/recordings",\n  subfolders_to_include = c("190", "201", "203"),\n  labels = c("BL", "Post", "Rec2")\n)\n\n')
#'
#'     cat("Useful inspection commands:\n")
#'     cat("1. class(sap)           # Shows the object type\n")
#'     cat("2. head(sap$metadata,5) # Displays first 5 rows of metadata\n")
#'     cat("3. sap$base_path        # Shows the base directory path\n")
#'   }
#'
#'   else if (topic == "clip") {
#'     cat("Audio Clip Creation:\n\n")
#'     cat("Example code:\n")
#'     cat('sap <- create_audio_clip(sap,\n  indices = c(1, 5),\n  start_time = c(1, 2),\n  end_time = c(2.5, 3.5),\n  clip_names = c("m1", "m2"))\n\n')
#'
#'     cat("Useful inspection commands:\n")
#'     cat("1. class(sap$template)                  # Shows template object type\n")
#'     cat("2. sap$template$template_info[,1:8]     # Displays template information\n")
#'   }
#'
#'   else if (topic == "template") {
#'     cat("Template Creation:\n\n")
#'     cat("Example code:\n")
#'     cat('sap <- create_template(sap,\n  template_name = "d",\n  clip_name = "m1",\n  start_time = 0.72,\n  end_time = 0.84,\n  freq_min = 1,\n  freq_max = 10,\n  write_template = TRUE)\n\n')
#'
#'     cat("Useful inspection commands:\n")
#'     cat("1. names(sap$templates$template_list)    # Shows available template names\n")
#'     cat("2. unique(sap$metadata$day_post_hatch)   # Shows available days\n")
#'     cat("3. sap$template$template_info[,1:8]      # Displays template information\n")
#'   }
#'
#'   else if (topic == "motif") {
#'     cat("Motif Detection:\n\n")
#'     cat("Example code:\n")
#'     cat('sap <- find_motifs(sap,\n  template_name = "d",\n  pre_time = 0.7,\n  lag_time = 0.5)\n\n')
#'
#'     cat("Useful inspection commands:\n")
#'     cat("1. head(sap$motifs, 5)    # Shows first 5 detected motifs\n")
#'   }
#'
#'   else {
#'     cat("Topic not found. Available topics: 'sap', 'clip', 'template', 'motif'\n")
#'   }
#' }
#'
#'
#'
#'
