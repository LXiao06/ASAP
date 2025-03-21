
# SAP Metadata --------------------------------------------------------------
# Update date: Feb.4, 2025

#' Create metadata for audio files recorded by SAP2011 (Sound Analysis Pro)
#'
#' @description
#' Creates a metadata data frame from WAV files in specified directories.
#' Parses filenames to extract information about bird ID, recording date, and time.
#'
#' @param base_path Character string specifying the base directory path
#' @param subfolders_to_include Character vector of subfolder names to include.
#'        If NULL, includes all subfolders except those in subfolders_to_exclude
#' @param subfolders_to_exclude Character vector of subfolder names to exclude.
#'        Default excludes 'templates' and 'temp_plots'
#' @param labels Character vector of labels corresponding to each subfolder.
#'        Must match the length of subfolders if provided
#'
#' @return A data frame containing metadata with columns:
#'   \item{filename}{Name of the WAV file}
#'   \item{bird_id}{Extracted bird identifier}
#'   \item{day_post_hatch}{Day post hatch (from subfolder name)}
#'   \item{recording_date}{Date of recording (MM-DD format)}
#'   \item{recording_time}{Time of recording (HH:MM:SS format)}
#'   \item{label}{Label assigned to the subfolder}
#'
#' @examples
#' \dontrun{
#' metadata <- create_sap_metadata(
#'   base_path = "path/to/recordings",
#'   subfolders_to_include = c("day1", "day2"),
#'   labels = c("pre", "post")
#' )
#' }
#'
#' @export
create_sap_metadata <- function(base_path,
                                subfolders_to_include = NULL,
                                subfolders_to_exclude = c("templates", "temp_plots"),
                                labels = NULL) {
  # Get all subfolders if not specified
  if (is.null(subfolders_to_include)) {
    subfolders <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)
    # Exclude specified subfolders
    subfolders <- subfolders[!subfolders %in% subfolders_to_exclude]
  } else {
    subfolders <- subfolders_to_include
  }

  # Validate input: check if labels match subfolder count
  if (!is.null(labels) && length(subfolders) != length(labels)) {
    stop("The number of subfolders (", length(subfolders),
         ") and labels (", length(labels), ") must be the same.")
  }

  # Assign NA labels if not provided
  if (is.null(labels)) {
    labels <- rep(NA, length(subfolders))
  }

  # Construct subfolder paths
  subfolder_paths <- file.path(base_path, subfolders)

  # Initialize metadata list
  metadata_list <- list()

  # Loop through subfolders and their labels
  for (i in seq_along(subfolders)) {
    # Get all .wav files from the current subfolder
    FileList <- list.files(
      path = subfolder_paths[i],
      pattern = "\\.wav$",
      full.names = TRUE,
      recursive = FALSE
    )

    # Create metadata for the current subfolder
    subfolder_metadata <- do.call(rbind, lapply(FileList, function(x) {
      parsed <- parse_filename(x)
      data.frame(
        filename = basename(x),
        bird_id = parsed$bird_id,
        day_post_hatch = basename(dirname(x)),
        recording_date = parsed$recording_date,
        recording_time = format(parsed$recording_time, "%H:%M:%S"),
        label = labels[i],
        stringsAsFactors = FALSE
      )
    }))

    # Append the subfolder's metadata to the metadata list
    metadata_list[[i]] <- subfolder_metadata
  }

  # Combine all metadata into a single dataframe
  metadata_df <- do.call(rbind, metadata_list)

  return(metadata_df)
}

#' Parse filename for SAP metadata
#'
#' @description
#' Internal function to parse WAV filenames and extract metadata components.
#'
#' @param filename Character string of the filename to parse
#' @return A list containing parsed components:
#'   \item{bird_id}{Extracted bird identifier}
#'   \item{recording_date}{Parsed recording date}
#'   \item{recording_time}{Parsed recording time}
#'
#' @keywords internal
parse_filename <- function(filename) {
  # Updated regular expression to match the actual filename pattern
  pattern <- "([^/]+)_(\\d+\\.\\d+)_(.+)\\.wav$"
  match <- regmatches(filename, regexec(pattern, filename))

  if (length(match[[1]]) == 0) {
    stop(paste("Filename format is not correct:", filename))
  }

  # Extract information
  bird_id <- match[[1]][2]
  date_time_str <- match[[1]][4]

  # Convert date and time
  recording_date <- format(as.Date(strptime(date_time_str, "%m_%d_%H_%M_%S")),"%m-%d")
  recording_time <- strptime(date_time_str, "%m_%d_%H_%M_%S")

  return(list(
    bird_id = bird_id,
    recording_date = recording_date,
    recording_time = recording_time
  ))
}


# SAP Object --------------------------------------------------------------
# Update date: Feb.4, 2025

#' Create a Sound Analysis Pro (SAP) Object from Audio Recordings
#'
#' @description
#' Creates a comprehensive SAP object from WAV files in specified directories,
#' with robust input validation and metadata extraction.
#'
#' @param base_path Character string specifying the base directory path containing audio recordings
#' @param subfolders_to_include Character vector of subfolder names to include.
#'        If NULL, includes all subfolders except those in subfolders_to_exclude
#' @param subfolders_to_exclude Character vector of subfolder names to exclude.
#'        Default excludes 'templates' and 'plots'
#' @param labels Character vector of labels corresponding to each subfolder.
#'        Must match the length of subfolders if provided
#'
#' @details
#' This function performs several key operations:
#' \itemize{
#'   \item Validates the input base path and its contents
#'   \item Creates metadata for WAV files using \code{create_sap_metadata()}
#'   \item Constructs a SAP object with metadata and additional tracking information
#'   \item Validates the created SAP object
#' }
#'
#' @return A SAP object containing:
#'   \item{metadata}{Data frame with file and recording metadata}
#'   \item{base_path}{Original base directory path}
#'   \item{misc}{List with creation details and timestamps}
#'
#' @examples
#' \dontrun{
#' # Create SAP object from all recordings
#' sap_obj <- create_sap_object("path/to/recordings")
#'
#' # Create SAP object with specific subfolders and labels
#' sap_obj <- create_sap_object(
#'   base_path = "path/to/recordings",
#'   subfolders_to_include = c("day1", "day2"),
#'   labels = c("pre", "post")
#' )
#' }
#'
#' @seealso
#' \code{\link{create_sap_metadata}} for metadata extraction
#' \code{\link{validate_sap}} for SAP object validation
#'
#' @export
create_sap_object <- function(base_path,
                              subfolders_to_include = NULL,
                              subfolders_to_exclude = c("templates", "plots"),
                              labels = NULL) {
  # Input validation
  if (!is.character(base_path) || length(base_path) != 1) {
    stop("base_path must be a single character string")
  }

  if (!is.null(subfolders_to_include) &&
      (!is.character(subfolders_to_include) || length(subfolders_to_include) == 0)) {
    stop("subfolders_to_include must be NULL or a character vector")
  }

  if (!is.null(subfolders_to_exclude) &&
      (!is.character(subfolders_to_exclude) || length(subfolders_to_exclude) == 0)) {
    stop("subfolders_to_exclude must be NULL or a character vector")
  }

  if (!is.null(labels) &&
      (!is.character(labels) ||
       !is.null(subfolders_to_include) && length(labels) != length(subfolders_to_include))) {
    stop("labels must be NULL or a character vector matching length of subfolders_to_include")
  }

  # Validate and normalize base_path
  tryCatch({
    # Check if path exists
    if (!dir.exists(base_path)) {
      stop("Directory does not exist: ", base_path)
    }

    # Normalize path
    base_path <- normalizePath(base_path, mustWork = TRUE)

    # Check if it's actually a directory
    if (!file.info(base_path)$isdir) {
      stop("Path exists but is not a directory: ", base_path)
    }

    # Check if directory is readable
    if (file.access(base_path, mode = 4) != 0) {
      stop("Directory exists but is not readable: ", base_path)
    }

    # Check if there are any files/subdirectories
    contents <- list.files(base_path, recursive = FALSE)
    if (length(contents) == 0) {
      stop("Directory is empty: ", base_path)
    }

  }, error = function(e) {
    stop("Invalid base_path: ", e$message)
  })

  # Create metadata
  tryCatch({
    metadata <- create_sap_metadata(
      base_path = base_path,
      subfolders_to_include = subfolders_to_include,
      subfolders_to_exclude = subfolders_to_exclude,
      labels = labels
    )

    if (nrow(metadata) == 0) {
      stop("No files found matching the specified criteria")
    }
  }, error = function(e) {
    stop("Error creating metadata: ", e$message)
  })

  # Create initial SAP object using internal constructor
  sap <- tryCatch({
    new_sap(
      metadata = metadata,
      base_path = base_path,
      misc = list(
        creation_date = Sys.time(),
        last_modified = Sys.time(),
        creation_args = list(
          subfolders_included = subfolders_to_include,
          subfolders_excluded = subfolders_to_exclude,
          labels = labels
        )
      )
    )
  }, error = function(e) {
    stop("Error creating SAP object: ", e$message)
  })

  # Validate SAP object
  tryCatch({
    validate_sap(sap)
  }, error = function(e) {
    stop("SAP object validation failed: ", e$message)
  })

  # Report creation success
  message(sprintf("SAP object created successfully with %d files", nrow(metadata)))
  if (!is.null(subfolders_to_include)) {
    message("Included subfolders: ", paste(subfolders_to_include, collapse = ", "))
  }

  return(sap)
}

#' Internal Constructor for Sound Analysis Pro (SAP) Object
#'
#' @description
#' An internal constructor function for creating a new SAP object with
#' predefined default components and structure.
#'
#' @param metadata A data frame containing file and recording metadata
#'        (default is an empty data frame)
#' @param base_path A character string representing the base directory path
#'        (default is an empty character vector)
#' @param templates A template_collection object for storing analysis templates
#'        (default is a new empty template collection)
#' @param motifs A segment object for motif-level analysis
#'        (default is a new empty segment object)
#' @param bouts A segment object for bout-level analysis
#'        (default is a new empty segment object)
#' @param syllables A segment object for syllable-level analysis
#'        (default is a new empty segment object)
#' @param segments A segment object for general segment analysis
#'        (default is a new empty segment object)
#' @param features A list containing feature categories for different analysis levels
#'        (default is an empty list with predefined categories)
#' @param misc A list for storing miscellaneous information
#'        (default is an empty list)
#' @param version A character string representing the SAP object version
#'        (default is "1.0.1")
#'
#' @details
#' This is an internal constructor function used to create a standardized SAP object.
#' It ensures that all SAP objects have a consistent structure with predefined
#' components and default empty values.
#'
#' The function creates a structure with the following components:
#' \itemize{
#'   \item metadata: Information about recorded files
#'   \item base_path: Directory path for the recordings
#'   \item motifs: Motif-level segment information
#'   \item bouts: Bout-level segment information
#'   \item syllables: Syllable-level segment information
#'   \item segments: General segment information
#'   \item templates: Collection of analysis templates
#'   \item features: Feature lists for different analysis levels
#'   \item misc: Miscellaneous metadata and information
#'   \item version: Package version
#' }
#'
#' @return
#' A SAP object of class "Sap" with the specified components
#'
#' @note
#' This is an internal function not intended to be called directly by users.
#' It is used within the package for creating SAP objects.
#'
#' @keywords internal
new_sap <- function(metadata = data.frame(),
                    base_path = character(),
                    templates = new_template_collection(),
                    motifs = new_segment(),
                    bouts = new_segment(),
                    syllables = new_segment(),
                    segments = new_segment(),
                    features = list(
                      motif = list(),
                      bout = list(),
                      syllable = list(),
                      segment = list()
                    ),
                    misc = list(),
                    version = "1.0.1") {

  # Create structure
  x <- structure(
    list(
      metadata = metadata,
      base_path = base_path,
      motifs = motifs,
      bouts = bouts,
      syllables = syllables,
      segments = segments,
      templates = templates,
      features = features,
      misc = misc,
      version = version
    ),
    class = "Sap"
  )

  return(x)
}


#' Validate a Sound Analysis Pro (SAP) Object
#'
#' @description
#' Internal function to perform comprehensive validation of a SAP object,
#' ensuring it meets the required structure and contains all necessary components.
#'
#' @param x An object to be validated as a SAP object
#'
#' @details
#' This is an internal validation function used primarily during SAP object creation
#' and manipulation. It performs rigorous checks on the object's structure:
#' \itemize{
#'   \item Verifies the object is of class 'Sap'
#'   \item Checks for the presence of all required components
#'   \item Validates the type and structure of each component
#'   \item Ensures specific objects have the correct class
#' }
#'
#' Required components include:
#' \itemize{
#'   \item metadata: A data frame with file and recording information
#'   \item base_path: A character string of the base directory
#'   \item motifs: A segment object for motif-level analysis
#'   \item bouts: A segment object for bout-level analysis
#'   \item syllables: A segment object for syllable-level analysis
#'   \item segments: A segment object for segment-level analysis
#'   \item templates: A template_collection object
#'   \item features: A list containing feature categories
#'   \item misc: A list of miscellaneous information
#'   \item version: A character string representing the SAP object version
#' }
#'
#' @return
#' Returns \code{TRUE} if the object passes all validation checks.
#' Throws an informative error if any validation fails.
#'
#' @note
#' This is an internal function not intended to be called directly by users.
#' It is used automatically during SAP object creation and manipulation.
#'
#' @keywords internal
validate_sap <- function(x) {
  # Check if it's a SAP object
  if (!inherits(x, "Sap")) {
    stop("Object must be of class 'Sap'")
  }

  # Check presence of all components
  required_components <- c("metadata", "base_path", "motifs", "bouts",
                           "syllables", "segments", "templates",
                           "features", "misc", "version")
  missing_components <- required_components[!required_components %in% names(x)]

  if (length(missing_components) > 0) {
    stop("Missing required components in SAP object: ",
         paste(missing_components, collapse = ", "))
  }

  # Check that metadata is a data frame
  if (!is.data.frame(x$metadata)) {
    stop("metadata must be a data frame")
  }

  # Check that base_path is a character string
  if (!is.character(x$base_path) || length(x$base_path) != 1) {
    stop("base_path must be a single character string")
  }

  # Check specialized components
  if (!inherits(x$templates, "template_collection")) {
    stop("templates must be a template_collection object")
  }

  # Check segment objects
  segment_components <- c("motifs", "bouts", "syllables", "segments")
  for (comp in segment_components) {
    if (!inherits(x[[comp]], "segment")) {
      stop(sprintf("%s must be a segment object", comp))
    }
  }

  # Check features structure
  required_features <- c("motif", "bout", "syllable", "segment")
  if (!all(required_features %in% names(x$features))) {
    stop("features must contain all required categories: ",
         paste(required_features, collapse = ", "))
  }

  # Check misc is a list
  if (!is.list(x$misc)) {
    stop("misc must be a list")
  }

  # Check version is character
  if (!is.character(x$version) || length(x$version) != 1) {
    stop("version must be a single character string")
  }

  return(TRUE)
}

# Template Object -----------------------------------------------------------
# Update date: Feb.4, 2025

#' Internal Constructor for Template Collection
#'
#' @description
#' An internal constructor function for creating a new template collection
#' with a standardized structure for storing and managing sound analysis templates.
#'
#' @details
#' Creates a template collection object with three key components:
#' \itemize{
#'   \item template_info: A data frame containing metadata about templates
#'   \item template_list: A list to store corTemplate objects
#'   \item template_matches: A list to store template detection results
#' }
#'
#' The template_info data frame includes the following columns:
#' \itemize{
#'   \item template_name: Name of the template
#'   \item start_time: Start time of the template
#'   \item end_time: End time of the template
#'   \item duration: Duration of the template
#'   \item freq_min: Minimum frequency
#'   \item freq_max: Maximum frequency
#'   \item threshold: Detection threshold
#'   \item clip_name: Name of the audio clip
#'   \item clip_path: Path to the audio clip
#'   \item source_file: Original source file name
#'   \item source_file_path: Path to the source file
#'   \item creation_date: Date and time of template creation
#' }
#'
#' @return
#' A template_collection object with an empty but structured template collection
#'
#' @note
#' This is an internal function not intended to be called directly by users.
#' It is used within the package for creating template collections.
#'
#' @keywords internal
new_template_collection <- function() {
  # Define required structure for template_info
  template_info <- data.frame(
    template_name = character(),
    start_time = numeric(),
    end_time = numeric(),
    duration = numeric(),
    freq_min = numeric(),
    freq_max = numeric(),
    threshold = numeric(),
    clip_name = character(),
    clip_path = character(),
    source_file = character(),
    source_file_path = character(),
    creation_date = as.POSIXct(character()),
    stringsAsFactors = FALSE
  )

  # Create structure
  x <- structure(
    list(
      template_info = template_info,
      template_list = list(),      # For storing corTemplate objects
      template_matches = list()     # For storing detection results
    ),
    class = "template_collection"
  )

  return(x)
}

#' Validate Template Collection Object
#'
#' @description
#' An internal validation function to ensure the integrity and
#' correct structure of a template collection object.
#'
#' @param x The template collection object to validate
#'
#' @details
#' Performs comprehensive checks on the template collection:
#' \itemize{
#'   \item Verifies the object is of class 'template_collection'
#'   \item Checks for the presence of required components
#'   \item Validates the structure of the template_info data frame
#'   \item Ensures template_list and template_matches are lists
#' }
#'
#' Required columns in template_info include:
#' \itemize{
#'   \item template_name
#'   \item start_time
#'   \item end_time
#'   \item duration
#'   \item freq_min
#'   \item freq_max
#'   \item threshold
#'   \item clip_name
#'   \item clip_path
#'   \item source_file
#'   \item source_file_path
#'   \item creation_date
#' }
#'
#' @return
#' Returns \code{TRUE} if the object passes all validation checks.
#' Throws an informative error if any validation fails.
#'
#' @note
#' This is an internal validation function not intended for direct user calls.
#' It is used automatically during template collection creation and manipulation.
#'
#' @keywords internal
validate_template_collection <- function(x) {
  # Check class
  if (!inherits(x, "template_collection")) {
    stop("Object must be of class 'template_collection'")
  }

  # Check structure
  if (!all(c("template_info", "template_list", "template_matches") %in% names(x))) {
    stop("template_collection must contain template_info, template_list, and template_matches")
  }

  # Check template_info structure
  required_cols <- c("template_name", "start_time", "end_time", "duration",
                     "freq_min", "freq_max", "threshold", "clip_name",
                     "clip_path", "source_file", "source_file_path", "creation_date")

  if (!all(required_cols %in% names(x$template_info))) {
    stop("template_info missing required columns")
  }

  # Check that template_list and template_matches are lists
  if (!is.list(x$template_list) || !is.list(x$template_matches)) {
    stop("template_list and template_matches must be lists")
  }

  return(TRUE)
}


# Segment Object ----------------------------------------------------------
# Update date: Feb.4, 2025

#' Internal Constructor for Segment Object
#'
#' @description
#' An internal constructor function for creating a segment object
#' with a standardized structure for storing acoustic segment information.
#'
#' @param x A data frame containing segment information (default is an empty data frame)
#'
#' @details
#' Creates a segment object with a predefined structure. If no data frame is provided,
#' it initializes an empty segment with the following columns:
#' \itemize{
#'   \item filename: Name of the source audio file
#'   \item day_post_hatch: Numeric day post-hatch
#'   \item label: Categorical label for the segment
#'   \item start_time: Numeric start time of the segment
#'   \item end_time: Numeric end time of the segment
#'   \item duration: Numeric duration of the segment
#' }
#'
#' When a data frame is provided, the function:
#' \itemize{
#'   \item Validates the presence of required columns
#'   \item Converts columns to appropriate data types
#'   \item Calculates segment duration if not provided
#' }
#'
#' @return
#' A segment object of class "segment" with the specified segment information
#'
#' @note
#' This is an internal function not intended to be called directly by users.
#' It is used within the package for creating segment objects.
#'
#' @keywords internal
new_segment <- function(x = data.frame()) {
  # Define required columns
  required_cols <- c("filename", "day_post_hatch", "label",
                     "start_time", "end_time")

  # Initialize empty data frame if x is empty
  if (nrow(x) == 0) {
    x <- data.frame(
      filename = character(),
      day_post_hatch = numeric(),
      label = character(),
      start_time = numeric(),
      end_time = numeric(),
      duration = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    # Basic validation
    if (!is.data.frame(x)) {
      stop("Input must be a data frame")
    }

    if (!all(required_cols %in% names(x))) {
      missing_cols <- setdiff(required_cols, names(x))
      stop("Missing required columns: ",
           paste(missing_cols, collapse = ", "))
    }

    # Ensure correct column types
    x$filename <- as.character(x$filename)
    x$day_post_hatch <- as.numeric(x$day_post_hatch)
    x$label <- as.character(x$label)
    x$start_time <- as.numeric(x$start_time)
    x$end_time <- as.numeric(x$end_time)

    # Add duration if not present
    if (!"duration" %in% names(x)) {
      x$duration <- x$end_time - x$start_time
    }
  }

  # Set class
  class(x) <- c("segment", "data.frame")

  return(x)
}

#' Validate Segment Object
#'
#' @description
#' An internal validation function to ensure the integrity and
#' correctness of a segment object.
#'
#' @param x The segment object to validate
#'
#' @details
#' Performs comprehensive checks on the segment object:
#' \itemize{
#'   \item Verifies the object is of class 'segment'
#'   \item Checks for invalid time representations
#'   \item Ensures segment times are non-negative
#'   \item Validates that end times are not before start times
#' }
#'
#' Specific validation checks include:
#' \itemize{
#'   \item Confirming the object is a segment class
#'   \item Checking that end times are greater than or equal to start times
#'   \item Ensuring all time values are non-negative
#' }
#'
#' @return
#' Returns \code{TRUE} if the object passes all validation checks.
#' Throws an informative error if any validation fails.
#'
#' @note
#' This is an internal validation function not intended for direct user calls.
#' It is used automatically during segment object creation and manipulation.
#'
#' @keywords internal
validate_segment <- function(x) {
  # Check class
  if (!inherits(x, "segment")) {
    stop("Object must be of class 'segment'")
  }

  # Check for negative durations
  if (any(x$end_time < x$start_time)) {
    stop("Found negative durations (end_time < start_time)")
  }

  # Check for negative times
  if (any(x$start_time < 0) || any(x$end_time < 0)) {
    stop("Found negative times")
  }

  return(TRUE)
}

# User-facing helper - handles all data cleaning and organization

#' Convert data to a segment object
#'
#' @param x A data frame to convert to a segment object
#'
#' @return A validated segment object
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   filename = "song.wav",
#'   day_post_hatch = 1,
#'   label = "a",
#'   start_time = 0,
#'   end_time = 1,
#'   duration = 1
#' )
#' segment <- as_segment(df)
#' }
#'
#' @export
as_segment <- function(x) {
  # If already a segment, check if needs processing
  if (inherits(x, "segment")) {
    needs_processing <- FALSE
  } else {
    needs_processing <- TRUE
  }

  # Create segment if needed
  if (needs_processing) {
    x <- new_segment(x)
  }

  # Validate
  validate_segment(x)

  # Check for NA values in required columns
  required_cols <- c("filename", "day_post_hatch", "label",
                     "start_time", "end_time", "duration")
  for (col in required_cols) {
    if (any(is.na(x[[col]]))) {
      stop("Found NA values in column: ", col)
    }
  }

  # Reorder columns
  standard_cols <- c("filename", "day_post_hatch", "label",
                     "start_time", "end_time", "duration")
  extra_cols <- setdiff(names(x), standard_cols)
  if (length(extra_cols) > 0) {
    x <- x[, c(standard_cols, extra_cols)]
  } else {
    x <- x[, standard_cols]
  }

  return(x)
}

#' Print method for segment objects
#'
#' @param x A segment object
#' @param ... Additional arguments passed to print
#'
#' @exportS3Method print segment
print.segment <- function(x, ...) {
  cat("Segment object with", nrow(x), "rows\n")
  cat("Time range:", round(min(x$start_time), 3), "to",
      round(max(x$end_time), 3), "seconds\n")
  cat("Days:", length(unique(x$day_post_hatch)),
      "Labels:", length(unique(x$label)), "\n")
  cat("\nFirst few rows:\n")
  print.data.frame(head(x[, c("filename", "day_post_hatch", "label",
                              "start_time", "end_time", "duration")]))
}

#' Summary method for segment objects
#'
#' @param object A segment object
#' @param ... Additional arguments passed to summary
#'
#' @exportS3Method summary segment
summary.segment <- function(object, ...) {
  cat("Segment Summary:\n")
  cat("-----------------\n")
  cat("Total segments:", nrow(object), "\n")
  cat("Unique files:", length(unique(object$filename)), "\n")
  cat("Days:", paste(unique(object$day_post_hatch), collapse = ", "), "\n")
  cat("Labels:", paste(unique(object$label), collapse = ", "), "\n")
  cat("\nDuration statistics (seconds):\n")
  print(summary(object$duration))
}
