# Create a Sound Analysis Pro (SAP) Object from Audio Recordings

Creates a comprehensive SAP object from WAV files in specified
directories, with robust input validation and metadata extraction.

## Usage

``` r
create_sap_object(
  base_path,
  subfolders_to_include = NULL,
  subfolders_to_exclude = c("templates", "plots"),
  labels
)
```

## Arguments

- base_path:

  Character string specifying the base directory path containing audio
  recordings

- subfolders_to_include:

  Character vector of subfolder names to include. If NULL, includes all
  subfolders except those in subfolders_to_exclude

- subfolders_to_exclude:

  Character vector of subfolder names to exclude. Default excludes
  'templates' and 'plots'

- labels:

  Character vector of labels corresponding to each subfolder. Must match
  the length of subfolders

## Value

A SAP object containing:

- metadata:

  Data frame with file and recording metadata

- base_path:

  Original base directory path

- misc:

  List with creation details and timestamps

## Details

This function performs several key operations:

- Validates the input base path and its contents

- Creates metadata for WAV files using
  [`create_sap_metadata()`](https://lxiao06.github.io/ASAP/reference/create_sap_metadata.md)

- Constructs a SAP object with metadata and additional tracking
  information

- Validates the created SAP object

## See also

[`create_sap_metadata`](https://lxiao06.github.io/ASAP/reference/create_sap_metadata.md)
for metadata extraction
[`validate_sap`](https://lxiao06.github.io/ASAP/reference/validate_sap.md)
for SAP object validation

## Examples

``` r
if (FALSE) { # \dontrun{
# Create SAP object from all recordings
sap_obj <- create_sap_object("path/to/recordings")

# Create SAP object with specific subfolders and labels
sap_obj <- create_sap_object(
  base_path = "path/to/recordings",
  subfolders_to_include = c("day1", "day2"),
  labels = c("pre", "post")
)
} # }
```
