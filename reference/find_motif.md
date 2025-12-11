# Find Motifs in Song Data

Identifies and extracts motifs from song recordings based on detection
times.

## Usage

``` r
find_motif(x, ...)

# Default S3 method
find_motif(
  x,
  pre_time = NULL,
  lag_time = NULL,
  wav_dir = NULL,
  add_path_attr = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
find_motif(
  x,
  template_name,
  pre_time = NULL,
  lag_time = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to process, either a data frame or SAP object

- ...:

  Additional arguments passed to specific methods

- pre_time:

  Time in seconds before detection point

- lag_time:

  Time in seconds after detection point

- wav_dir:

  For default method: Directory containing WAV files

- add_path_attr:

  For default method: Add wav_dir as attribute (default: TRUE)

- verbose:

  Whether to print processing information (default: TRUE)

- template_name:

  For SAP objects: Name of template to process

## Value

For default method: Data frame containing:

- filename: Source WAV file name

- detection_time: Original detection time

- start_time, end_time: Motif boundaries

- duration: Motif duration

For SAP objects: Updated object with motifs stored in motifs slot

## Details

For detection data frames:

- Requires columns: filename, time

- Validates motif boundaries against audio duration

- Processes each unique audio file

- Returns combined results with metadata

For SAP objects:

- Processes template-based detections

- Organizes results by recording day

- Validates motif boundaries

- Updates object with extracted motifs

## See also

[`detect_template`](https://lxiao06.github.io/ASAP/reference/detect_template.md)
for template detection

## Examples

``` r
if (FALSE) { # \dontrun{
# Find motifs from detection data frame
motifs <- find_motif(detections,
                     pre_time = 0.1,
                     lag_time = 0.2,
                     wav_dir = "path/to/wavs")

# Process with path attribute
motifs <- find_motif(detections,
                     pre_time = 0.1,
                     lag_time = 0.2,
                     wav_dir = "path/to/wavs",
                     add_path_attr = TRUE)

# Find motifs in SAP object
sap_obj <- find_motif(sap_object,
                      template_name = "template1",
                      pre_time = 0.7,
                      lag_time = 0.5)

# Process with custom timing
sap_obj <- find_motif(sap_object,
                      template_name = "template2",
                      pre_time = 0.5,
                      lag_time = 0.3,
                      verbose = TRUE)
} # }
```
