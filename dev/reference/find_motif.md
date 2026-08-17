# Find Motifs in Song Data

Identifies and extracts motifs from song recordings based on detection
times. Supports single or multiple templates with template-specific pre-
and lag-times.

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
  template_name = NULL,
  pre_time = NULL,
  lag_time = NULL,
  day_post_hatch = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to process, either a data frame of detections or a SAP
  object.

- ...:

  Additional arguments passed to specific methods.

- pre_time:

  Time in seconds before detection point. Can be a single numeric value
  (broadcast to all templates) or a numeric vector matching
  `template_name` in length (optionally named).

- lag_time:

  Time in seconds after detection point. Can be a single numeric value
  (broadcast to all templates) or a numeric vector matching
  `template_name` in length (optionally named).

- wav_dir:

  For default method: Directory containing WAV files.

- add_path_attr:

  For default method: Add wav_dir as attribute (default: TRUE).

- verbose:

  Whether to print processing information (default: TRUE).

- template_name:

  For SAP objects: Character vector of template name(s) to process. If
  `NULL` (default), processes all templates found in `template_matches`.

- day_post_hatch:

  For SAP objects: Numeric value(s) to use for `day_post_hatch` when it
  cannot be determined from subfolder names (e.g. non-numeric folder
  names like `"FD_661_667"`). Can be a single numeric value (applied to
  all folders) or a named numeric vector mapping subfolder names to
  day-post-hatch values (e.g.
  `c(FD_661_667 = 61, PD_661_667 = 65, UD = 90)`). If `NULL` (default),
  the value is parsed from the subfolder name.

## Value

For default method: Data frame containing:

- filename: Source WAV file name

- detection_time: Original detection time

- start_time, end_time: Motif boundaries

- duration: Motif duration

For SAP objects: Updated SAP object with motifs stored in the `motifs`
slot as a validated `segment` object.

## Details

For detection data frames:

- Requires columns: filename, time

- Validates motif boundaries against audio duration

- Processes each unique audio file

- Returns combined results with metadata

For SAP objects:

- Processes template-based detections for one or multiple templates

- Applies template-specific `pre_time` and `lag_time` boundaries

- Organizes results by recording day / subfolder

- Validates motif boundaries against recording durations

- Stores all extracted motifs in a unified `x$motifs` segment data frame
  with a `template_name` column identifying the source template

## See also

[`detect_template`](https://lxiao06.github.io/ASAP/dev/reference/detect_template.md)
for template detection

## Examples

``` r
if (FALSE) { # \dontrun{
# Find motifs from detection data frame
motifs <- find_motif(detections,
                     pre_time = 0.1,
                     lag_time = 0.2,
                     wav_dir = "path/to/wavs")

# Find motifs in SAP object with single template
sap_obj <- find_motif(sap_object,
                      template_name = "template1",
                      pre_time = 0.7,
                      lag_time = 0.5)

# Find motifs with multiple templates and template-specific boundaries
sap_obj <- find_motif(sap_object,
                      template_name = c("c", "c1", "c2"),
                      pre_time = c(0.58, 0.78, 1.68),
                      lag_time = c(0.70, 0.87, 1.80))

# Broadcast same pre_time and lag_time across multiple templates
sap_obj <- find_motif(sap_object,
                      template_name = c("c", "c1", "c2"),
                      pre_time = 0.5,
                      lag_time = 0.5)
} # }
```
