# Detect Templates in Song Data

Performs template matching on audio files using correlation-based
detection.

## Usage

``` r
detect_template(x, ...)

# Default S3 method
detect_template(
  x,
  template,
  cor.method = "pearson",
  save_plot = FALSE,
  plot_dir = NULL,
  proximity_window = NULL,
  ...
)

# S3 method for class 'Sap'
detect_template(
  x,
  day = NULL,
  indices = NULL,
  template_name,
  threshold = NULL,
  cores = NULL,
  cor.method = "pearson",
  save_plot = FALSE,
  plot_percent = 10,
  verbose = TRUE,
  proximity_window = NULL,
  ...
)
```

## Arguments

- x:

  An object to process, either a file path or SAP object

- ...:

  Additional arguments passed to specific methods

- template:

  For default method: A template object created by create_template()

- cor.method:

  Correlation method ("pearson" or "spearman")

- save_plot:

  Whether to save detection plots

- plot_dir:

  For default method: Directory to save plots

- proximity_window:

  Time window in seconds to filter nearby detections (NULL to disable
  filtering). Only the detection with the highest score within each
  window is retained.

- day:

  For SAP objects: Numeric vector of days to process

- indices:

  For SAP objects: Numeric vector of indices to process

- template_name:

  For SAP objects: Name of template to use

- threshold:

  For SAP objects: New threshold value

- cores:

  For SAP objects: Number of cores for parallel processing

- plot_percent:

  For SAP objects: Percentage of files to plot (default: 10)

- verbose:

  For SAP objects: Whether to print progress messages

## Value

For default method: A data frame containing detection results with
columns:

- filename: Name of the processed file

- time: Time point of detection

- score: Correlation score

For SAP objects: Updated SAP object with detection results stored in
template_matches

## Details

For WAV files:

- Validates input file and template

- Performs correlation matching

- Finds peaks in correlation scores

- Optionally filters nearby detections (if proximity_window is set)

- Optionally saves detection plots

For SAP objects:

- Parallel processing support

- Day-specific processing

- Optional threshold adjustment

- Progress tracking and reporting

- Selective plot generation

- Filtering of nearby detections when proximity_window is specified

## Proximity Filtering

When `proximity_window` is specified, the function will filter
detections that occur within the specified time window (in seconds). For
each group of detections within the window, only the one with the
highest score is kept. This is useful for removing false positive
detections.

## See also

[`create_template`](https://lxiao06.github.io/ASAP/reference/create_template.md)
for creating templates

## Examples

``` r
if (FALSE) { # \dontrun{
# Detect template in single WAV file
detections <- detect_template("path/to/song.wav",
                             template = template_obj,
                             save_plot = TRUE)

# Detect template in SAP object
sap_obj <- detect_template(sap_object,
                          template_name = "template1",
                          day = c(30, 40),
                          threshold = 0.7,
                          cores = 4)

# Process specific indices with plots
sap_obj <- detect_template(sap_object,
                          template_name = "template1",
                          indices = 1:10,
                          save_plot = TRUE)

# Filter nearby detections within 0.5 seconds
sap_obj <- detect_template(sap_object,
                          template_name = "template1",
                          proximity_window = 0.5)
} # }
```
