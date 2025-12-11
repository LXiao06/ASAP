# Plot Traces of Song Features

Creates line plot visualizations of various acoustic features from audio
segments, including amplitude envelope, fundamental frequency, pitch
goodness, and Wiener entropy. The function supports multiple
visualization options.

## Usage

``` r
plot_traces(x, ...)

# Default S3 method
plot_traces(
  x,
  labels = NULL,
  plot_type = c("combined", "individual", "average", "cv"),
  feature = c("env", "pitch", "goodness", "entropy"),
  alpha = 0.2,
  ncol = 1,
  palette = "Set1",
  ...
)

# S3 method for class 'Sap'
plot_traces(
  x,
  segment_type = c("motifs", "syllables", "segments"),
  feature = c("env", "pitch", "goodness", "entropy"),
  labels = NULL,
  plot_type = c("combined", "individual", "average", "cv"),
  alpha = 0.2,
  ncol = 1,
  palette = "Set1",
  ...
)
```

## Arguments

- x:

  An object to visualize (matrix or SAP object)

- ...:

  Additional arguments passed to specific methods

- labels:

  Optional vector of labels to include (default: NULL, uses all labels)

- plot_type:

  Type of plot: "individual", "average", "cv", or "combined" (default:
  "combined")

- feature:

  Type of feature to plot: "env" (amplitude envelope), "pitch"
  (fundamental frequency), "goodness" (pitch goodness), or "entropy"
  (Wiener entropy)

- alpha:

  Transparency for individual traces (default: 0.2)

- ncol:

  Number of columns for facet wrapping (default: 1)

- palette:

  Color palette for plotting (default: "Set1")

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'segments')

## Value

A ggplot object with the specified trace visualization. For SAP objects,
the function returns the input object invisibly after printing the plot.

## Details

For matrices:

- Requires matrix with 'time_window' attribute

- Column names used as labels

- Supports direct visualization of pre-computed traces

For SAP objects:

- Supports multiple segment types(motifs, syllables, segments)

- Can plot multiple acoustic features (amplitude envelope, fundamental
  frequency, pitch goodness, and Wiener entropy)

- Flexible visualization options

Features:

- env:

  Amplitude envelope - represents the overall loudness/intensity of
  sound over time

- pitch:

  Fundamental frequency (in kHz) - represents the pitch contour

- goodness:

  Pitch goodness - represents the periodicity/quality of pitch
  estimation

- entropy:

  Wiener entropy - measures the width and uniformity of the power
  spectrum

Plot Types:

- individual:

  Shows each rendition's trace, faceted by label

- average:

  Displays mean trace with standard error

- cv:

  Shows coefficient of variation (CV) over time - displays relative
  variability as a percentage at each time point for each label

- combined:

  Shows both individual traces and mean trace

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot amplitude envelope traces from a matrix
plot_traces(sap_obj$features$motif$amp_env,
            feature = "env",
            plot_type = "combined")

# Plot fundamental frequency traces from a SAP object
plot_traces(sap_obj,
            segment_type = "motifs",
            feature = "pitch",
            plot_type = "individual")

# Plot Wiener entropy traces
plot_traces(sap_obj,
            segment_type = "syllables",
            feature = "entropy",
            plot_type = "average")

# Plot coefficient of variation to examine variability
plot_traces(sap_obj,
            segment_type = "motifs",
            feature = "pitch",
            plot_type = "cv")

# Customize visualization
plot_traces(sap_obj$features$motif$amp_env,
            labels = c("BL", "Post"),
            feature = "env",
            plot_type = "average",
            alpha = 0.1,
            ncol = 2)
} # }
```
