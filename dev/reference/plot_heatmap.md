# Plot Heatmap of Amplitude Envelopes

Creates heatmap visualizations of amplitude envelopes from audio
segments, supporting multiple data sources and visualization options.

## Usage

``` r
plot_heatmap(x, ...)

# Default S3 method
plot_heatmap(
  x,
  wav_dir = NULL,
  msmooth = c(256, 50),
  color_palette = NULL,
  n_colors = 500,
  contrast = 3,
  ...
)

# S3 method for class 'matrix'
plot_heatmap(
  x,
  labels = NULL,
  color_palette = NULL,
  n_colors = 500,
  contrast = 3,
  reference_lines = NULL,
  reference_line_color = "black",
  main = "Amplitude Envelope Heatmap",
  ylabel = "Labels",
  ...
)

# S3 method for class 'Sap'
plot_heatmap(
  x,
  segment_type = c("motifs", "bouts"),
  sample_percent = NULL,
  balanced = FALSE,
  labels = NULL,
  clusters = NULL,
  cores = NULL,
  seed = 222,
  msmooth = c(256, 50),
  color_palette = NULL,
  n_colors = 500,
  contrast = 3,
  ordered = FALSE,
  descending = TRUE,
  padding_quantile = 0.9,
  window = NULL,
  reference_lines = NULL,
  reference_line_color = "white",
  ylabel = "Labels",
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to visualize (data frame, SAP object, or matrix)

- ...:

  Additional arguments passed to specific methods

- wav_dir:

  For default method: Path to WAV files directory

- msmooth:

  Smoothing parameters c(window_length, overlap_percentage)

- color_palette:

  Function generating color palette

- n_colors:

  Number of colors in heatmap (default: 500)

- contrast:

  Contrast factor for visualization (default: 3)

- labels:

  Optional vector of labels to include

- reference_lines:

  Numeric vector specifying time positions (in seconds) relative to the
  alignment point where vertical dashed reference lines should be drawn.
  Can include both positive and negative values.

- reference_line_color:

  Character vector specifying the color(s) of the reference lines. Can
  be a single color (applied to all lines) or a vector of colors (one
  for each line). If fewer colors than lines are provided, colors will
  be recycled.

- main:

  For matrix method: Plot title

- ylabel:

  Character string specifying the label for the y-axis. (default:
  "Labels")

- segment_type:

  For SAP objects: Type of segments ('motifs', 'bouts', 'syllables',
  'segments')

- sample_percent:

  For SAP objects: Percentage to sample

- balanced:

  For SAP objects: Balance across labels

- clusters:

  Numeric vector of cluster IDs to filter

- cores:

  For SAP objects: Number of processing cores

- seed:

  For SAP objects: Random seed (default: 222)

- ordered:

  For SAP objects: Order by embeddings

- descending:

  For SAP objects: Direction of ordering

- padding_quantile:

  For SAP objects: Quantile for bout padding (default: 0.9)

- window:

  For SAP objects: Numeric vector of length 2 (pre, post) specifying
  time windows in seconds around the alignment point for bouts. When
  NULL (default), windows are auto-calculated. Ignored for motifs.
  Example: `c(0.5, 1.5)`.

- verbose:

  For SAP objects: Print progress messages

## Value

For default method: List containing segments, matrix, and plot For SAP
objects: Updated object with amplitude features For matrices: Lattice
plot object

## Details

For data frames:

- Requires columns: filename, start_time, end_time

- Calculates envelopes for each segment

- Creates matrix of aligned envelopes

For SAP objects:

- Supports multiple segment types

- Optional balanced sampling

- Parallel processing support

- Ordering by feature embeddings

For matrices:

- Direct visualization of pre-computed envelopes

- Label-based organization

- Visual separation between groups

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with data frame
plot_heatmap(segments, wav_dir = "path/to/wavs")

# SAP object with options
plot_heatmap(sap_obj,
             segment_type = "motifs",
             balanced = TRUE,
             ordered = TRUE)

# Matrix with specific labels
plot_heatmap(amp_matrix,
             labels = c("a", "b"),
             contrast = 2)

# Advanced SAP object usage
plot_heatmap(sap_obj,
             segment_type = "bouts",
             sample_percent = 80,
             cores = 4,
             ordered = TRUE,
             descending = FALSE)
} # }
```
