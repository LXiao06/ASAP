# Visualize Motif Boundaries in Audio Data

Generates a heatmap visualization of motif boundaries with optional
cluster filtering and UMAP-based ordering. Incorporates amplitude
envelopes and boundary markers for precise temporal analysis.

## Usage

``` r
plot_motif_boundaries(
  x,
  sample_percent = NULL,
  balanced = FALSE,
  labels = NULL,
  clusters = NULL,
  ordered = FALSE,
  descending = TRUE,
  cores = NULL,
  seed = 222,
  msmooth = c(256, 50),
  color_palette = NULL,
  n_colors = 500,
  contrast = 3,
  marginal_window = 0.1,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A SAP object containing motifs and metadata

- sample_percent:

  Percentage of data to sample from each group (0-100)

- balanced:

  Logical indicating whether to balance group sizes

- labels:

  Character vector of specific labels to include

- clusters:

  Numeric vector of cluster IDs to filter

- ordered:

  Logical indicating UMAP-based ordering

- descending:

  Logical for descending order in UMAP sorting

- cores:

  Number of processing cores (NULL for auto-detection)

- seed:

  Random seed for reproducibility

- msmooth:

  Vector of smoothing parameters for amplitude envelopes

- color_palette:

  Color palette function for heatmap

- n_colors:

  Number of color gradations in palette

- contrast:

  Contrast adjustment for color scaling

- marginal_window:

  Time window extension for motif margins (seconds)

- verbose:

  Logical flag for progress messages

- ...:

  Additional parameters passed to lattice::levelplot

## Value

A lattice plot object showing motif boundaries heatmap

## Details

Key features:

- Integrates motif boundary data from refine_motif_boundaries()

- Supports cluster filtering using feature embeddings

- Enables UMAP-based ordering of motifs

- Visualizes boundaries with colored markers (green = onset, cyan =
  offset)

- Includes automatic duration calculation and marginal window extension

Requires prior execution of refine_motif_boundaries() for boundary
detection. Uses parallel processing for efficient amplitude envelope
calculation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with default parameters
plot_motif_boundaries(sap_object)

# Cluster-specific plot with custom colors
plot_motif_boundaries(sap_object,
                      clusters = c(1,3),
                      ordered = TRUE)
} # }
```
