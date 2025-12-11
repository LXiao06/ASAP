# Plot UMAP Visualization for Trajectory Analysis

Creates UMAP visualizations optimized for trajectory analysis, with
support for continuous color mapping and overlay comparisons.

## Usage

``` r
plot_umap2(x, ...)

# Default S3 method
plot_umap2(
  x,
  dims = c("UMAP1", "UMAP2"),
  color.by = ".time",
  split.by = "label",
  order.by = "day_post_hatch",
  pt.size = 1.2,
  alpha_range = c(0.1, 0.5),
  ncol = NULL,
  title = NULL,
  overlay_mode = FALSE,
  base_label = NULL,
  compare_labels = NULL,
  base_color = "steelblue",
  compare_color = "orangered",
  ...
)

# S3 method for class 'Sap'
plot_umap2(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  data_type = c("feat.embeds", "traj.embeds"),
  dims = c("UMAP1", "UMAP2"),
  color.by = ".time",
  split.by = "label",
  order.by = "day_post_hatch",
  pt.size = 1.2,
  alpha_range = c(0.1, 0.5),
  ncol = NULL,
  title = NULL,
  overlay_mode = FALSE,
  base_label = NULL,
  compare_labels = NULL,
  base_color = "steelblue",
  compare_color = "orangered",
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to visualize, either a data frame or SAP object

- ...:

  Additional arguments passed to specific methods

- dims:

  UMAP dimensions to plot (default: c("UMAP1", "UMAP2"))

- color.by:

  Column for continuous color mapping (default: ".time")

- split.by:

  Column for faceting (default: "label")

- order.by:

  Column for ordering facets (default: "day_post_hatch")

- pt.size:

  Point size (default: 1.2)

- alpha_range:

  Range for alpha transparency (default: c(0.1, 0.5))

- ncol:

  Number of columns in layout

- title:

  Plot title

- overlay_mode:

  Whether to create overlay comparisons (default: FALSE)

- base_label:

  Base label for overlay comparison

- compare_labels:

  Labels to compare against base

- base_color:

  Color for base label (default: "steelblue")

- compare_color:

  Color for comparison labels (default: "orangered")

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'bouts',
  'segments')

- data_type:

  For SAP objects: Type of embedding data ('feat.embeds', 'traj.embeds')

- verbose:

  For SAP objects: Whether to print progress messages

## Value

For default method: A ggplot object For SAP objects: Updated SAP object
with plot as side effect

## Details

Supports two visualization modes:

- Standard mode: Continuous color mapping for trajectory visualization

- Overlay mode: Direct comparison between trajectory patterns

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic trajectory plot
plot_umap2(traj_df, color.by = ".time")

# Overlay comparison plot
plot_umap2(traj_df,
           overlay_mode = TRUE,
           base_label = "a",
           compare_labels = c("b", "c"))

# Plot motif trajectories from SAP object
plot_umap2(sap_obj,
           segment_type = "motifs",
           data_type = "traj.embeds",
           color.by = ".time")

# Compare trajectories between labels
plot_umap2(sap_obj,
           segment_type = "motifs",
           data_type = "traj.embeds",
           overlay_mode = TRUE,
           base_label = "pre")
} # }
```
