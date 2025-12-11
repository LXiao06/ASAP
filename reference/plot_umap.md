# Plot UMAP Visualization

Creates customizable UMAP visualizations with options for grouping,
highlighting, and customization.

## Usage

``` r
plot_umap(x, ...)

# Default S3 method
plot_umap(
  x,
  dims = c("UMAP1", "UMAP2"),
  subset.by = NULL,
  subset.value = NULL,
  group.by = NULL,
  split.by = NULL,
  cols = NULL,
  pt.size = 0.5,
  stroke = 0.5,
  alpha = 0.3,
  highlight.alpha = NULL,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  highlight.by = NULL,
  highlight.value = NULL,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  background.value = NULL,
  na.value = "grey80",
  ncol = NULL,
  combine = TRUE,
  ...
)

# S3 method for class 'Sap'
plot_umap(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  dims = c("UMAP1", "UMAP2"),
  group.by = NULL,
  split.by = NULL,
  subset.by = NULL,
  subset.value = NULL,
  cols = NULL,
  pt.size = 0.5,
  stroke = 0.5,
  alpha = 0.3,
  highlight.alpha = NULL,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  highlight.by = NULL,
  highlight.value = NULL,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  background.value = NULL,
  na.value = "grey80",
  ncol = NULL,
  combine = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to visualize, either a data frame with UMAP coordinates or a
  SAP object

- ...:

  Additional arguments passed to specific methods

- dims:

  UMAP dimensions to plot (default: c("UMAP1", "UMAP2"))

- subset.by:

  Column name for subsetting data

- subset.value:

  Values to subset by

- group.by:

  Column name for grouping points

- split.by:

  Column name for faceting plots

- cols:

  Custom colors for groups

- pt.size:

  Point size (default: 0.5)

- stroke:

  Point stroke width (default: 0.5)

- alpha:

  Point transparency (default: 0.3)

- highlight.alpha:

  Transparency for highlighted points

- label:

  Whether to add labels (default: FALSE)

- label.size:

  Size of labels (default: 4)

- repel:

  Whether to use repelling labels (default: FALSE)

- highlight.by:

  Column name for highlighting

- highlight.value:

  Values to highlight

- cols.highlight:

  Colors for highlighted points (default: '#DE2D26')

- sizes.highlight:

  Size for highlighted points (default: 1)

- background.value:

  Background group value

- na.value:

  Color for NA values (default: 'grey80')

- ncol:

  Number of columns in multi-plot layout

- combine:

  Whether to combine multiple plots (default: TRUE)

- segment_type:

  For SAP objects: Type of segments to visualize ('motifs', 'syllables',
  'bouts', 'segments')

- verbose:

  For SAP objects: Whether to print progress messages

## Value

For default method: A ggplot object or list of plots For SAP objects:
Updated SAP object with plot as side effect

## Details

This function creates UMAP visualizations with the following features:

- Flexible grouping and highlighting

- Customizable point appearance

- Optional labels and faceting

- Multiple plot combinations

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic UMAP plot from data frame
plot_umap(umap_df, group.by = "cluster")

# Plot with highlighting
plot_umap(umap_df,
          group.by = "cluster",
          highlight.by = "label",
          highlight.value = "a")

# Plot with faceting
plot_umap(umap_df,
          group.by = "cluster",
          split.by = "day_post_hatch")

# Plot from SAP object
plot_umap(sap_obj,
          segment_type = "motifs",
          group.by = "label")

# SAP object plot with custom grouping and highlighting
plot_umap(sap_obj,
          segment_type = "syllables",
          group.by = "label",
          highlight.by = "cluster",
          highlight.value = c(1, 2))
} # }
```
