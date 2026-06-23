# Plot Trajectory Variability Results

A unified plotting function for results produced by
[`trajectory_dispersion`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_dispersion.md),
[`trajectory_path_deviation`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_path_deviation.md),
or
[`trajectory_umap_occupancy`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_umap_occupancy.md).
Dispatches to the appropriate panel layout based on the `type` field
embedded in `result` or retrieved from the SAP object and draws
significance brackets when statistical tests are present.

## Usage

``` r
plot_trajectory_variability(x, ...)

# Default S3 method
plot_trajectory_variability(x, palette = "Set1", max_annotations = 10, ...)

# S3 method for class 'Sap'
plot_trajectory_variability(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  variability_type = c("dispersion", "path_deviation", "umap_occupancy"),
  palette = "Set1",
  max_annotations = 10,
  ...
)
```

## Arguments

- x:

  A list returned by
  [`trajectory_dispersion()`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_dispersion.md),
  [`trajectory_path_deviation()`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_path_deviation.md),
  or
  [`trajectory_umap_occupancy()`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_umap_occupancy.md)
  (for the default method); or a SAP object (for the Sap method).

- ...:

  Additional arguments passed to specific methods.

- palette:

  RColorBrewer palette name (default: `"Set1"`). When the number of
  labels exceeds the palette's maximum, colours are interpolated
  automatically via `colorRampPalette`.

- max_annotations:

  Maximum number of pairwise significance brackets to draw per panel
  (default: `10`). When more comparisons exist, the most significant
  pairs are retained and a message is issued.

- segment_type:

  For SAP objects: Type of segments to visualize ('motifs', 'syllables',
  'bouts', 'segments')

- variability_type:

  For SAP objects: Which computed variability type to plot
  ('dispersion', 'path_deviation', 'umap_occupancy')

## Value

The assembled patchwork object, printed as a side-effect and returned
invisibly so the caller can save or further modify it.

## Details

Panel layouts by result type:

- `"dispersion"`:

  3 panels: Mean Pairwise Distance · Centroid Dispersion · Path Length

- `"path_deviation"`:

  3 panels: Total RMS · Orthogonal RMS · Parallel RMS

- `"umap_occupancy"`:

  4 panels: Occupied Fraction · Occupancy Entropy · Peripheral Fraction
  · kNN Dispersion

Each panel displays a violin + box plot coloured by label. When
statistical tests are not `NULL`, Kruskal-Wallis p-values are shown as
subtitles and pairwise Wilcoxon p-values appear as significance brackets
above the data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plotting directly from a result list
result <- trajectory_dispersion(sap$features$motif$traj.embeds, dims = c("PC1", "PC2"))
plot_trajectory_variability(result)

# Plotting from a SAP object with pre-computed results
sap <- trajectory_dispersion(sap)
plot_trajectory_variability(sap, segment_type = "motifs", variability_type = "dispersion")
} # }
```
