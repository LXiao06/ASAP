# Plot Trajectory Variability Results

A unified plotting function for results produced by
[`trajectory_variability`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_variability.md),
[`trajectory_width_variability`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_width_variability.md),
or
[`trajectory_umap_occupancy`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_umap_occupancy.md).
Dispatches to the appropriate panel layout based on the `type` field
embedded in `result` and draws significance brackets when statistical
tests are present.

## Usage

``` r
plot_trajectory_variability(
  result,
  palette = "Set1",
  max_annotations = 10,
  ...
)
```

## Arguments

- result:

  A list returned by
  [`trajectory_variability()`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_variability.md),
  [`trajectory_width_variability()`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_width_variability.md),
  or
  [`trajectory_umap_occupancy()`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_umap_occupancy.md).
  Must contain a `type` element (one of `"variability"`,
  `"width_variability"`, `"umap_occupancy"`).

- palette:

  RColorBrewer palette name (default: `"Set1"`). When the number of
  labels exceeds the palette's maximum, colours are interpolated
  automatically via `colorRampPalette`.

- max_annotations:

  Maximum number of pairwise significance brackets to draw per panel
  (default: `10`). When more comparisons exist, the most significant
  pairs are retained and a message is issued.

- ...:

  Currently unused; reserved for future extensions.

## Value

The assembled patchwork object, printed as a side-effect and returned
invisibly so the caller can save or further modify it.

## Details

Panel layouts by result type:

- `"variability"`:

  3 panels: Mean Pairwise Distance · Centroid Dispersion · Path Length

- `"width_variability"`:

  3 panels: Total RMS · Orthogonal RMS · Parallel RMS

- `"umap_occupancy"`:

  4 panels: Occupied Fraction · Occupancy Entropy · Peripheral Fraction
  · kNN Dispersion

Each panel displays a violin + box plot coloured by label. When
`result$tests` is not `NULL`, Kruskal-Wallis p-values are shown as
subtitles and pairwise Wilcoxon p-values appear as significance brackets
above the data.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- trajectory_variability(sap)
plot_trajectory_variability(result)

result2 <- trajectory_width_variability(sap)
plot_trajectory_variability(result2, palette = "Dark2")

result3 <- trajectory_umap_occupancy(sap)
p <- plot_trajectory_variability(result3, max_annotations = 6)
} # }
```
