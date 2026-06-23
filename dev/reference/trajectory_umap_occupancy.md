# Trajectory UMAP Occupancy Analysis

Quantifies rendition-to-rendition diversity in local UMAP occupancy
using a shared 2D grid and same-label neighborhood structure. This is
intended as a complementary analysis to trajectory width in PCA space
when the biological question is about local latent-state exploration
rather than calibrated geometric variance.

## Usage

``` r
trajectory_umap_occupancy(x, ...)

# Default S3 method
trajectory_umap_occupancy(
  x,
  dims = c("UMAP1", "UMAP2"),
  grid_n = 40,
  k = 15,
  peripheral_quantile = 0.2,
  labels = NULL,
  stats = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
trajectory_umap_occupancy(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  dims = c("UMAP1", "UMAP2"),
  grid_n = 40,
  k = 15,
  peripheral_quantile = 0.2,
  labels = NULL,
  stats = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to analyze: a trajectory embeddings data frame or SAP object

- ...:

  Additional arguments

- dims:

  Character vector of UMAP columns to use (default: c("UMAP1", "UMAP2"))

- grid_n:

  Number of bins per UMAP axis for occupancy calculations (default: 40)

- k:

  Number of same-label nearest neighbors used for local dispersion
  (default: 15)

- peripheral_quantile:

  Quantile of global occupied-bin density used to define peripheral bins
  (default: 0.2)

- labels:

  Optional character vector of labels to include

- stats:

  Logical. If `TRUE` (default), run Kruskal-Wallis and pairwise Wilcoxon
  tests. Set to `FALSE` to skip statistical testing and return `NULL`
  for `tests`.

- verbose:

  Whether to print progress messages (default: TRUE)

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'bouts',
  'segments')

## Value

A list (returned invisibly) with the following elements:

- `type`: Character string `"umap_occupancy"`, used by
  [`plot_trajectory_variability()`](https://lxiao06.github.io/ASAP/dev/reference/plot_trajectory_variability.md)
  for dispatch

- `occupancy`: Per-rendition occupancy metrics

- `summary`: Summary table with mean and SD for each metric per label

- `annotated_points`: Original data with occupancy annotations

- `bin_counts`: Shared-grid counts per label and bin

- `grid_info`: Grid settings and peripheral threshold metadata

- `tests`: Kruskal-Wallis and pairwise Wilcoxon tests (`NULL` when
  `stats = FALSE` or only one label is present)

Use
[`plot_trajectory_variability`](https://lxiao06.github.io/ASAP/dev/reference/plot_trajectory_variability.md)`(result)`
to visualise the output.

## Details

The function computes four per-rendition metrics in a shared UMAP grid:

- Occupied Fraction:

  Fraction of grid bins visited by the rendition

- Occupancy Entropy:

  Shannon entropy of the rendition's occupancy over grid bins,
  normalized to `[0, 1]` using the full shared grid

- Peripheral Fraction:

  Fraction of the rendition's points falling in globally low-density
  bins of the shared UMAP manifold

- Same-Label kNN Dispersion:

  Average distance from each point to its same-label nearest neighbors
  in UMAP space

## Examples

``` r
if (FALSE) { # \dontrun{
result <- trajectory_umap_occupancy(sap)
result$summary
head(result$occupancy)
} # }
```
