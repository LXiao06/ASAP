# Trajectory Variability Analysis

Quantifies trajectory variability across experimental conditions using
three complementary metrics computed in PCA or UMAP embedding space.

## Usage

``` r
trajectory_variability(x, ...)

# Default S3 method
trajectory_variability(
  x,
  dims = c("PC1", "PC2"),
  labels = NULL,
  max_pairs = 5000,
  seed = 222,
  stats = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
trajectory_variability(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  dims = c("PC1", "PC2"),
  labels = NULL,
  max_pairs = 5000,
  seed = 222,
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

  Character vector of dimension columns to use (default: c("PC1",
  "PC2"))

- labels:

  Optional character vector of labels to include

- max_pairs:

  Maximum number of pairwise comparisons per label (default: 5000)

- seed:

  Random seed for reproducible pair sampling (default: 222)

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

- `pairwise`: Data frame of pairwise distances (label, pair_id,
  mean_dist)

- `dispersion`: Data frame of centroid dispersion (label, rendition,
  dispersion)

- `path_length`: Data frame of path lengths (label, rendition,
  path_length)

- `summary`: Summary table with mean and SD for each metric per label

- `tests`: List of statistical test results (`NULL` when
  `stats = FALSE`)

- `type`: Character string `"variability"`, used by
  [`plot_trajectory_variability()`](https://lxiao06.github.io/ASAP/dev/reference/plot_trajectory_variability.md)
  for dispatch

Use
[`plot_trajectory_variability`](https://lxiao06.github.io/ASAP/dev/reference/plot_trajectory_variability.md)`(result)`
to visualise the output.

## Details

Three metrics are computed:

- Mean Pairwise Distance:

  Average Euclidean distance between all pairs of trajectories within
  each condition, measured at matched time points and averaged across
  time. Higher values indicate greater spread among renditions.

- Centroid Dispersion:

  For each rendition, the mean distance to the condition's centroid
  trajectory at each time point. Provides a per-rendition measure of how
  far each trajectory deviates from the "average" path.

- Path Length:

  Total Euclidean distance traveled along each rendition's trajectory
  through the embedding space. Captures the overall complexity and
  extent of each trajectory, independent of the group mean.

Statistical testing is performed automatically:

- Kruskal-Wallis test for overall group differences

- Pairwise Wilcoxon rank-sum tests with Bonferroni correction

## Examples

``` r
if (FALSE) { # \dontrun{
# From SAP object using PC dimensions
result <- trajectory_variability(sap)

# Using UMAP dimensions
result <- trajectory_variability(sap, dims = c("UMAP1", "UMAP2"))

# From trajectory embeddings data frame directly
result <- trajectory_variability(sap$features$motif$traj.embeds,
  dims = c("PC1", "PC2")
)

# Access results
result$summary # summary table
result$tests # statistical tests
result$dispersion # per-rendition dispersion values
} # }
```
