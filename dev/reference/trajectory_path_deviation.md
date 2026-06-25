# Trajectory Path Deviation Analysis

Quantifies rendition-to-rendition trajectory deviation by measuring
residual spread around each label's mean trajectory and decomposing that
spread into orthogonal and parallel components relative to the local
trajectory tangent.

## Usage

``` r
trajectory_path_deviation(x, ...)

# Default S3 method
trajectory_path_deviation(
  x,
  dims = c("PC1", "PC2"),
  trim_fraction = 0.1,
  min_coverage = 0.5,
  time_digits = 6,
  labels = NULL,
  stats = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
trajectory_path_deviation(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  dims = c("PC1", "PC2"),
  trim_fraction = 0.1,
  min_coverage = 0.5,
  time_digits = 6,
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

  Character vector of dimension columns to use (default: c("PC1",
  "PC2"))

- trim_fraction:

  Trim fraction for robust mean trajectory estimation (default: 0.1)

- min_coverage:

  Minimum fraction of renditions that must cover a time step for it to
  contribute to the reference trajectory (default: 0.5)

- time_digits:

  Number of decimal places used to bin `.time` before grouping and
  matching trajectories (default: `6`).

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

For default method: A list (returned invisibly) with the following
elements:

- `type`: Character string `"path_deviation"`, used by
  [`plot_trajectory_variability()`](https://lxiao06.github.io/ASAP/dev/reference/plot_trajectory_variability.md)
  for dispatch

- `width`: Per-rendition width metrics

- `summary`: Summary table with mean and SD for each metric per label

- `mean_trajectories`: Label-specific mean trajectories

- `tangent_vectors`: Label-specific unit tangent vectors

- `tests`: Kruskal-Wallis and pairwise Wilcoxon tests (`NULL` when
  `stats = FALSE` or only one label is present)

For SAP objects: The updated SAP object with results stored in
`x$features[[feature_type]][["trajectory_path_deviation"]]` (returned
invisibly).

Use
[`plot_trajectory_variability`](https://lxiao06.github.io/ASAP/dev/reference/plot_trajectory_variability.md)`(result)`
to visualise the output.

## Details

For each label, the function builds a robust mean trajectory in the
requested dimensions after binning `.time` to `time_digits`, estimates a
local tangent vector at each retained time step, and decomposes each
rendition's residual into:

- Total RMS Residual:

  Overall deviation from the label-specific mean trajectory

- Orthogonal RMS Residual:

  Deviation perpendicular to the local tangent; a direct measure of
  trajectory width / jitter around the backbone

- Parallel RMS Residual:

  Deviation along the local tangent; often reflects phase or advance-lag
  variability

## Examples

``` r
if (FALSE) { # \dontrun{
result <- trajectory_path_deviation(sap, dims = c("PC1", "PC2"))
result$summary
result$width
} # }
```
