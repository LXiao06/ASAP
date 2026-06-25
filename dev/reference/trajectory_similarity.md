# Trajectory Similarity to Reference Path

Measures trajectory similarity by quantifying how individual trial PC
trajectories compare to a reference path (e.g., adult/crystallized
song). Uses complementary distance metrics (RMS, Fr\\echet, DTW) with
unified variability scaling and correlation analysis.

## Usage

``` r
trajectory_similarity(x, ...)

# Default S3 method
trajectory_similarity(
  x,
  dims = c("PC1", "PC2"),
  reference_label = NULL,
  labels = NULL,
  metrics = c("rms", "correlation"),
  interpolate_n = NULL,
  trim_fraction = 0.1,
  min_coverage = 0.5,
  time_digits = 6,
  similarity_baseline = c("reference", "zero"),
  similarity_scale_multiplier = 1,
  stats = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
trajectory_similarity(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  dims = c("PC1", "PC2"),
  reference_label = NULL,
  labels = NULL,
  metrics = c("rms", "correlation"),
  interpolate_n = NULL,
  trim_fraction = 0.1,
  min_coverage = 0.5,
  time_digits = 6,
  similarity_baseline = c("reference", "zero"),
  similarity_scale_multiplier = 1,
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

- reference_label:

  Character. Label to use as reference trajectory (default: NULL, uses
  last sorted label)

- labels:

  Optional character vector of labels to include

- metrics:

  Character vector of metrics to compute, or `"all"` for all four
  metrics (default: c("rms", "correlation"))

- interpolate_n:

  Optional integer. If NULL (default), use matched time points. If
  integer, resample paths to interpolate_n equally spaced points

- trim_fraction:

  Numeric. Trim fraction for robust reference path building (default:
  0.1, removes 10% from each tail)

- min_coverage:

  Minimum fraction of reference-label renditions that must cover a
  binned time step for it to contribute to the reference trajectory
  (default: 0.5).

- time_digits:

  Number of decimal places used to bin `.time` before grouping and
  matching trajectories (default: `6`).

- similarity_baseline:

  How distance similarities are transformed. `"reference"` (default)
  treats distances within the metric-specific reference-label median as
  converged before applying the exponential transform. `"zero"` uses the
  raw normalized distance.

- similarity_scale_multiplier:

  Multiplier applied to metric-specific reference scales before
  transforming distance similarities (default: `1`).

- stats:

  Logical. If `TRUE` (default), run Kruskal-Wallis and pairwise Wilcoxon
  tests on non-reference labels. Set to `FALSE` to skip statistical
  testing

- verbose:

  Whether to print progress messages (default: TRUE)

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'bouts',
  'segments')

## Value

For default method: A list (returned invisibly) with:

- `type`: Character string `"similarity"`

- `dims`: Requested dimensions

- `reference_label`: Reference label used

- `trim_fraction`: Trim fraction applied

- `min_coverage`: Minimum coverage threshold applied to the reference
  path

- `time_digits`: Decimal places used to bin `.time`

- `reference_scale`: Scalar variability scale

- `reference_scales`: Metric-specific variability scales

- `similarity_baseline`: Similarity transform baseline

- `similarity_scale_multiplier`: Similarity scale multiplier

- `reference_path`: Data frame of reference trajectory

- `metrics`: Metrics computed

- `interpolate_n`: Interpolation info

- `similarity`: Per-rendition results with all metric columns

- `summary`: Per-label summary statistics

- `tests`: Statistical test results (`NULL` if stats=FALSE or only one
  non-reference label)

For SAP objects: Updated object with results stored in
`x$features[[feature_type]][["trajectory_similarity"]]`

## Details

**Similarity Metrics:**

- RMS distance: Pointwise Euclidean distance, emphasizes persistent
  deviations across the full trajectory

- Fr\\echet distance: Curve shape distance; handles variable-length
  paths

- DTW distance: Dynamic time warping; captures shape similarity under
  timing shifts

- Correlation: Per-dimension Pearson r, averaged; intrinsically
  normalized (scale-free, range \\\[-1, 1\]\\)

**Reference Building:** The reference trajectory is the robust
trimmed-mean path of the reference label. Trim fraction (default 0.1)
removes outlier reference renditions at each time point to avoid
distortion. Time values are first rounded to `time_digits`, and only
binned time points covered by at least `min_coverage` of reference
renditions are retained.

**Variability Scaling:** Metric-specific reference scales are computed
as the median reference-label rendition distance to the reference path
for RMS, Fr\\echet, and DTW. These scales are used to normalize distance
metrics (not correlation, which is already scale-free).

**Similarity Scores:** Distance metrics are normalized by their
metric-specific reference scale. By default, the reference scale is
multiplied by `similarity_scale_multiplier = 1` before transformation.
By default, `similarity_baseline = "reference"` transforms excess
distance beyond the scaled reference threshold:
`similarity = exp(-pmax(normalized_distance - 1, 0))`. With
`similarity_baseline = "zero"`, the transform is:
`similarity = exp(-normalized_distance)`.

Interpretation with the default reference baseline: similarity
approximately equals 1 means within the scaled reference threshold,
similarity approximately equals 0.37 means one additional scaled
reference unit beyond that threshold, and similarity less than 0.1 means
far from the reference path.

**Statistical Testing:** Reference label is included in results and
plots but excluded from Kruskal-Wallis and pairwise Wilcoxon tests (it
defines the target, not a similarity target).

## See also

[`plot_trajectory_similarity`](https://lxiao06.github.io/ASAP/dev/reference/plot_trajectory_similarity.md)
for visualization,
[`trajectory_path_deviation`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_path_deviation.md)
for width-based trajectory analysis

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute similarity from trajectory embeddings
result <- trajectory_similarity(sap$features$motif$traj.embeds,
  dims = c("PC1", "PC2")
)

# From SAP object
sap <- trajectory_similarity(sap)

# Custom reference and metrics
result <- trajectory_similarity(sap,
  reference_label = "Adult",
  metrics = c("rms", "correlation")
)

# All four metrics (RMS, correlation, Frechet, DTW)
result_all <- trajectory_similarity(sap,
  metrics = "all"
)

# Access results
result$summary # per-label statistics
result$tests # statistical tests
} # }
```
