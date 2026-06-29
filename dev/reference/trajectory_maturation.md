# Compute Trajectory Maturation Scores

Calculate composite developmental maturity scores by combining
trajectory similarity (similarity to reference) and within-rendition
variability (consistency). Higher scores indicate more mature,
crystallized trajectories.

## Usage

``` r
trajectory_maturation(x, ...)

# Default S3 method
trajectory_maturation(
  x,
  similarity_metric = c("rms", "frechet", "dtw", "correlation"),
  variability_metric = c("dispersion", "orthogonal_rms", "parallel_rms"),
  score_type = c("maturation", "stability", "both"),
  invert_variability = TRUE,
  epsilon = 0.01,
  normalize_variability = c("none", "reference"),
  reference_label = NULL,
  norm_epsilon = 1e-06,
  scale_method = c("minmax", "zscore", "none"),
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
trajectory_maturation(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  similarity_metric = c("rms", "frechet", "dtw", "correlation"),
  variability_metric = c("dispersion", "orthogonal_rms", "parallel_rms"),
  score_type = c("maturation", "stability", "both"),
  invert_variability = TRUE,
  epsilon = 0.01,
  normalize_variability = c("none", "reference"),
  reference_label = NULL,
  norm_epsilon = 1e-06,
  scale_method = c("minmax", "zscore", "none"),
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to analyze: a data.frame with similarity and variability
  results, or a SAP object

- ...:

  Additional arguments

- similarity_metric:

  Character. Which similarity metric to use: "rms" (default), "frechet",
  "dtw", or "correlation"

- variability_metric:

  Character. Which variability metric to use: "dispersion" (default),
  "orthogonal_rms", or "parallel_rms"

- score_type:

  Character vector. Which scores to compute: "maturation" (default),
  "stability", or "both"

- invert_variability:

  Logical. If TRUE (default), higher variability reduces the score

- epsilon:

  Numeric. Small constant to avoid division by zero in stability_index
  (default: 0.01)

- normalize_variability:

  Character. How to normalize variability for cross-animal comparison:
  "none" (default), or "reference" (normalize to reference label). When
  "reference", variability is divided by the mean variability at the
  reference label, making values interpretable as "X times more variable
  than reference"

- reference_label:

  Character. Label to use as normalization reference. If NULL (default),
  uses the last label (typically adult). Only used when
  normalize_variability = "reference"

- norm_epsilon:

  Numeric. Small constant added to reference variability to avoid
  division by zero (default: 1e-6)

- scale_method:

  Character. How to scale variability: "minmax" (default), "zscore", or
  "none"

- verbose:

  Logical. Print progress messages (default: TRUE)

- segment_type:

  For SAP objects: Type of segments ("motifs", "syllables", "bouts", or
  "segments")

## Value

For data.frame: A data.frame with original data plus score columns:

- `variability_raw`: Original variability metric value (renamed from
  input)

- `variability_scaled`: Scaled variability (for ML)

- `variability_normalized`: Normalized to reference label (if
  normalize_variability = "reference"). Values represent "X times the
  reference variability". For example, 2.0 means twice as variable as
  the reference label (typically adult)

- `maturation_score`: Computed maturation score (if requested)

- `stability_index`: Computed stability index (if requested)

Metadata attributes include: variability_metric, scale_method,
similarity_metric, score_type, invert_variability,
normalize_variability, reference_label

For SAP objects: Updated object with scores stored in
`sap$features[[feature_type]]$maturation_scores`

## Details

**Score Types:**

- Maturation score: similarity × (1 - scaled_variability). Combines
  convergence to reference with consistency into a unified developmental
  score

- Stability index: similarity / (scaled_variability + epsilon). Ratio
  form that shows relative stability independent of scale

**Prerequisites:** Must first run
[`trajectory_similarity`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_similarity.md)
and either
[`trajectory_dispersion`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_dispersion.md)
or
[`trajectory_path_deviation`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_path_deviation.md).

**Variability Metrics:**

- `dispersion`: Mean distance to centroid trajectory (from
  trajectory_dispersion)

- `orthogonal_rms`: Perpendicular deviation from mean trajectory (from
  trajectory_path_deviation)

- `parallel_rms`: Along-trajectory timing/phase deviation (from
  trajectory_path_deviation)

## See also

[`plot_trajectory_maturation`](https://lxiao06.github.io/ASAP/dev/reference/plot_trajectory_maturation.md)
for visualization,
[`trajectory_similarity`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_similarity.md)
for similarity computation,
[`trajectory_dispersion`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_dispersion.md)
for dispersion-based variability

## Examples

``` r
if (FALSE) { # \dontrun{
# After running trajectory_similarity and trajectory_dispersion
sap <- trajectory_maturation(sap,
  segment_type = "motifs",
  similarity_metric = "rms",
  score_type = "both"
)

# With reference-based normalization (for cross-animal ML)
sap <- trajectory_maturation(sap,
  segment_type = "motifs",
  normalize_variability = "reference",
  reference_label = "90"  # Adult as reference
)

# Access results
scores <- sap$features$motif$maturation_scores
# variability_normalized = 2.0 means "2x adult variability"
} # }
```
