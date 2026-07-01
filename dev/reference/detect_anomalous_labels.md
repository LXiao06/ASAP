# Detect Anomalous Time Labels

Identifies labels with abnormal variability, low maturation scores, or
unexpected deviations from trajectory trends. Useful for detecting
developmental anomalies or time points requiring further investigation.

## Usage

``` r
detect_anomalous_labels(
  sap,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  methods = c("trajectory", "multivariate", "lof"),
  extreme_threshold = 2,
  isolation_quantile = 0.99,
  min_sample_size = 3,
  lof_k = NULL,
  cov_regularization = 0.5,
  verbose = TRUE
)
```

## Arguments

- sap:

  SAP object with computed maturation scores.

- segment_type:

  Type of segments: `"motifs"`, `"syllables"`, `"bouts"`, or
  `"segments"`.

- methods:

  Vector of detection methods to use:

  - `"trajectory"`: Fits a loess trend to each metric and z-scores the
    residuals. Flags labels with \|z-score\| \> extreme_threshold.

  - `"multivariate"`: Mahalanobis distance on the loess residual
    z-scores. Falls back to raw feature means when trajectory is not
    run. Requires MASS package.

  - `"lof"`: Local Outlier Factor (non-parametric). Requires the
    `dbscan` package.

- extreme_threshold:

  Z-score threshold for detecting extreme values (default: 2). Flags
  when \|z-score\| \> extreme_threshold.

- isolation_quantile:

  Chi-squared quantile for multivariate threshold (default: 0.99).
  Increase to flag fewer labels; decrease to flag more.

- min_sample_size:

  Minimum samples per label to include (default: 3).

- lof_k:

  Integer. Number of nearest neighbours for LOF. Defaults to
  `max(3, min(10, floor(n/4)))` where `n` is the number of labels.

- cov_regularization:

  Numeric. Regularization constant for covariance shrinkage (default:
  0.5).

- verbose:

  Logical. Print progress messages (default: TRUE).

## Value

A list of class `"label_anomaly"` with:

- `label_scores`: Data frame with anomaly scores per label.

- `anomalous_labels`: Vector of flagged labels.

- `method_flags`: Data frame with per-method flags (`trajectory_flag`,
  `multivariate_flag`, `lof_flag`, `total_flags`, `is_anomalous`).

- `parameters`: List of detection parameters used.

## See also

[`plot_label_anomaly`](https://lxiao06.github.io/ASAP/dev/reference/plot_label_anomaly.md)
for visualising results.

## Examples

``` r
if (FALSE) { # \dontrun{
# Detect anomalous labels using all three methods (default)
anomaly_results <- detect_anomalous_labels(sap, segment_type = "motifs")

# View flagged labels
print(anomaly_results$anomalous_labels)

# View detailed scores
View(anomaly_results$label_scores)

# Trajectory method only
anomaly_results <- detect_anomalous_labels(
  sap,
  segment_type = "motifs",
  methods = "trajectory"
)

# Custom LOF neighbourhood size
anomaly_results <- detect_anomalous_labels(
  sap,
  segment_type = "motifs",
  lof_k = 5
)
} # }
```
