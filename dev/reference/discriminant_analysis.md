# Canonical Discriminant Analysis with Robustness Testing

Performs Canonical Discriminant Analysis (CDA/LDA) on multivariate
feature data (such as trajectory similarities, acoustic feature
statistics, or PCA embeddings) to quantify and visualize group
separation.

## Usage

``` r
discriminant_analysis(
  data,
  group_col = "label",
  feature_cols = NULL,
  reference_group = NULL,
  scale = TRUE,
  n_perm = 500,
  cv = TRUE,
  seed = 42,
  plot = TRUE,
  save_plot = FALSE,
  plot_dir = NULL,
  palette = "Set1",
  ellipse = TRUE,
  ellipse_level = 0.95,
  alpha = 0.4,
  point_size = 2,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame containing group labels and numeric feature columns

- group_col:

  Character. Name of the grouping variable (default: `"label"`)

- feature_cols:

  Character vector of feature column names. If `NULL` (default), all
  numeric columns (excluding common metadata columns) are automatically
  used

- reference_group:

  Optional character. Reference group label (e.g., `"tutor"`) to compute
  multivariate Mahalanobis distances from all groups to the reference

- scale:

  Logical. Whether to standardize features to zero-mean and
  unit-variance prior to analysis (default: `TRUE`)

- n_perm:

  Integer. Number of label permutations for empirical robustness testing
  (default: `500`). Set to `0` to skip permutation testing

- cv:

  Logical. Whether to perform Leave-One-Out Cross-Validation (default:
  `TRUE`)

- seed:

  Integer. Random seed for reproducible permutation testing (default:
  `42`)

- plot:

  Logical. Whether to generate and display the canonical separation plot
  (default: `TRUE`)

- save_plot:

  Logical. Whether to save the generated plot to disk (default: `FALSE`)

- plot_dir:

  Optional character. Directory to save the plot when `save_plot = TRUE`

- palette:

  Character. RColorBrewer palette name for group colors (default:
  `"Set1"`)

- ellipse:

  Logical. Whether to draw confidence ellipses around group clusters
  (default: `TRUE`)

- ellipse_level:

  Numeric. Confidence level for ellipses (default: `0.95`)

- alpha:

  Numeric. Transparency for individual data points (default: `0.4`)

- point_size:

  Numeric. Size of individual data points (default: `2`)

- verbose:

  Logical. Whether to print analytical summary to console (default:
  `TRUE`)

- ...:

  Additional arguments passed to internal plot helper

## Value

A list of class `"DiscriminantAnalysis"` containing:

- `scores`:

  Data frame of projected canonical variate scores per observation.

- `centroids`:

  Data frame of group centroids in canonical space.

- `loadings`:

  Standardized canonical coefficients for each feature.

- `prop_trace`:

  Proportion of between-group variance explained by each canonical axis.

- `mahalanobis`:

  Matrix of pairwise Mahalanobis distances between group centroids.

- `reference_distances`:

  Distances from each group to `reference_group` (if specified).

- `loocv_accuracy`:

  Overall cross-validated classification accuracy.

- `confusion_matrix`:

  Confusion matrix of observed versus LOOCV-predicted classes.

- `permutation`:

  List containing `n_perm`, empirical `p_value`, and null distribution.

- `manova`:

  Multivariate ANOVA (Wilks' Lambda) summary.

- `plot`:

  The generated `ggplot` object (or `NULL` if `plot = FALSE`).

- `features_used`:

  Character vector of feature column names used in the analysis.

## Details

Canonical Discriminant Analysis finds linear combinations of features
(Canonical Variates: `CV1`, `CV2`, ...) that maximize between-group
variance relative to within-group variance.

**Robustness and Validation:**

- **Leave-One-Out Cross-Validation (LOOCV):** Evaluates unbiased
  classification accuracy without sample-size inflation.

- **Permutation Test:** Randomly permutes group labels `n_perm` times to
  construct an empirical null distribution and compute a non-parametric
  permutation p-value for group discriminability.

- **Mahalanobis Distance:** Computes multivariate distances between
  group centroids accounting for feature covariance.

## See also

[`trajectory_similarity`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_similarity.md),
[`feature_stats`](https://lxiao06.github.io/ASAP/dev/reference/feature_stats.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Using trajectory similarity results
sim_data <- sap$features$motif$trajectory_similarity$similarity
res <- discriminant_analysis(
  data = sim_data,
  group_col = "label",
  reference_group = "tutor",
  n_perm = 500
)

# Inspect classification accuracy and confusion matrix
res$loocv_accuracy
res$confusion_matrix

# Access or modify the plot
res$plot + ggplot2::theme_minimal()
} # }
```
