# Run Principal Component Analysis

A generic function to perform PCA with support for multiple methods and
large-scale data processing.

## Usage

``` r
run_pca(x, ...)

# Default S3 method
run_pca(
  x,
  method = c("irlba", "base", "parallel"),
  n_components = 50,
  n_cores = NULL,
  diagnostic_plots = TRUE,
  output_scores = TRUE,
  scale_scores = FALSE,
  ...
)

# S3 method for class 'Sap'
run_pca(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  data_type = "traj_mat",
  method = c("irlba", "base", "parallel"),
  n_components = 50,
  n_cores = NULL,
  diagnostic_plots = TRUE,
  scale_scores = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to analyze, either a matrix/data frame or a SAP object

- ...:

  Additional arguments passed to specific methods

- method:

  PCA method ('irlba', 'base', or 'parallel')

- n_components:

  Number of principal components

- n_cores:

  Number of cores for parallel processing

- diagnostic_plots:

  Whether to create diagnostic plots

- output_scores:

  Whether to return PC scores

- scale_scores:

  Whether to scale PC scores

- segment_type:

  For SAP objects: Type of segments to analyze ('motifs', 'syllables',
  'bouts', 'segments')

- data_type:

  For SAP objects: Type of data to analyze

- verbose:

  For SAP objects: Whether to print progress

## Value

For default method: PCA results or PC scores matrix For SAP objects:
Updated SAP object with PCA results stored in features slot

## Details

This generic function supports PCA through two methods:

- Default method for matrices with multiple PCA implementations

- SAP object method for organized trajectory data

## Examples

``` r
if (FALSE) { # \dontrun{
# Run PCA on matrix
pca <- run_pca(matrix, method = "irlba", n_components = 50)

# Run PCA on SAP object
sap_obj <- run_pca(sap_object,
                   segment_type = "motifs",
                   data_type = "traj_mat")
} # }
```
