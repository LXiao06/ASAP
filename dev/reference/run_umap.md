# Run UMAP Dimensionality Reduction

A generic function to perform UMAP dimensionality reduction on feature
data.

## Usage

``` r
run_umap(x, ...)

# Default S3 method
run_umap(
  x,
  metadata_cols = NULL,
  scale = TRUE,
  n_neighbors = 15,
  n_components = 2,
  min_dist = 0.1,
  seed = 222,
  n_threads = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
run_umap(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  data_type = c("spectral_feature", "spectrogram", "traj_mat"),
  label = NULL,
  scale = TRUE,
  n_neighbors = 20,
  n_components = 2,
  min_dist = 0.1,
  seed = 222,
  n_threads = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to analyze, either a data frame or SAP object

- ...:

  Additional arguments passed to specific methods

- metadata_cols:

  Column indices for metadata (for default method)

- scale:

  Whether to scale features before UMAP

- n_neighbors:

  Number of neighbors (default: 15)

- n_components:

  Number of output dimensions (default: 2)

- min_dist:

  Minimum distance parameter (default: 0.1)

- seed:

  Random seed for reproducibility

- n_threads:

  Number of computation threads

- verbose:

  Whether to print progress messages

- segment_type:

  For SAP objects: Type of segments to analyze ('motifs', 'syllables',
  'bouts', 'segments')

- data_type:

  For SAP objects: Type of feature data
  ('spectral_feature','spectrogram', 'traj_mat')

- label:

  For SAP objects: Specific label to filter data

## Value

For default method: Matrix of UMAP coordinates For SAP objects: Updated
SAP object with UMAP coordinates stored in features slot

## Details

This generic function supports UMAP analysis through two methods:

- Default method for feature data frames

- SAP object method for organized song features

## Examples

``` r
if (FALSE) { # \dontrun{
# Run UMAP on feature data frame
coords <- run_umap(features,
                   metadata_cols = c(1:5),
                   n_neighbors = 15)

# Run UMAP on SAP object
sap_obj <- run_umap(sap_object,
                    segment_type = "motifs",
                    data_type = "spectral_feature")

# UMAP with specific parameters
coords <- run_umap(features,
                   metadata_cols = 1:3,
                   scale = TRUE,
                   n_neighbors = 20,
                   seed = 123)

# UMAP with label filtering
sap_obj <- run_umap(sap_obj,
                    segment_type = "syllables",
                    data_type = "spectral_feature",
                    label = "a")
} # }
```
