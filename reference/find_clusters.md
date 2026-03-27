# Find Clusters in Feature Data

Performs cluster analysis on feature data using shared nearest neighbor
(SNN) clustering.

## Usage

``` r
find_clusters(x, ...)

# Default S3 method
find_clusters(
  x,
  metadata_cols,
  k.param = 20,
  prune.SNN = 1/15,
  n.pcs = 20,
  resolution = 0.2,
  n.start = 10,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
find_clusters(
  x,
  segment_type = c("motifs", "syllables", "segments"),
  data_type = c("spectral_feature", "spectrogram"),
  label = NULL,
  k.param = 20,
  prune.SNN = 1/15,
  n.pcs = NULL,
  resolution = NULL,
  n.start = 10,
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

  For default method: Column indices for metadata

- k.param:

  Number of nearest neighbors (default: 20)

- prune.SNN:

  Pruning threshold for SNN graph (default: 1/15)

- n.pcs:

  Number of principal components to use

- resolution:

  Resolution parameter for clustering

- n.start:

  Number of random starts (default: 10)

- verbose:

  Whether to print progress messages (default: TRUE)

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'segments')

- data_type:

  For SAP objects: Type of feature data ('spectral_feature',
  "spectrogram")

- label:

  For SAP objects: Specific label to filter data

## Value

For default method: Data frame with metadata columns and cluster
assignments For SAP objects: Updated object with clustering results in
features slot

## Details

For feature data frames:

- Separates metadata and feature columns

- Finds nearest neighbors using PCA

- Constructs SNN graph

- Performs community detection

For SAP objects:

- Supports multiple segment types

- Optional label filtering

- Stores results in features slot

- Updates feature embeddings

The clustering approach is similar to that used in Seurat V3,
implementing:

- PCA-based neighbor finding

- SNN graph construction

- Louvain community detection

## See also

[`run_umap`](https://lxiao06.github.io/ASAP/reference/run_umap.md) for
visualization of clusters

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic clustering of feature data
clusters <- find_clusters(features,
                         metadata_cols = c(1:5))

# Clustering with custom parameters
clusters <- find_clusters(features,
                         metadata_cols = c(1:5),
                         k.param = 30,
                         resolution = 0.3)

# Cluster SAP object features
sap_obj <- find_clusters(sap_object,
                        segment_type = "motifs",
                        data_type = "spectral_feature")

# Label-specific clustering
sap_obj <- find_clusters(sap_object,
                        segment_type = "syllables",
                        label = "a",
                        resolution = 0.4)
} # }
```
