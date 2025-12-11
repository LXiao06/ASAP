# Plot Clusters in Time Series Data

Generic function for plotting clusters in time series data, primarily
designed for visualizing bird song syllable clusters and their manual
labels. Supports both numeric cluster IDs (pre manual labeling) and
alphabetic syllable labels (post manual labeling).

## Usage

``` r
plot_clusters(x, ...)

# S3 method for class 'matrix'
plot_clusters(x, labels = NULL, cluster_colors = NULL, main = NULL, ...)

# S3 method for class 'Sap'
plot_clusters(
  x,
  data_type = c("segment", "syllable"),
  label_type = c("pre", "post"),
  time_resolution = 1000,
  cluster_colors = NULL,
  sample_percent = NULL,
  balanced = FALSE,
  labels = NULL,
  motif_clusters = NULL,
  ordered = FALSE,
  descending = TRUE,
  seed = 222,
  cores = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  Object to plot clusters from. Can be either:

  - A matrix with cluster/syllable data

  - A Sap object containing bird song analysis data at motif and
    syllable levels

- ...:

  Additional arguments passed between methods

- labels:

  Optional vector of experimental condition labels to subset (e.g.,
  "BL", "Post", "Rec")

- cluster_colors:

  Optional named vector of colors for clusters

- main:

  Plot title (for matrix method)

- data_type:

  Type of data to analyze (for Sap objects):

  - "segment": Uses segment-level features and clusters

  - "syllable": Uses syllable-level features and clusters

- label_type:

  Type of labels to visualize (for Sap objects):

  - "pre": Shows numerical cluster IDs from automatic clustering

  - "post": Shows alphabetic syllable labels assigned via manual_label()

- time_resolution:

  Number of time points for plotting (default: 1000)

- sample_percent:

  Optional percentage of data to sample

- balanced:

  Whether to balance samples across labels

- motif_clusters:

  Optional vector of specific clusters to include

- ordered:

  Whether to order motifs by UMAP coordinates

- descending:

  Direction of UMAP-based ordering

- seed:

  Random seed for reproducibility

- cores:

  Number of cores for parallel processing

- verbose:

  Whether to print progress messages

## Value

Invisibly returns x

## Details

The function provides two main visualization approaches:

1.  Pre-labeling visualization (label_type = "pre"):

    - Shows numerical cluster IDs from automatic clustering

    - Uses data from x\$features\$segment/syllable\$feat.embeds

    - Clusters are represented by numbers (1, 2, 3, etc.)

    - Useful for evaluating automatic clustering results

2.  Post-labeling visualization (label_type = "post"):

    - Shows alphabetic syllable labels from manual annotation

    - Uses data from x\$syllables (created by manual_label())

    - Syllables are represented by letters (a, b, c, etc.)

    - Useful for viewing manual syllable classifications

For matrix input:

- Must have column names as labels

- Requires a 'time_window' attribute

- Can contain either numeric clusters or character syllable labels

For Sap objects:

- data_type determines the feature level for analysis:

  - "segment": Fine-grained analysis of song segments

  - "syllable": Analysis at the syllable level

- Supports experimental condition labels (BL, Post, Rec)

- Can order motifs by UMAP coordinates for pattern visualization

- Allows balanced sampling across conditions

## See also

- manual_label() for creating syllable labels

- auto_label() for automatic cluster generation

## Examples

``` r
if (FALSE) { # \dontrun{
# Matrix method example
mat <- matrix(sample(1:5, 1000, replace = TRUE), ncol = 10)
attr(mat, "time_window") <- 1.2
colnames(mat) <- rep(c("BL", "Post"), each = 5)
plot_clusters(mat)

# Sap object - Pre-labeling (numerical clusters)
# Basic cluster visualization
plot_clusters(sap,
            data_type = "syllable",
            label_type = "pre")

# Ordered by UMAP with specific conditions
plot_clusters(sap,
            data_type = "syllable",
            label_type = "pre",
            ordered = TRUE,
            labels = c("BL", "Post"))

# Sap object - Post-labeling (syllable letters)
# After running manual_label()
plot_clusters(sap,
            data_type = "syllable",
            label_type = "post")

# Balanced sampling across conditions
plot_clusters(sap,
            data_type = "syllable",
            label_type = "post",
            balanced = TRUE,
            labels = c("BL", "Post", "Rec"))
} # }
```
