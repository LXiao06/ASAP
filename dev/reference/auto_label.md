# Automatic Syllable Labeling in Song Motif

Performs two-stage clustering of bird song segments to identify
syllables:

1.  Temporal clustering using weighted time and duration features

2.  UMAP-based refinement of temporal clusters Finally merges similar
    clusters based on UMAP proximity.

## Usage

``` r
auto_label(
  x,
  eps_time = 0.1,
  eps_umap = 0.6,
  min_pts = 5,
  outlier_threshold = 0.01,
  umap_threshold = 1,
  weight_time = 4,
  ...
)
```

## Arguments

- x:

  A Sap object containing bird song analysis data

- eps_time:

  Epsilon parameter for temporal DBSCAN clustering (default: 0.1)

- eps_umap:

  Epsilon parameter for UMAP-based DBSCAN clustering (default: 0.6)

- min_pts:

  Minimum points parameter for DBSCAN clustering (default: 5)

- outlier_threshold:

  Threshold for removing small clusters (default: 0.01)

- umap_threshold:

  Distance threshold for merging similar clusters (default: 1)

- weight_time:

  Weight factor for temporal features (default: 4)

- ...:

  Additional arguments (not currently used)

## Value

Returns the input Sap object with updated syllable clusters in:
x\$features\$syllable\$feat.embeds

## Details

The function performs clustering in multiple stages:

1.  Temporal Clustering:

    - Combines time position and duration information

    - Uses weighted features (controlled by weight_time):

      - weight_time \> 1: Emphasizes temporal position over duration

      - weight_time = 1: Equal weighting

      - weight_time \< 1: Emphasizes duration over temporal position

2.  UMAP-based Refinement:

    - Further splits temporal clusters based on UMAP coordinates

    - Helps distinguish syllables with similar timing but different
      acoustic features

3.  Cluster Cleaning:

    - Removes small clusters (controlled by outlier_threshold)

    - Merges similar clusters based on UMAP proximity (controlled by
      umap_threshold)

Parameter Tuning Guidelines:

- eps_time: Controls temporal separation sensitivity

  - Smaller values create more temporal splits

  - Typical range: 0.05-0.2

- eps_umap: Controls acoustic feature sensitivity

  - Smaller values create more acoustic splits

  - Typical range: 0.4-0.8

- weight_time: Controls temporal vs duration importance

  - Default (4) weights time 4x more than duration

  - Increase for more temporal separation

  - Decrease for more duration-based separation

## See also

- manual_label() for manual refinement of automatic clusters

- plot_cluster() for visualizing clustering results

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with default parameters
sap <- auto_label(sap)

# Emphasize temporal separation
sap <- auto_label(sap,
                 eps_time = 0.08,
                 weight_time = 6)

# More acoustic feature sensitivity
sap <- auto_label(sap,
                 eps_umap = 0.4,
                 min_pts = 3)

# Stricter cluster cleaning
sap <- auto_label(sap,
                 outlier_threshold = 0.02,
                 umap_threshold = 0.8)
} # }
```
