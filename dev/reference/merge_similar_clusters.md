# Merge Similar Clusters Based on UMAP Proximity

Internal function to merge clusters that are close in UMAP space.

## Usage

``` r
merge_similar_clusters(data, umap_threshold = NULL)
```

## Arguments

- data:

  Data frame containing cluster information

- umap_threshold:

  Distance threshold for merging clusters

## Value

Data frame with merged cluster assignments
