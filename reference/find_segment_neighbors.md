# Find Nearest Neighbors for Segments

Internal function to compute nearest neighbors and SNN graph.

## Usage

``` r
find_segment_neighbors(
  segments_data,
  k.param = 20,
  prune.SNN = 1/15,
  compute.SNN = TRUE,
  n.pcs = 10
)
```

## Arguments

- segments_data:

  Matrix or data frame of segment features

- k.param:

  Number of nearest neighbors

- prune.SNN:

  SNN pruning threshold

- compute.SNN:

  Whether to compute SNN graph

- n.pcs:

  Number of PCs to use

## Value

List containing neighbor graphs
