# Find Clusters Using Modularity Optimization

Internal function to perform cluster detection.

## Usage

``` r
find_segment_clusters(
  neighbor_graphs,
  modularity.fxn = 1,
  resolution = 0.8,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  verbose = TRUE
)
```

## Arguments

- neighbor_graphs:

  List of neighbor graphs

- modularity.fxn:

  Modularity function type

- resolution:

  Resolution parameter

- algorithm:

  Clustering algorithm

- n.start:

  Number of random starts

- n.iter:

  Maximum iterations

- random.seed:

  Random seed

- group.singletons:

  Whether to group singleton clusters

- verbose:

  Print progress messages

## Value

Data frame with cluster assignments
