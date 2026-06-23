# Plot PCA coordinates of feature embeddings

Wrapper for
[`plot_umap`](https://lxiao06.github.io/ASAP/dev/reference/plot_umap.md)
to plot PCA coordinates (PC1 and PC2) instead of UMAP. Inherits all
arguments from
[`plot_umap`](https://lxiao06.github.io/ASAP/dev/reference/plot_umap.md).

## Usage

``` r
plot_PC(x, ...)

# Default S3 method
plot_PC(x, dims = c("PC1", "PC2"), ...)

# S3 method for class 'Sap'
plot_PC(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  dims = c("PC1", "PC2"),
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A data frame or a SAP object.

- ...:

  Additional arguments passed to
  [`plot_umap`](https://lxiao06.github.io/ASAP/dev/reference/plot_umap.md).

- dims:

  Character vector of length 2 specifying the PCA dimensions to plot.
  Default is `c("PC1", "PC2")`.

- segment_type:

  For SAP objects: Type of segments to visualize ('motifs', 'syllables',
  'bouts', 'segments')

- verbose:

  For SAP objects: Whether to print progress messages
