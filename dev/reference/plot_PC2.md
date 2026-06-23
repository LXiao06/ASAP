# Plot PCA coordinates of trajectories over time

Wrapper for
[`plot_umap2`](https://lxiao06.github.io/ASAP/dev/reference/plot_umap2.md)
to plot PCA coordinates (PC1 and PC2) instead of UMAP. Inherits all
arguments from
[`plot_umap2`](https://lxiao06.github.io/ASAP/dev/reference/plot_umap2.md).

## Usage

``` r
plot_PC2(x, ...)

# Default S3 method
plot_PC2(x, dims = c("PC1", "PC2"), ...)

# S3 method for class 'Sap'
plot_PC2(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  data_type = c("feat.embeds", "traj.embeds"),
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
  [`plot_umap2`](https://lxiao06.github.io/ASAP/dev/reference/plot_umap2.md).

- dims:

  Character vector of length 2 specifying the PCA dimensions to plot.
  Default is `c("PC1", "PC2")`.

- segment_type:

  For SAP objects: Type of segments to visualize ('motifs', 'syllables',
  'bouts', 'segments')

- data_type:

  For SAP objects: Type of embedding data ('feat.embeds', 'traj.embeds')

- verbose:

  For SAP objects: Whether to print progress messages
