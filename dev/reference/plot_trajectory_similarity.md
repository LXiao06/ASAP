# Plot Trajectory Similarity Results

Creates multi-panel visualization of trajectory similarity metrics
showing how individual trial trajectories compare to a reference path
across labels.

## Usage

``` r
plot_trajectory_similarity(x, ...)

# Default S3 method
plot_trajectory_similarity(
  x,
  palette = "Set1",
  max_annotations = 10,
  similarity_baseline = c("result", "reference", "zero"),
  similarity_scale_multiplier = NULL,
  ...
)

# S3 method for class 'Sap'
plot_trajectory_similarity(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  palette = "Set1",
  max_annotations = 10,
  similarity_baseline = c("result", "reference", "zero"),
  similarity_scale_multiplier = NULL,
  ...
)
```

## Arguments

- x:

  A list returned by
  [`trajectory_similarity()`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_similarity.md)
  or a SAP object with trajectory similarity results

- ...:

  Additional arguments passed to specific methods

- palette:

  Character. Color palette name for ggplot2 (default: "Set1")

- max_annotations:

  Integer. Maximum number of pairwise significance annotations to
  display per panel (default: 10)

- similarity_baseline:

  Similarity transform to use for plotted distance metrics. `"result"`
  uses the baseline stored in `x`; `"reference"` treats distances within
  reference-label variability as converged; `"zero"` uses
  `exp(-normalized_distance)`.

- similarity_scale_multiplier:

  Optional multiplier for metric-specific reference scales when plotting
  distance similarities. Defaults to the multiplier stored in `x`, or
  `1` when unavailable.

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'bouts',
  'segments')

## Value

A patchwork object combining all metric panels (returned invisibly)

## Details

Creates one panel for each computed metric. Panels are ordered as
Pointwise Similarity, Shape Correlation, Path-Shape Similarity, and
Timing-Adjusted Similarity when those metrics are present. For small
label sets (\\\le\\6 labels), uses violin + box plots. For larger label
sets, uses mean \\\pm\\ SE trend lines. Distance similarities are
recalculated from raw distances using metric-specific reference scales
when available.

Reference label is included in results and plots for context but
excluded from statistical test annotations (reference defines target,
not similarity).

## See also

[`trajectory_similarity`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_similarity.md)
for analysis computation

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot trajectory similarity
result <- trajectory_similarity(sap)
plot_trajectory_similarity(result)

# From SAP object directly
plot_trajectory_similarity(sap,
  segment_type = "motifs"
)

# With all four metrics
result_all <- trajectory_similarity(sap, metrics = "all")
plot_trajectory_similarity(result_all)
} # }
```
