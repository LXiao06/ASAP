# Plot Trajectory Convergence Results

Creates multi-panel visualization of trajectory convergence metrics
showing how individual trial trajectories compare to a reference path
across labels.

## Usage

``` r
plot_trajectory_convergence(x, ...)

# Default S3 method
plot_trajectory_convergence(
  x,
  palette = "Set1",
  max_annotations = 10,
  similarity_baseline = c("result", "reference", "zero"),
  similarity_scale_multiplier = NULL,
  ...
)

# S3 method for class 'Sap'
plot_trajectory_convergence(
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
  [`trajectory_convergence()`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_convergence.md)
  or a SAP object with trajectory convergence results

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
Pointwise Convergence, Shape Correlation, Path-Shape Convergence, and
Timing-Adjusted Convergence when those metrics are present. For small
label sets (≤6 labels), uses violin + box plots. For larger label sets,
uses mean ± SE trend lines. Distance similarities are recalculated from
raw distances using metric-specific reference scales when available.

Reference label is included in results and plots for context but
excluded from statistical test annotations (reference defines target,
not convergence).

## See also

[`trajectory_convergence`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_convergence.md)
for analysis computation

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot trajectory convergence
result <- trajectory_convergence(sap)
plot_trajectory_convergence(result)

# From SAP object directly
plot_trajectory_convergence(sap,
  segment_type = "motifs"
)
} # }
```
