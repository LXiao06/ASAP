# Plot Trajectory Maturation Scores

Visualize developmental maturation scores across labels or age groups.
Uses ASAP's standard trajectory plotting style (violin + boxplot
panels).

## Usage

``` r
plot_trajectory_maturation(x, ...)

# Default S3 method
plot_trajectory_maturation(
  x,
  score_cols = NULL,
  palette = "Set1",
  max_annotations = 10,
  ...
)

# S3 method for class 'Sap'
plot_trajectory_maturation(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  score_cols = NULL,
  palette = "Set1",
  max_annotations = 10,
  ...
)
```

## Arguments

- x:

  An object to plot: a data.frame with maturation scores or a SAP object

- ...:

  Additional arguments

- score_cols:

  Character vector. Which score columns to plot (default: auto-detect
  from available scores)

- palette:

  Character. RColorBrewer palette name (default: "Set1")

- max_annotations:

  Numeric. Maximum number of significance brackets to show (default: 10)

- segment_type:

  For SAP objects: Type of segments ("motifs", "syllables", "bouts", or
  "segments")

## Value

A ggplot2 or patchwork object

## Details

This function uses ASAP's existing plotting infrastructure from
trajectory.R:

- `sort_labels()`: Intelligent numeric-aware label sorting

- `make_pal()`: RColorBrewer palette generation

- `panel()`: Violin + boxplot panel creation

- `trend_panel()`: Line plot for \>6 groups

- `fmt_p()`: P-value formatting

- `brackets()`: Significance bracket computation

- `add_brackets()`: Significance annotation overlay

The resulting plots match the style of other ASAP trajectory
visualizations.

## See also

[`trajectory_maturation`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_maturation.md)
for score computation

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot all computed scores
plot_trajectory_maturation(sap, segment_type = "motifs")

# Plot specific score
plot_trajectory_maturation(scores_df, score_cols = "maturation_score")
} # }
```
