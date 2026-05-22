# Feature Statistics Summary and Analysis

A generic function to calculate, print, and plot feature statistics. The
default method accepts a pre-computed segment statistics data frame,
while the `Sap` method extracts feature statistics from a `Sap` object
and forwards them to the default method.

## Usage

``` r
feature_stats(x, ...)

# Default S3 method
feature_stats(
  x,
  feature_name = "Feature",
  labels = NULL,
  run_anova = TRUE,
  plot = TRUE,
  plot_type = c("both", "boxplot", "violin"),
  palette = "Set1",
  jitter = TRUE,
  point_alpha = 0.35,
  facet_segments = FALSE,
  ncol = NULL,
  ...
)

# S3 method for class 'Sap'
feature_stats(
  x,
  feature = c("fund_freq", "wiener_entropy"),
  segment_type = c("motifs", "syllables", "segments"),
  ...
)
```

## Arguments

- x:

  An object containing statistics: either a data frame (e.g.
  `sap$features$motif$stats_fund_freq`) or a `Sap` object.

- ...:

  Additional arguments passed to methods.

- feature_name:

  For default method: Character string used as the y-axis label and plot
  title (default: `"Feature"`).

- labels:

  Optional character vector of labels to filter the data.

- run_anova:

  Logical. Whether to perform a one-way ANOVA (default: `TRUE`).

- plot:

  Logical. Whether to generate the ggplot visualization (default:
  `TRUE`).

- plot_type:

  Character. Style of visual layer: `"both"` (default), `"boxplot"`, or
  `"violin"`.

- palette:

  Character. RColorBrewer palette name for fills/colors (default:
  `"Set1"`).

- jitter:

  Logical. Overlay individual rendition data points (default: `TRUE`).

- point_alpha:

  Numeric `[0, 1]`. Transparency of jittered points (default: `0.35`).

- facet_segments:

  Logical. Facet the plot by `segment_id` (default: `FALSE`).

- ncol:

  Integer. Number of columns for facets if `facet_segments = TRUE`.

- feature:

  For SAP objects: Feature to analyze, either `"fund_freq"` or
  `"wiener_entropy"`.

- segment_type:

  For SAP objects: Type of segments to extract, either `"motifs"`,
  `"syllables"`, or `"segments"`.

## Value

A list (invisibly) with components:

- `summary`:

  A tibble of per-label summary statistics.

- `anova`:

  ANOVA result tibble, or `NULL` if not run.

- `plot`:

  The ggplot object, or `NULL` if `plot = FALSE`.

## See also

[`anova_analysis`](https://lxiao06.github.io/ASAP/dev/reference/anova_analysis.md),
[`refine_FF`](https://lxiao06.github.io/ASAP/dev/reference/refine_FF.md),
[`refine_sh`](https://lxiao06.github.io/ASAP/dev/reference/refine_sh.md)
