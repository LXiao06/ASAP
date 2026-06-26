# Plot Individual Trial Scores

Creates a dot plot of per-trial similarity or maturation scores, ordered
by `.source_row` within each label/day. Helps visualise the sequential
progression of scores across individual renditions within and across
days. An optional running-mean trace (default window: 50 trials)
highlights the local trend.

## Usage

``` r
plot_trajectory_trials(x, ...)

# Default S3 method
plot_trajectory_trials(
  x,
  score_col = "rms_similarity",
  labels = NULL,
  running_mean = TRUE,
  window_size = 50,
  palette = "Set1",
  ...
)

# S3 method for class 'Sap'
plot_trajectory_trials(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  data_type = c("similarity", "maturation"),
  score_col = "rms_similarity",
  labels = NULL,
  running_mean = TRUE,
  window_size = 50,
  palette = "Set1",
  ...
)
```

## Arguments

- x:

  An object to plot: a data frame with per-trial scores, or a SAP object
  with pre-computed trajectory results.

- ...:

  Additional arguments passed to specific methods.

- score_col:

  Character. Column name of the score to plot (default:
  "rms_similarity").

- labels:

  Character vector of labels (days) to include. Default NULL shows the
  first and last label when labels are sorted numerically.

- running_mean:

  Logical. If TRUE (default), overlay a running-mean trace.

- window_size:

  Integer. Number of trials for the running-mean window (default: 50).

- palette:

  Character. RColorBrewer palette name (default: "Set1").

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'bouts',
  'segments').

- data_type:

  For SAP objects: Which results to plot ('similarity' or 'maturation').

## Value

A ggplot2 object, printed as a side-effect and returned invisibly.

## Details

**Data source for SAP objects:**

- `data_type = "similarity"`: Uses
  `x$features[[feature_type]]$trajectory_similarity$similarity`. Common
  score columns: "rms_similarity", "correlation", "frechet_similarity",
  "dtw_similarity".

- `data_type = "maturation"`: Uses
  `x$features[[feature_type]]$maturation_scores`. Common score columns:
  "maturation_score", "stability_index".

**Label selection:** When `labels = NULL`, the function selects the
first and last labels from the sorted unique labels in the data
(numeric-aware sorting). This typically corresponds to the earliest and
latest days in a developmental series. Pass an explicit vector to
customise.

**Running mean:** The running mean is computed with
[`stats::filter()`](https://rdrr.io/r/stats/filter.html) using a centred
window of `window_size` trials. It is only drawn when
`running_mean = TRUE` and there are at least `window_size` trials in a
label.

## See also

[`trajectory_similarity`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_similarity.md)
for similarity computation,
[`trajectory_maturation`](https://lxiao06.github.io/ASAP/dev/reference/trajectory_maturation.md)
for maturation score computation

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot RMS similarity for first and last labels
plot_trajectory_trials(sap,
  data_type = "similarity",
  score_col = "rms_similarity"
)

# Plot maturation scores for specific days without running mean
plot_trajectory_trials(sap,
  data_type = "maturation",
  score_col = "maturation_score",
  labels = c("60", "80", "100"),
  running_mean = FALSE
)
} # }
```
