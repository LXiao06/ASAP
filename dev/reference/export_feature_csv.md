# Export a Stored Feature Data Frame to CSV

Extracts a pre-computed feature data frame stored in
`x$features[[feature_type]][[slot_name]]` and writes it to a CSV file.
Works universally with any feature slot produced by ASAP analysis
functions, including spectral features, trajectory
similarity/variability/maturation scores, and any other data frame or
list-of-data-frames stored in the `features` slot of a SAP object.

When the stored object is a plain data frame (e.g., `spectral_feature`,
`maturation_scores`), it is exported directly.

When the stored object is a list (e.g., `trajectory_similarity`,
`trajectory_dispersion`, `trajectory_path_deviation`), use the
`sub_table` argument to select which element of the list to export
(e.g., `sub_table = "similarity"`, `"dispersion"`, `"deviation"`).

## Usage

``` r
export_feature_csv(
  x,
  feature_type,
  slot_name,
  sub_table = NULL,
  merge_segments = TRUE,
  segment_type = NULL,
  time_match_digits = 3,
  balance_labels = FALSE,
  balance_threshold = 0.5,
  exclude_anomalies = TRUE,
  csv_filename = NULL,
  output_dir = ".",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  A SAP object containing the features to export

- feature_type:

  Character. The top-level feature key inside `x$features`, e.g.
  `"motif"`, `"syllable"`, `"bout"`, or `"segment"`

- slot_name:

  Character. The sub-slot key inside `x$features[[feature_type]]`, e.g.
  `"spectral_feature"`, `"maturation_scores"`,
  `"trajectory_similarity"`, `"trajectory_dispersion"`, or
  `"trajectory_path_deviation"`

- sub_table:

  Character or `NULL`. When the stored object is a list (trajectory
  results), the name of the list element to extract as a data frame
  (default: `NULL`)

- merge_segments:

  Logical. If `TRUE` (default) and the extracted data frame contains
  `filename`, `start_time`, and `end_time` columns, it is left-joined
  against the canonical segment metadata stored in `x[[segment_type]]`
  to enrich the export with per-segment columns. When those columns are
  absent, the merge is skipped

- segment_type:

  Character or `NULL`. The segment slot to use for metadata merging (one
  of `"motifs"`, `"syllables"`, `"bouts"`, `"segments"`). If `NULL`
  (default), derived automatically from `feature_type` by appending
  `"s"`

- time_match_digits:

  Integer. Number of decimal places used when building the merge key
  (default: `3`)

- balance_labels:

  Logical. If `TRUE`, labels whose row count falls below threshold are
  dropped before writing the CSV (default: `FALSE`)

- balance_threshold:

  Numeric in `(0, 1]`. Fraction of the median label count used as the
  minimum retention threshold when `balance_labels = TRUE` (default:
  `0.5`)

- exclude_anomalies:

  Logical. If `TRUE`, excludes rows belonging to anomalous labels stored
  under `x$features[[feature_type]][["anomalous_labels"]]` if available
  (default: `TRUE`)

- csv_filename:

  Character. Output file name (required)

- output_dir:

  Character. Directory where the CSV is written (default: `"."`)

- overwrite:

  Logical. Whether to overwrite an existing file (default: `FALSE`)

- verbose:

  Logical. Print progress and summary messages (default: `TRUE`)

## Value

The exported data frame, invisibly.

## Details

When the stored object is a list, the following sub-tables are typical
choices for `sub_table`:

- `trajectory_similarity`: `"similarity"`, `"summary"`

- `trajectory_dispersion`: `"dispersion"`, `"pairwise"`,
  `"path_length"`, `"summary"`

- `trajectory_path_deviation`: `"deviation"`, `"deviation_stats"`,
  `"summary"`

## Examples

``` r
if (FALSE) { # \dontrun{
# Export spectral features stored from a previous analyze_spectral() run
export_feature_csv(sap,
  feature_type = "motif",
  slot_name    = "spectral_feature",
  csv_filename = "spectral_features.csv"
)

# Export trajectory similarity scores
export_feature_csv(sap,
  feature_type   = "motif",
  slot_name      = "trajectory_similarity",
  sub_table      = "similarity",
  merge_segments = FALSE,
  csv_filename   = "traj_similarity.csv"
)
} # }
```
