# Flag Anomalous Labels in a SAP Object

Assigns a vector of anomalous labels to a SAP object for a given feature
type. These anomalous labels can then be excluded when calling
[`export_feature_csv`](https://lxiao06.github.io/ASAP/dev/reference/export_feature_csv.md).

## Usage

``` r
flag_anomalous_labels(
  x,
  feature_type,
  labels,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  A SAP object

- feature_type:

  Character. The feature type under which to store the anomaly list
  (e.g., `"motif"`, `"syllable"`)

- labels:

  Character vector. The labels to designate as anomalous

- overwrite:

  Logical. If `TRUE`, replaces any existing anomalous labels for this
  feature type. If `FALSE` (default), appends to the existing anomalous
  labels

- verbose:

  Logical. Print progress and summary messages (default: `TRUE`)

## Value

The updated SAP object, invisibly.
