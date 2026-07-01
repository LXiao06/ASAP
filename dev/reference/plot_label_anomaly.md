# Plot Label Anomaly Detection Results

Visualises label anomaly detection results, highlighting time points
flagged as having abnormal variability or maturation patterns. When
`heatmap_only = FALSE` (default), a multi-panel figure of individual
metric traces is returned (heatmap excluded). When
`heatmap_only = TRUE`, a standalone summary heatmap is returned instead.

## Usage

``` r
plot_label_anomaly(
  detection_result,
  sap = NULL,
  highlight_color = "red",
  heatmap_only = FALSE
)
```

## Arguments

- detection_result:

  Output from
  [`detect_anomalous_labels`](https://lxiao06.github.io/ASAP/dev/reference/detect_anomalous_labels.md).

- sap:

  Original SAP object (reserved for future enhancements).

- highlight_color:

  Color for anomalous labels (default: `"red"`).

- heatmap_only:

  Logical. If `TRUE`, returns only the flag heatmap as a standalone
  summary figure. If `FALSE` (default), returns a multi-panel figure of
  individual metric traces without the heatmap.

## Value

A `ggplot` or `patchwork` object.

## See also

[`detect_anomalous_labels`](https://lxiao06.github.io/ASAP/dev/reference/detect_anomalous_labels.md)

## Examples

``` r
if (FALSE) { # \dontrun{
anomaly_results <- detect_anomalous_labels(sap, segment_type = "motifs")

# Multi-panel metric traces
plot_label_anomaly(anomaly_results, sap)

# Standalone heatmap summary
plot_label_anomaly(anomaly_results, sap, heatmap_only = TRUE)
} # }
```
