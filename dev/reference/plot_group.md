# Plot Group of Segments

An internal helper function to plot spectrograms for a group of
segments.

## Usage

``` r
plot_group(
  segments,
  base_path,
  n_samples,
  fft_window_size,
  overlap,
  dark_mode,
  legend,
  by_column,
  label,
  cluster = NULL,
  metrics_labels = NULL,
  metrics_label_cex = 0.9,
  metrics_label_font = 2,
  metrics_label_family = "mono",
  metrics_label_location = c("overlay", "margin")
)
```

## Arguments

- segments:

  A data frame of segment information

- base_path:

  Base directory path for audio files

- n_samples:

  Number of samples to plot

- fft_window_size:

  Size of FFT window

- overlap:

  Overlap between windows

- dark_mode:

  Use dark theme

- legend:

  Show spectrogram legend

- by_column:

  Arrange plots by columns

- label:

  Label for the group of segments

- cluster:

  Optional cluster identifier

- metrics_labels:

  Named character vector of formatted metric labels keyed by segment row
  name, or NULL (default) to omit

- metrics_label_cex:

  Numeric character expansion for metric labels

- metrics_label_font:

  Font face for metric labels

- metrics_label_family:

  Font family for metric labels

- metrics_label_location:

  Location for metric labels ("overlay" or "margin")

## Value

Generates a plot of segment spectrograms as a side effect
