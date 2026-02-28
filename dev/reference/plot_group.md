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
  cluster = NULL
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

## Value

Generates a plot of segment spectrograms as a side effect
