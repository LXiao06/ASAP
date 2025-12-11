# Internal function to calculate the fundamental frequency of time series data derived from a single audio segment

Estimates the fundamental frequency (F0) of an audio segment using
either cepstral analysis or YIN method for a specified audio segment.

## Usage

``` r
FF_single_row(
  segment_row,
  wav_dir = NULL,
  method = c("cepstrum", "yin"),
  wl = 512,
  ovlp = 50,
  fmax = 1400,
  threshold = 10,
  plot = FALSE,
  ...
)
```

## Arguments

- segment_row:

  A single-row data frame containing segment information

- wav_dir:

  Directory containing wav files.

- method:

  Pitch estimation method. Either "cepstrum" or "yin" (default:
  "cepstrum")

- wl:

  Window length for spectral analysis (default: 512)

- ovlp:

  Overlap percentage between windows (default: 80)

- fmax:

  Maximum frequency to consider (default: 1400 Hz)

- threshold:

  Threshold Amplitude threshold for cepstral method in % (default = 10).

- plot:

  Logical, whether to plot the pitch estimation (default: FALSE)

## Value

A matrix with two columns:

- First column: Time points

- Second column: Fundamental frequency values
