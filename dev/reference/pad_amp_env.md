# Pad Amplitude Envelope for Alignment

Internal function to pad and align amplitude envelopes for bout
analysis.

## Usage

``` r
pad_amp_env(
  segment_row,
  wav_dir = NULL,
  msmooth = NULL,
  max_pre_window = NULL,
  max_post_window = NULL,
  samples_per_second = NULL,
  pad_value = 0,
  ...
)
```

## Arguments

- segment_row:

  Single row of segment data

- wav_dir:

  Directory containing WAV files

- msmooth:

  Smoothing parameters

- max_pre_window:

  Maximum pre-window time

- max_post_window:

  Maximum post-window time

- samples_per_second:

  Sampling rate

- pad_value:

  Value to use for padding

- ...:

  Additional arguments

## Value

Padded and aligned amplitude envelope vector
