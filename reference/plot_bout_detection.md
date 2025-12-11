# Plot Bout Detection Results

Creates a two-panel plot showing RMS envelope and spectrogram with bout
boundaries.

## Usage

``` r
plot_bout_detection(
  wav_file,
  time_points,
  rms_env,
  rms_threshold,
  bout_df,
  bout_onsets,
  bout_offsets,
  wl,
  ovlp
)
```

## Arguments

- rms_threshold:

  Threshold for bout detection (default: 0.1)

- wl:

  Window length for RMS calculation (default: 1024)

- ovlp:

  Overlap percentage between windows (default: 50)
