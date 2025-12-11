# Plot Syllable Detection Results

Internal function to create visualization of syllable detection results.

## Usage

``` r
plot_syllable_detection(
  times,
  frequencies,
  sp,
  syllables,
  max_level_db,
  ref_level_db,
  final_threshold,
  final_envelope,
  silence_threshold,
  smooth = FALSE
)
```

## Arguments

- times:

  Time points vector

- frequencies:

  Frequency points vector

- sp:

  Spectrogram matrix

- syllables:

  Detected syllables data frame

- max_level_db:

  Maximum threshold level

- ref_level_db:

  Reference level

- final_threshold:

  Final detection threshold

- final_envelope:

  Final detection envelope

- silence_threshold:

  Silence threshold

- smooth:

  Whether to smooth visualization

## Value

Function that creates the plot

## Details

Creates visualization with:

- Spectrogram with syllable boundaries

- Detection envelope trace

- Threshold indicators

- Optional smoothing
