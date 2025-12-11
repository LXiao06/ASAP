# Extract and pad spectrograms from audio segments

Extract and pad spectrograms from audio segments

Default method for extract_spec

## Usage

``` r
extract_spec(x, ...)

# Default S3 method
extract_spec(
  x,
  wav_dir = NULL,
  cores = NULL,
  wl = 256,
  ovlp = 20,
  wn = "hanning",
  freq_range = c(1, 10),
  fftw = TRUE,
  ...
)

# S3 method for class 'Sap'
extract_spec(
  x,
  segment_type = c("segments", "syllables"),
  sample_percent = NULL,
  balanced = FALSE,
  labels = NULL,
  seed = 222,
  cores = NULL,
  wl = 256,
  ovlp = 20,
  wn = "hanning",
  freq_range = c(1, 10),
  fftw = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  Data frame with audio segment information or a SAP object

- ...:

  Additional parameters passed to spectro

- wav_dir:

  Directory containing WAV files (default: NULL)

- cores:

  Number of CPU cores to use for parallel processing (default: NULL,
  uses available cores - 1)

- wl:

  Window length for spectrogram analysis (default: 256 ) Time
  resolution: wl/sampling rate per window Frequency resolution: sampling
  rate/wl Hz per bin

- ovlp:

  Overlap percentage between windows (default: 20 )

- wn:

  Window function name (default: "hanning")

- freq_range:

  Frequency range to analyze in kHz (default: c(1,10) )

- fftw:

  Logical, use FFTW or not (default: TRUE)

- segment_type:

  Type of segments to process: "segments" or "syllables" (default:
  "segments")

- sample_percent:

  Percentage of segments to sample (default: NULL)

- balanced:

  Whether to balance samples across labels (default: FALSE)

- labels:

  Specific labels to include (default: NULL)

- seed:

  Random seed for reproducible sampling (default: 222)

- verbose:

  Whether to display progress messages (default: TRUE)

## Value

For default method: data frame with original metadata and flattened,
padded spectrograms For SAP method: updated SAP object with spectrograms
added to features
