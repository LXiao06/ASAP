# Internal Spectral Analysis Function

Internal function to perform spectral analysis on a single audio
segment.

Internal functions adapted from the warbleR package for spectral
analysis. These functions are modified versions of the original warbleR
code.

## Usage

``` r
spectral_analysis(
  x,
  wav_dir = NULL,
  wl = 512,
  ovlp = 50,
  wn = "hanning",
  fftw = TRUE,
  freq_range = NULL,
  threshold = 15,
  fsmooth = 0.1,
  fast = TRUE,
  ...
)
```

## Arguments

- x:

  Single row data frame with segment information

- wav_dir:

  Path to WAV files directory

- wl:

  Window length for analysis

- ovlp:

  Overlap percentage

- wn:

  Window name

- freq_range:

  Frequency range

- threshold:

  Detection threshold

- fsmooth:

  Smoothing parameter

- fast:

  Skip peak frequency calculation

- ...:

  Additional arguments

## Value

Data frame with spectral features

## Details

Original source: warbleR package Citation: Araya-Salas, M. and
Smith-Vidaurre, G. (2017), warbleR: an r package to streamline analysis
of animal acoustic signals. Methods Ecol Evol. 8, 184-191.
