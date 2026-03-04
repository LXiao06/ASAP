# Calculate Amplitude Envelope for Audio Segment

Calculates the amplitude envelope for a specified segment of an audio
file.

## Usage

``` r
amp_env(
  segment_row,
  wav_dir = NULL,
  msmooth = NULL,
  amp_normalize = c("none", "peak", "rms"),
  plot = FALSE
)
```

## Arguments

- segment_row:

  A row from either find_motif output or segment class in SAP object

- wav_dir:

  Optional path to WAV files directory (default: NULL)

- msmooth:

  Numeric vector of length 2 for envelope smoothing:

  - First value: Window length in number of points

  - Second value: Overlap between windows (percentage)

  If NULL, no smoothing is applied

- amp_normalize:

  Waveform amplitude normalization before envelope extraction: one of
  "none", "peak", or "rms" (default: "none")

- plot:

  Logical, whether to plot the envelope (default: FALSE)

## Value

A numeric vector containing the amplitude envelope values

## Details

The function processes audio segments with the following steps:

- Validates input segment data

- Constructs correct file path

- Reads specified portion of audio file

- Calculates amplitude envelope

- Optionally smooths and normalizes the envelope

## See also

[`find_motif`](https://lxiao06.github.io/ASAP/dev/reference/find_motif.md)
for creating segment data

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic envelope calculation
env <- amp_env(segments[1,], wav_dir = "path/to/wavs")

# With smoothing and normalization
env <- amp_env(segments[1,],
               wav_dir = "path/to/wavs",
               msmooth = c(256, 50))

# With waveform peak normalization
env <- amp_env(segments[1,],
               wav_dir = "path/to/wavs",
               msmooth = c(256, 50),
               amp_normalize = "peak")

# With waveform RMS normalization
env <- amp_env(segments[1,],
               wav_dir = "path/to/wavs",
               msmooth = c(256, 50),
               amp_normalize = "rms")

# With plot
env <- amp_env(segments[1,],
               wav_dir = "path/to/wavs",
               msmooth = c(256, 50),
               plot = TRUE)
} # }
```
