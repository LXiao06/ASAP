# Visualize Song Data

Creates spectrograms and visualizes acoustic data from WAV files or SAP
objects.

## Usage

``` r
visualize_song(x, ...)

# Default S3 method
visualize_song(
  x,
  start_time_in_second = NULL,
  end_time_in_second = NULL,
  fft_window_size = 1024,
  overlap = 0.5,
  dark_mode = TRUE,
  legend = FALSE,
  keep.par = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
visualize_song(
  x,
  template_clips = FALSE,
  indices = NULL,
  n_samples = NULL,
  random = TRUE,
  start_time_in_second = NULL,
  end_time_in_second = NULL,
  fft_window_size = 1024,
  overlap = 0.75,
  keep.par = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  An object to visualize, either a file path or SAP object

- ...:

  Additional arguments passed to specific methods

- start_time_in_second:

  Numeric start time in seconds

- end_time_in_second:

  Numeric end time in seconds

- fft_window_size:

  Size of FFT window (default: 512 for default, 1024 for SAP)

- overlap:

  Overlap between windows (default: 0.5 for default, 0.75 for SAP)

- dark_mode:

  For default method: Use dark theme (default: TRUE)

- legend:

  For default method: Show spectrogram legend (default: FALSE)

- keep.par:

  Preserve plotting parameters

- verbose:

  Print processing messages

- template_clips:

  Logical. For SAP objects: whether to visualize original songs (FALSE)
  or template clips (TRUE) (default: FALSE)

- indices:

  For SAP objects: Numeric vector of specific indices to visualize

- n_samples:

  For SAP objects: Number of samples to visualize if indices is NULL.
  Default is 6 or max available

- random:

  For SAP objects: Randomly sample songs if TRUE

## Value

Generates spectrogram plot(s) and returns the input object invisibly.

## Details

For WAV files:

- Creates single spectrogram using FFmpeg's FFT

- Customizable time range and FFT settings

- Optional dark mode and legend

For SAP objects:

- Creates multi-panel spectrograms

- Supports random or sequential sampling

- Maintains plotting state for sequential viewing

- Adds day and label information to plots

## Examples

``` r
if (FALSE) { # \dontrun{
# Visualize a single WAV file
visualize_song("path/to/song.wav",
               start_time_in_second = 10,
               end_time_in_second = 20)

# Basic visualization from SAP object
visualize_song(sap_object, n_sample = 4)

# Visualize specific indices with custom FFT settings
visualize_song(sap_object,
               indices = c(1, 3, 5),
               fft_window_size = 2048,
               overlap = 0.8)

# Sequential visualization with time ranges
visualize_song(sap_object,
               n_sample = 6,
               random = FALSE,
               start_time_in_second = rep(0, 6),
               end_time_in_second = rep(5, 6))
} # }
```
