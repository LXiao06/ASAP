# Visualize Song Segments

Creates multi-panel spectrogram visualizations of audio segments from
various sources.

## Usage

``` r
visualize_segments(x, ...)

# Default S3 method
visualize_segments(
  x,
  wav_dir,
  n_samples = NULL,
  seed = NULL,
  fft_window_size = 1024,
  overlap = 0.75,
  dark_mode = TRUE,
  legend = FALSE,
  ...
)

# S3 method for class 'Sap'
visualize_segments(
  x,
  segment_type = c("motifs", "bouts", "segments"),
  labels = NULL,
  clusters = NULL,
  n_samples = NULL,
  by_column = TRUE,
  seed = NULL,
  fft_window_size = 1024,
  overlap = 0.75,
  dark_mode = TRUE,
  legend = FALSE,
  show_maturation = FALSE,
  label_metrics = c("correlation", "maturation_score"),
  maturation_label_cex = 0.9,
  maturation_label_font = 2,
  maturation_label_family = "mono",
  maturation_label_location = c("overlay", "margin"),
  ...
)
```

## Arguments

- x:

  An object to visualize, either a data frame or SAP object

- ...:

  Additional arguments passed to specific methods

- wav_dir:

  For default method: Directory containing WAV files

- n_samples:

  Number of samples to display

- seed:

  Random seed for sample selection

- fft_window_size:

  Size of FFT window (default: 1024)

- overlap:

  Overlap between windows (default: 0.75)

- dark_mode:

  Use dark theme (default: TRUE)

- legend:

  Show spectrogram legend (default: FALSE)

- segment_type:

  For SAP objects: Type of segments ('motifs', 'bouts', 'segments')

- labels:

  For SAP objects: Labels to include

- clusters:

  For SAP objects: Specific clusters to visualize

- by_column:

  For SAP objects: Arrange by columns (default: TRUE)

- show_maturation:

  Logical. If TRUE and segment_type == "motifs", only show segments that
  have maturation scores (default: FALSE)

- label_metrics:

  Character vector of metric column names from the maturation scores to
  display at the top-right of each spectrogram panel (default:
  c("correlation", "maturation_score"))

- maturation_label_cex:

  Numeric character expansion for maturation metric labels (default:
  0.9)

- maturation_label_font:

  Font face for maturation metric labels, passed to base R text
  functions (default: 2, bold)

- maturation_label_family:

  Font family for maturation metric labels (default: "mono")

- maturation_label_location:

  Location for maturation metric labels. Use "overlay" to place labels
  inside each spectrogram with a contrast box, or "margin" to place
  labels above each panel (default: "overlay")

## Value

Generates a multi-panel spectrogram plot and returns the input object
invisibly.

## Details

For data frames:

- Requires columns: filename, start_time, end_time

- Optional column: day_post_hatch for hierarchical file structure

- Supports random sampling of segments

For SAP objects:

- Supports visualization by labels and/or clusters

- Flexible sampling within groups

- Customizable layout (by row or column)

- Automatic handling of file paths

## See also

[`visualize_song`](https://lxiao06.github.io/ASAP/dev/reference/visualize_song.md)
for single file visualization

## Examples

``` r
if (FALSE) { # \dontrun{
# Visualize from data frame
song_df <- data.frame(
  filename = c("song1.wav", "song2.wav"),
  start_time = c(0, 10),
  end_time = c(30, 40)
)
visualize_segments(song_df,
                   wav_dir = "path/to/wav/files",
                   n_samples = 5)

# Basic SAP object visualization
visualize_segments(sap_object,
                   segment_type = "motifs",
                   n_samples = 3)

# Cluster-specific visualization
visualize_segments(sap_object,
                   segment_type = "segments",
                   clusters = c(1, 2),
                   labels = c("a", "b"),
                   n_samples = 4)

# Custom layout
visualize_segments(sap_object,
                   segment_type = "motifs",
                   by_column = FALSE,
                   fft_window_size = 2048)
} # }
```
