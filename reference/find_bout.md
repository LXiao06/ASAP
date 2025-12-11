# Detect Song Bouts in Audio Recordings

Detects and analyzes song bouts using RMS amplitude thresholding with
bandpass filtering.

## Usage

``` r
find_bout(x, ...)

# Default S3 method
find_bout(
  x,
  wl = 1024,
  ovlp = 50,
  norm_method = c("quantile", "max"),
  rms_threshold = 0.1,
  min_duration = 0.5,
  gap_duration = 0.3,
  edge_window = 0.05,
  freq_range = c(3, 5),
  plot = TRUE,
  save_plot = FALSE,
  plot_dir = NULL,
  ...
)

# S3 method for class 'Sap'
find_bout(
  x,
  day = NULL,
  indices = NULL,
  segment_type = "motifs",
  cores = NULL,
  save_plot = FALSE,
  plot_percent = 10,
  wl = 1024,
  ovlp = 50,
  norm_method = c("quantile", "max"),
  rms_threshold = 0.1,
  min_duration = 0.5,
  gap_duration = 0.3,
  edge_window = 0.05,
  freq_range = c(3, 5),
  summary = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to analyze, either a file path or SAP object

- ...:

  Additional arguments passed to specific methods

- wl:

  Window length for RMS calculation (default: 1024)

- ovlp:

  Overlap percentage between windows (default: 50)

- norm_method:

  Method for normalizing RMS values ("quantile" or "max")

- rms_threshold:

  Threshold for bout detection (default: 0.1)

- min_duration:

  Minimum bout duration in seconds (default: 0.5)

- gap_duration:

  Minimum gap between bouts (default: 0.3)

- edge_window:

  Time window for edge effects (default: 0.05)

- freq_range:

  Frequency range for bandpass filter (default: c(3, 5))

- plot:

  Whether to display visualization (default: TRUE)

- save_plot:

  Whether to save plots to file (default: FALSE)

- plot_dir:

  Directory for saving plots

- day:

  For SAP objects: Days to process

- indices:

  For SAP objects: Specific indices to process

- segment_type:

  For SAP objects: Type of segments (default: "motifs")

- cores:

  For SAP objects: Number of processing cores

- plot_percent:

  For SAP objects: Percentage of files to plot (default: 10)

- summary:

  For SAP objects: Include additional statistics (default: FALSE)

- verbose:

  For SAP objects: Print progress messages (default: TRUE)

## Value

For default method: Data frame containing:

- filename: Name of WAV file

- selec: Bout number

- start_time: Onset time

- end_time: Offset time

For SAP objects: Updated object with bout information in bouts slot

## Details

For WAV files:

- Applies bandpass filtering to focus on vocalization frequencies

- Calculates RMS envelope with specified window parameters

- Detects bouts using adaptive thresholding

- Handles edge cases and minimum duration constraints

- Creates optional visualizations

For SAP objects:

- Processes multiple recordings in parallel

- Validates bouts against existing motif detections

- Provides optional summary statistics

- Maintains metadata relationships

- Supports selective plotting

When summary = TRUE for SAP objects with motifs:

- n_motifs: Count of motifs per bout

- align_time: First motif time for alignment

- bout_number_day: Sequential numbering

- bout_gap: Time from previous bout

## See also

[`segment`](https://lxiao06.github.io/ASAP/reference/segment.md) for
syllable-level segmentation

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic bout detection from file
bouts <- find_bout("song.wav",
                   rms_threshold = 0.1,
                   min_duration = 0.7)

# Custom parameters with visualization
bouts <- find_bout("song.wav",
                   freq_range = c(2, 8),
                   plot = TRUE,
                   save_plot = TRUE)

# Process SAP object with summary
sap_obj <- find_bout(sap_object,
                     segment_type = "motifs",
                     day = c(30, 40),
                     summary = TRUE)

# Process specific files with plots
sap_obj <- find_bout(sap_object,
                     indices = 1:5,
                     save_plot = TRUE,
                     cores = 4)
} # }
```
