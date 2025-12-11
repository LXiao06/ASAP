# Fundamental Frequency Analysis and Visualization

Creates fundamental frequency analyses and visualizations from audio
segments, supporting multiple data types (WAV file, data frames, SAP
objects) and visualization options. Provides both analytical results and
optional single trial or heatmap visualizations.

## Usage

``` r
FF(x, ...)

# Default S3 method
FF(
  x,
  start_time = NULL,
  end_time = NULL,
  wl = 512,
  ovlp = 50,
  fmax = 1400,
  threshold = 10,
  method = c("cepstrum", "yin"),
  plot = TRUE,
  ...
)

# S3 method for class 'data.frame'
FF(
  x,
  wav_dir = NULL,
  wl = 512,
  ovlp = 50,
  fmax = 1400,
  threshold = 10,
  method = c("cepstrum", "yin"),
  plot = TRUE,
  plot_freq_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  cores = NULL,
  ...
)

# S3 method for class 'Sap'
FF(
  x,
  segment_type = c("motifs", "syllables", "segments"),
  sample_percent = NULL,
  balanced = FALSE,
  labels = NULL,
  clusters = NULL,
  cores = NULL,
  seed = 222,
  wl = 512,
  ovlp = 50,
  fmax = 1400,
  threshold = 10,
  method = c("cepstrum", "yin"),
  plot = TRUE,
  plot_freq_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  ordered = FALSE,
  descending = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'matrix'
FF(
  x,
  labels = NULL,
  plot_freq_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  main = "Fundamental Frequency Heatmap",
  ...
)
```

## Arguments

- x:

  Input object:

  - character: path to WAV file (default method)

  - data frame: containing segment information

  - SAP object

  - Pre-computed F0 matrix

- ...:

  Additional arguments passed to methods or lattice::levelplot

- start_time:

  Numeric, start time in seconds (for default method)

- end_time:

  Numeric, end time in seconds (for default method)

- wl:

  Window length for spectral analysis (default: 512)

- ovlp:

  Overlap percentage between windows (default: 50)

- fmax:

  Maximum frequency to consider (default: 1400 Hz)

- threshold:

  Amplitude threshold for pitch detection in % (default: 10)

- method:

  Pitch estimation method ("cepstrum" or "yin")

- plot:

  Logical, whether to generate visualization(default: TRUE)

- wav_dir:

  Directory containing WAV files (for data frame methods)

- plot_freq_lim:

  Optional vector of length 2 specifying frequency limits for plotting

- color_palette:

  Function generating color palette (default: black-to-red spectrum)

- n_colors:

  Number of colors in heatmap (default: 500)

- cores:

  Number of cores for parallel processing

- segment_type:

  For SAP objects: Type of segments (currently only 'motifs')

- sample_percent:

  For SAP objects: Percentage of segments to sample

- balanced:

  For SAP objects: Whether to balance samples across labels

- labels:

  Optional vector of labels for subsetting or grouping

- clusters:

  For SAP objects: Numeric vector of cluster IDs to include

- seed:

  Random seed for sampling (default: 222)

- ordered:

  For SAP objects: Whether to order by feature embeddings

- descending:

  For SAP objects: Direction of embedding-based ordering

- verbose:

  Logical, whether to print progress messages

- main:

  Plot title (for matrix method)

## Value

Returns an object depending on the method used:

- Default method: Matrix of time and frequency values

- Data frame method: List with F0 matrix and metadata

- SAP method: Updated SAP object with F0 features

- Matrix method: Lattice plot object (invisibly)

## Details

The function provides different methods depending on the input type:

Default method (WAV file method):

- Analyzes a single WAV file within specified time window

- Supports both cepstrum and YIN-based pitch detection

- Creates spectrogram with F0 overlay visualization

Data frame method:

- Processes multiple segments in parallel

- Normalizes time series across renditions

- Creates aligned F0 matrix

- Generates heatmap visualization for multiple segments

SAP object method:

- Integrates with SAP feature analysis pipeline

- Supports segment filtering and ordering

- Stores results in features component

Matrix method:

- Visualizes pre-computed F0 data

- Supports label-based organization

- Adds visual separators between groups

## See also

[`levelplot`](https://rdrr.io/pkg/lattice/man/levelplot.html) for
underlying plotting function
[`fund`](https://rdrr.io/pkg/seewave/man/fund.html) for cepstral
analysis YIN algorithm implementation based on the approach used in the
Python librosa package
(<https://librosa.org/doc/main/generated/librosa.yin.html>)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single WAV file analysis
FF("path/to/sound.wav", start_time = 1, end_time = 2)

# Multiple segment analysis
FF(segments_df, wav_dir = "path/to/wavs",
   plot = TRUE, plot_freq_lim = c(0.5, 1.5))

# SAP object analysis
FF(sap_obj, ordered = TRUE, balanced = TRUE,
   sample_percent = 80, seed = 123)

# Direct matrix visualization
FF(f0_matrix, labels = c("a", "b", "c"),
   main = "Custom Title")
} # }
```
