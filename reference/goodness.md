# Pitch Goodness Analysis

Calculates and visualizes pitch goodness for audio segments using
cepstral analysis. Supports single WAV file analysis, batch processing,
and multiple visualization options.

## Usage

``` r
goodness(x, ...)

# Default S3 method
goodness(
  x,
  start_time = NULL,
  end_time = NULL,
  wl = 512,
  ovlp = 50,
  fmax = 1500,
  plot = TRUE,
  ...
)

# S3 method for class 'data.frame'
goodness(
  x,
  wav_dir = NULL,
  wl = 512,
  ovlp = 50,
  fmax = 1500,
  plot = TRUE,
  plot_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  cores = NULL,
  ...
)

# S3 method for class 'Sap'
goodness(
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
  fmax = 1500,
  plot = TRUE,
  plot_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  ordered = FALSE,
  descending = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'matrix'
goodness(
  x,
  labels = NULL,
  plot_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  main = "Pitch Goodness Heatmap",
  ...
)
```

## Arguments

- x:

  Input object:

  - character: path to WAV file (default method)

  - data frame: containing segment information

  - SAP object

  - Pre-computed goodness matrix

- ...:

  Additional arguments passed to methods or plotting functions

- start_time:

  Numeric, start time in seconds (for default method)

- end_time:

  Numeric, end time in seconds (for default method)

- wl:

  Window length for spectral analysis (default: 512)

- ovlp:

  Overlap percentage between windows (default: 50)

- fmax:

  Maximum frequency to consider (default: 1500 Hz)

- plot:

  Logical, whether to generate visualization (default: TRUE)

- wav_dir:

  Directory containing WAV files (for data frame methods)

- plot_lim:

  Optional vector of length 2 specifying goodness limits for plotting

- color_palette:

  Function generating color palette (default: black-to-white spectrum)

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

  Logical, whether to print progress messages(default: TRUE)

- main:

  Plot title (for matrix method)

## Value

Depending on the method used:

- Default method (WAV file):

  - Matrix:

    With columns: time, goodness

- Data frame method:

  - List:

    With components: goodness_matrix, reference_time, original_times,
    plot (if plot=TRUE)

- SAP method:

  - SAP object:

    Updated with goodness features

- Matrix method:

  - Lattice plot:

    (invisibly)

## Details

The function provides different methods depending on the input type:

Default method (WAV file):

- Analyzes a single WAV file within specified time window

- Uses cepstral analysis for pitch goodness calculation

- Creates spectrogram with goodness overlay visualization

Data frame method:

- Processes multiple segments in parallel

- Supports batch processing

- Provides heatmap visualization for multiple segments

SAP object method:

- Integrates with SAP feature analysis pipeline

- Supports segment filtering and ordering

- Stores results in features component

Matrix method:

- Visualizes pre-computed F0 data

- Supports label-based organization

- Adds visual separators between groups

SAP object method:

- Integrates with SAP feature analysis pipeline

- Supports segment filtering and ordering

- Stores results in features component

Matrix method:

- Visualizes pre-computed pitch goodness data

- Supports label-based organization

- Adds visual separators between groups

The pitch goodness measure:

- Based on cepstral peak prominence

- Higher values indicate stronger pitch periodicity

- Computed across time windows

- Normalized for comparison across segments

## See also

[`FF`](https://lxiao06.github.io/ASAP/reference/Fundamental_Frequency.md)
for fundamental frequency analysis
[`spectro`](https://rdrr.io/pkg/seewave/man/spectro.html) for underlying
spectral analysis

## Examples

``` r
if (FALSE) { # \dontrun{
# Single WAV file analysis
goodness("path/to/sound.wav", start_time = 1, end_time = 2)

# Multiple segment analysis
goodness(segments_df, wav_dir = "path/to/wavs",
         plot = TRUE, plot_lim = c(0, 1))

# SAP object analysis
goodness(sap_obj, ordered = TRUE, balanced = TRUE,
   sample_percent = 80, seed = 123)

# Direct matrix visualization
goodness(goodness_matrix, labels = c("a", "b", "c"),
   main = "Custom Title")
} # }
```
