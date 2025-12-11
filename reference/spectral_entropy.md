# Calculate Spectral Entropy for Audio Segments

Calculates spectral entropy (Wiener or Shannon) for audio segments,
providing measures of spectral uniformity or complexity in sound
signals.

## Usage

``` r
spectral_entropy(x, ...)

# Default S3 method
spectral_entropy(
  x,
  start_time = NULL,
  end_time = NULL,
  wl = 512,
  wn = "hanning",
  ovlp = 50,
  fftw = TRUE,
  freq_range = c(500, 15000),
  threshold = 10,
  method = c("weiner", "shannon"),
  normalize = FALSE,
  plot = TRUE,
  ...
)

# S3 method for class 'data.frame'
spectral_entropy(
  x,
  wav_dir = NULL,
  wl = 512,
  wn = "hanning",
  ovlp = 50,
  fftw = TRUE,
  freq_range = c(500, 15000),
  threshold = 10,
  method = c("weiner", "shannon"),
  normalize = FALSE,
  plot = TRUE,
  plot_entropy_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  cores = NULL,
  ...
)

# S3 method for class 'Sap'
spectral_entropy(
  x,
  segment_type = c("motifs", "syllables", "segments"),
  sample_percent = NULL,
  balanced = FALSE,
  labels = NULL,
  clusters = NULL,
  cores = NULL,
  seed = 222,
  wl = 512,
  wn = "hanning",
  ovlp = 50,
  fftw = TRUE,
  freq_range = c(500, 15000),
  threshold = 10,
  method = c("weiner", "shannon"),
  normalize = FALSE,
  plot = TRUE,
  plot_entropy_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  ordered = FALSE,
  descending = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'matrix'
spectral_entropy(
  x,
  labels = NULL,
  plot_entropy_lim = NULL,
  color_palette = NULL,
  n_colors = 500,
  main = "Spectral Entropy Heatmap",
  ...
)
```

## Arguments

- x:

  Input object:

  - character: path to WAV file (default method)

  - data frame: containing segment information

  - SAP object

  - Pre-computed matrix for spectral entropy

- ...:

  Additional arguments passed to levelplot

- start_time:

  Numeric, start time in seconds (for default method)

- end_time:

  Numeric, end time in seconds (for default method)

- wl:

  Window length for FFT (default: 512)

- wn:

  Window name for FFT (default: "hanning")

- ovlp:

  Overlap percentage between windows (default: 50)

- fftw:

  Logical, use FFTW or not (default: TRUE)

- freq_range:

  Frequency range for analysis c(min, max) in Hz (default: c(500,
  15000))

- threshold:

  Amplitude threshold for power spectrum (default: 10)

- method:

  Entropy type ("weiner" or "shannon")

- normalize:

  Logical, whether to normalize entropy values (default: FALSE)

- plot:

  Logical, whether to plot results (default: TRUE)

- wav_dir:

  Directory containing WAV files (for data frame methods)

- plot_entropy_lim:

  Optional limits for entropy plot c(min, max)

- color_palette:

  Custom color palette function

- n_colors:

  Number of colors in palette (default: 500)

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

Depending on the method used:

- Default method (WAV file):

  - Matrix:

    With columns: time, entropy

- Data frame method:

  - List:

    With components: entropy_matrix, reference_time, original_times,
    plot (if plot=TRUE)

- SAP method:

  - SAP object:

    Updated with entropy features

- Matrix method:

  - Lattice plot:

    (invisibly)

## Details

Wiener entropy measures the uniformity of power distribution:

- Non-normalized (default): returns log-scaled values from -Inf to 0

  - 0 indicates uniform distribution (white noise)

  - Large negative values indicate structured sound (pure tones)

- Normalized: returns values from 0 to 1

  - 1 indicates uniform distribution

  - Values close to 0 indicate structured sound

Shannon entropy measures information content:

- Non-normalized: returns values \>= 0

  - Higher values indicate more uniformity/randomness

  - Lower values indicate more structure/predictability

- Normalized: returns values from 0 to 1

  - Scale relative to maximum possible entropy for vector length

## Examples

``` r
if (FALSE) { # \dontrun{
# Calculate entropy for a single WAV file
entropy <- spectral_entropy("path/to/sound.wav",
                          start_time = 1,
                          end_time = 2,
                          method = "weiner",
                          normalize = FALSE)

# Calculate normalized Shannon entropy for multiple segments using data frame method
entropy_list <- spectral_entropy(segments,
                               wav_dir = "path/to/wavs",
                               method = "shannon",
                               normalize = TRUE,
                               freq_range = c(1000, 10000))

# Method for Sap objects
sap <- spectral_entropy(sap,
                       segment_type = "motifs",
                       method = "weiner",
                       normalize = FALSE,
                       plot = TRUE)

# Access the entropy matrix
wiener_entropy <- sap$features$motif$weiner_entropy

# Method for entropy matrices
# Plot existing entropy matrix with custom settings
spectral_entropy(wiener_entropy,
                plot_entropy_lim = c(-3, 0),
                color_palette = colorRampPalette(c("purple", "blue",
                                                  "cyan", "yellow",
                                                  "orange", "red")),
                main = "Custom Entropy Visualization")

# Plot subset of entropy matrix using labels
spectral_entropy(wiener_entropy,
                labels = c("a", "b", "c"),
                plot_entropy_lim = c(0, 1),
                n_colors = 1000)
} # }
```
