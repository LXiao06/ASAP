# Create Spectrogram Matrices for Song Trajectory Analysis

Creates spectrogram matrices from sliding windows for song trajectory
analysis.

## Usage

``` r
create_trajectory_matrix(x, ...)

# Default S3 method
create_trajectory_matrix(
  x,
  wav_dir,
  window_size = 0.1,
  step_size = 0.005,
  wl = 128,
  ovlp = 50,
  fftw = TRUE,
  flim = c(1, 12),
  cores = NULL,
  ...
)

# S3 method for class 'Sap'
create_trajectory_matrix(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  data_type = NULL,
  clusters = NULL,
  sample_percent = NULL,
  balanced = FALSE,
  labels = NULL,
  seed = 222,
  window_size = 0.1,
  step_size = 0.005,
  wl = 128,
  ovlp = 50,
  fftw = TRUE,
  flim = c(1, 12),
  cores = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to process, either a data frame or SAP object

- ...:

  Additional arguments passed to specific methods

- wav_dir:

  Directory containing WAV files (for default method)

- window_size:

  Size of sliding window in seconds (default: 0.1)

- step_size:

  Step size between windows (default: 0.005)

- wl:

  Window length for spectrogram (default: 128)

- ovlp:

  Overlap percentage (default: 50)

- fftw:

  Logical, use FFTW or not (default: TRUE)

- flim:

  Frequency limits (default: c(1, 12))

- cores:

  Number of processing cores

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'bouts',
  'segments')

- data_type:

  For SAP objects: Type of data to analyze

- clusters:

  For SAP objects: Specific clusters to include

- sample_percent:

  For SAP objects: Percentage to sample

- balanced:

  For SAP objects: Whether to balance across groups

- labels:

  For SAP objects: Specific labels to include

- seed:

  For SAP objects: Random seed

- verbose:

  For SAP objects: Whether to print progress

## Value

For default method: A list containing:

- spectrogram_matrix: Matrix of spectrogram vectors

- sliding_windows: Data frame of window information

For SAP objects: Updated SAP object with trajectory matrix stored in
features slot

## Details

Creates trajectory matrix with the following steps:

- Generates sliding windows for each segment

- Computes spectrograms for each window

- Combines results into matrix form

For SAP objects, additional features include:

- Support for different segment types

- Optional cluster/label filtering

- Balanced sampling options

- Results storage in features slot

## Examples

``` r
if (FALSE) { # \dontrun{
# Create trajectory matrix from segments
matrix <- create_trajectory_matrix(segments,
                                  wav_dir = "path/to/wavs",
                                  window_size = 0.1,
                                  step_size = 0.005)

# Create trajectory matrix from SAP object
sap_obj <- create_trajectory_matrix(sap_object,
                                   segment_type = "motifs",
                                   balanced = TRUE,
                                   sample_percent = 80)

# Create matrix with specific clusters
sap_obj <- create_trajectory_matrix(sap_object,
                                   segment_type = "syllables",
                                   clusters = c(1, 2),
                                   labels = c("a", "b"))
} # }
```
