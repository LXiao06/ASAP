# Analyze Spectral Features of Audio Segments

Calculates comprehensive spectral features from audio segments with
support for parallel processing and organized data structures.

## Usage

``` r
analyze_spectral(x, ...)

# Default S3 method
analyze_spectral(
  x,
  wav_dir = NULL,
  cores = NULL,
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

# S3 method for class 'Sap'
analyze_spectral(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  sample_percent = NULL,
  balanced = FALSE,
  labels = NULL,
  seed = 222,
  cores = NULL,
  wl = 512,
  ovlp = 50,
  wn = "hanning",
  fftw = TRUE,
  freq_range = NULL,
  threshold = 15,
  fsmooth = 0.1,
  fast = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to analyze, either a data frame or SAP object

- ...:

  Additional arguments passed to specific methods

- wav_dir:

  Directory containing WAV files

- cores:

  Number of cores for parallel processing

- wl:

  Window length for spectral analysis (default: 512)

- ovlp:

  Overlap percentage (0-100) (default: 50)

- wn:

  Window name ("hanning", "hamming", etc.)

- fftw:

  Logical, use FFTW or not (default: TRUE)

- freq_range:

  Frequency range c(min, max) in kHz

- threshold:

  Threshold for frequency tracking (default: 15)

- fsmooth:

  Frequency smoothing parameter (default: 0.1)

- fast:

  Whether to skip peak frequency calculation (default: TRUE)

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'bouts',
  'segments')

- sample_percent:

  For SAP objects: Percentage of segments to sample

- balanced:

  For SAP objects: Whether to balance groups across labels

- labels:

  For SAP objects: Specific labels to include

- seed:

  For SAP objects: Random seed for sampling (default: 222)

- verbose:

  For SAP objects: Whether to print progress messages

## Value

For default method: A data frame containing spectral features for all
segments For SAP objects: Updated SAP object with spectral features
stored in features slot

## Details

For data frames:

- Requires columns: filename, start_time, end_time

- Processes segments in parallel

- Calculates comprehensive spectral features

- Returns combined results

For SAP objects:

- Supports multiple segment types

- Optional balanced sampling

- Stores results in features slot

- Preserves segment metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze segments from data frame
features <- analyze_spectral(segments,
                            wav_dir = "path/to/wavs",
                            cores = 4,
                            freq_range = c(1, 10))

# Basic analysis from SAP object
sap_obj <- analyze_spectral(sap_object,
                           segment_type = "motifs",
                           cores = 4)

# Balanced sampling with specific labels
sap_obj <- analyze_spectral(sap_object,
                           segment_type = "syllables",
                           sample_percent = 80,
                           balanced = TRUE,
                           labels = c("a", "b"))

# Custom spectral parameters
sap_obj <- analyze_spectral(sap_object,
                           segment_type = "motifs",
                           wl = 1024,
                           ovlp = 75,
                           freq_range = c(2, 8))
} # }

```
