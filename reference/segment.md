# Segment Audio into Syllables

Segments audio recordings into syllables using dynamic thresholding of
spectrograms.

## Usage

``` r
segment(x, ...)

# Default S3 method
segment(
  x,
  start_time = NULL,
  end_time = NULL,
  wl = 256,
  ovlp = 80,
  fftw = TRUE,
  flim = c(1, 10),
  silence_threshold = 0.05,
  min_syllable_ms = 50,
  max_syllable_ms = 200,
  min_level_db = 0,
  max_level_db = 40,
  db_delta = 5,
  search_direction = c("up", "down"),
  verbose = TRUE,
  plot = TRUE,
  smooth = FALSE,
  save_plot = FALSE,
  plot_dir = NULL,
  ...
)

# S3 method for class 'Sap'
segment(
  x,
  day = NULL,
  indices = NULL,
  segment_type = c("bouts", "motifs"),
  cores = NULL,
  save_plot = FALSE,
  plot_percent = 10,
  wl = 256,
  ovlp = 80,
  fftw = TRUE,
  flim = c(1, 10),
  silence_threshold = 0.05,
  min_syllable_ms = 50,
  max_syllable_ms = 200,
  min_level_db = 0,
  max_level_db = 40,
  db_delta = 5,
  search_direction = c("up", "down"),
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to segment, either a file path or SAP object

- ...:

  Additional arguments passed to specific methods

- start_time:

  Start time in seconds

- end_time:

  End time in seconds

- wl:

  Window length for spectrogram (default: 256)

- ovlp:

  Overlap percentage (0-100) (default: 80)

- fftw:

  Logical, use FFTW or not (default: TRUE)

- flim:

  Frequency limits in kHz (default: c(1, 10))

- silence_threshold:

  Threshold for silence detection (0-1) (default: 0.05)

- min_syllable_ms:

  Minimum syllable length in milliseconds (default: 50)

- max_syllable_ms:

  Maximum syllable length in milliseconds (default: 200)

- min_level_db:

  Minimum threshold level in dB (Default is 0)

- max_level_db:

  Maximum threshold level in dB (Default is 40)

- db_delta:

  Step size for threshold search in dB (default: 5)

- search_direction:

  Direction for threshold search: "up" or "down".

  - "up": Starts from min_level_db and increases (recommended for quiet
    or variable recordings)

  - "down": Starts from max_level_db and decreases (recommended for
    loud, clear recordings)

- verbose:

  Print progress messages (default: TRUE)

- plot:

  Display detection plot (default: TRUE)

- smooth:

  Smooth spectrogram visualization (default: FALSE)

- save_plot:

  Save detection plot (default: FALSE)

- plot_dir:

  Directory to save plots

- day:

  For SAP objects: Days to process

- indices:

  For SAP objects: Specific indices to process

- segment_type:

  For SAP objects: Type of segments ('bouts', 'motifs')

- cores:

  For SAP objects: Number of processing cores

- plot_percent:

  For SAP objects: Percentage of files to plot (default: 10)

## Value

Returns a data frame containing syllable information:

- filename: Name of the audio file

- selec: Sequential number for each syllable

- threshold: Final threshold used for detection

- .start: Start time relative to the analyzed segment

- .end: End time relative to the analyzed segment

- start_time: Absolute start time in the original audio file

- end_time: Absolute end time in the original audio file

- duration: Duration of syllable in seconds

- silence_gap: Gap to next syllable in seconds (NA for last syllable)

If no syllables are detected, returns NULL.

For SAP objects: Updated object with syllable information in segments
slot

## Details

For WAV files:

- Reads and validates the audio file

- Computes spectrogram using Short-Time Fourier Transform

- Performs adaptive thresholding to detect syllables

- Validates detected segments against duration constraints

For SAP objects:

- Supports batch processing with parallel execution

- Processes specific days or indices

- Organizes results by source type

- Maintains metadata relationships

dB Scale Conversion:

- User input: 0 to 40 dB (intuitive positive scale)

- Internal conversion: Subtracts reference level (20 dB)

- Actual dBFS: -60 to -20 dB relative to full scale

Search Direction Guidelines:

- Use "up" when:

  - Recording has variable amplitude

  - Background noise is significant

  - Want to detect quieter syllables

- Use "down" when:

  - Recording is clean with good SNR

  - Want to avoid false positives

  - Syllables are consistently loud

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic segmentation of WAV file
syllables <- segment("song.wav")

# Custom parameters for clean recording
syllables <- segment("clean_song.wav",
                     search_direction = "down",
                     min_syllable_ms = 30,
                     max_syllable_ms = 150)

# Process specific days in SAP object
sap_obj <- segment(sap_object,
                   segment_type = "bouts",
                   day = c(30, 40),
                   cores = 4)

# Process with custom detection parameters
sap_obj <- segment(sap_object,
                   segment_type = "motifs",
                   min_level_db = 10,
                   max_level_db = 30,
                   save_plot = TRUE)
} # }
```
