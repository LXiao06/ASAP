# Compute WAV File Durations

Calculates durations for all WAV files referenced in a SAP object's
metadata using parallel processing. Handles missing files gracefully by
returning NA durations and providing warnings.

## Usage

``` r
compute_wav_durations(x, cores = NULL, verbose = TRUE)
```

## Arguments

- x:

  A SAP object containing metadata and base path to audio files

- cores:

  Number of cores to use for parallel processing (NULL for
  auto-detection: total cores - 1)

- verbose:

  Logical flag to control progress messages and warnings (default: TRUE)

## Value

Returns the modified SAP object with added duration column in metadata
containing wave file durations in seconds. Missing files will have NA
durations.

## Details

Key features:

- Parallel processing implementation using `parallel_apply`

- Automatic core detection with fallback to single-core processing

- Progress tracking and missing file warnings

- Preserves original object structure while adding duration information
