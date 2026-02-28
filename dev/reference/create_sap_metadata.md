# Create metadata for audio files recorded by SAP2011 (Sound Analysis Pro)

Creates a metadata data frame from WAV files in specified directories.
Parses filenames to extract information about bird ID, recording date,
and time.

## Usage

``` r
create_sap_metadata(
  base_path,
  subfolders_to_include = NULL,
  subfolders_to_exclude = c("templates", "temp_plots"),
  labels
)
```

## Arguments

- base_path:

  Character string specifying the base directory path

- subfolders_to_include:

  Character vector of subfolder names to include. If NULL, includes all
  subfolders except those in subfolders_to_exclude

- subfolders_to_exclude:

  Character vector of subfolder names to exclude. Default excludes
  'templates' and 'temp_plots'

- labels:

  Character vector of labels corresponding to each subfolder. Must match
  the length of subfolders

## Value

A data frame containing metadata with columns:

- filename:

  Name of the WAV file

- bird_id:

  Extracted bird identifier

- day_post_hatch:

  Day post hatch (from subfolder name)

- recording_date:

  Date of recording (MM-DD format)

- recording_time:

  Time of recording (HH:MM:SS format)

- label:

  Label assigned to the subfolder

## Examples

``` r
if (FALSE) { # \dontrun{
metadata <- create_sap_metadata(
  base_path = "path/to/recordings",
  subfolders_to_include = c("day1", "day2"),
  labels = c("pre", "post")
)
} # }
```
