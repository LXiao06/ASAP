# Parse filename for SAP metadata

Internal function to parse WAV filenames and extract metadata
components.

## Usage

``` r
parse_filename(filename)
```

## Arguments

- filename:

  Character string of the filename to parse

## Value

A list containing parsed components:

- bird_id:

  Extracted bird identifier

- recording_date:

  Parsed recording date

- recording_time:

  Parsed recording time
