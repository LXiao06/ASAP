# Parse filename for SAP metadata

Internal function to parse WAV filenames and extract metadata
components. Supports the standard SAP2011 filename format and falls back
to a generic parser for non-standard filenames, using file modification
time for date/time.

## Usage

``` r
parse_filename(filename)
```

## Arguments

- filename:

  Character string of the full file path to parse

## Value

A list containing parsed components:

- bird_id:

  Extracted bird identifier

- recording_date:

  Parsed recording date

- recording_time:

  Parsed recording time

- parse_method:

  "standard" or "fallback"
