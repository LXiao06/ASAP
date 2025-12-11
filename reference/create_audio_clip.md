# Create Audio Clips from Sound Files

Creates audio clips from WAV files or SAP objects by extracting
specified time segments.

## Usage

``` r
create_audio_clip(x, ...)

# Default S3 method
create_audio_clip(
  x,
  ...,
  start_time,
  end_time,
  clip_name = NULL,
  unit = "second"
)

# S3 method for class 'Sap'
create_audio_clip(
  x,
  indices,
  start_time,
  end_time,
  clip_names,
  unit = "second",
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to process, either a file path or SAP object

- ...:

  Additional arguments passed to specific methods

- start_time:

  Numeric start time(s) of the clip(s)

- end_time:

  Numeric end time(s) of the clip(s)

- clip_name, clip_names:

  Name(s) for the output clip(s)

- unit:

  Time unit ("second" or "millisecond")

- indices:

  For SAP objects: Numeric vector of indices to process

- verbose:

  For SAP objects: Whether to print progress messages

## Value

For default method: Character string containing path to created audio
clip For SAP objects: Updated SAP object with new template information

Updated SAP object with new template information

## Details

For single WAV files:

- Validates input file and time parameters

- Creates templates directory if needed

- Extracts specified segment from audio file

For SAP objects:

- Creates clips for specified indices

- Updates template information in SAP object

- Maintains metadata about created clips

## Examples

``` r
if (FALSE) { # \dontrun{
# Create clip from single WAV file
create_audio_clip("path/to/song.wav",
                  start_time = 10,
                  end_time = 20,
                  clip_name = "song_clip")

# Create multiple clips from SAP object
create_audio_clip(sap_object,
                  indices = c(1, 2),
                  start_time = c(10, 20),
                  end_time = c(20, 30),
                  clip_names = c("clip1", "clip2"))

# Create clip with millisecond units
create_audio_clip("song.wav",
                  start_time = 10000,
                  end_time = 20000,
                  unit = "millisecond")
} # }
```
