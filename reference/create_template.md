# Create Correlation Templates for Song Analysis

Creates correlation templates from WAV files or SAP objects for song
detection and analysis.

## Usage

``` r
create_template(x, ...)

# Default S3 method
create_template(
  x,
  template_name,
  start_time = NULL,
  end_time = NULL,
  freq_min = 0,
  freq_max = 15,
  threshold = 0.6,
  write_template = FALSE,
  ...
)

# S3 method for class 'Sap'
create_template(
  x,
  template_name,
  clip_name,
  start_time = NULL,
  end_time = NULL,
  freq_min = 0,
  freq_max = 15,
  threshold = 0.6,
  write_template = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to process, either a file path or SAP object

- ...:

  Additional arguments passed to specific methods

- template_name:

  Character name for the template

- start_time:

  Numeric start time of template segment

- end_time:

  Numeric end time of template segment

- freq_min:

  Numeric minimum frequency in kHz (default: 0)

- freq_max:

  Numeric maximum frequency in kHz (default: 15)

- threshold:

  Numeric correlation threshold (default: 0.6)

- write_template:

  Logical whether to write template to disk

- clip_name:

  For SAP objects: Character name of the clip to use

- verbose:

  For SAP objects: Whether to print progress messages

## Value

For default method: A correlation template object from monitoR package
For SAP objects: Updated SAP object with new template information

## Details

For WAV files:

- Validates input parameters

- Creates template using monitoR package

- Optionally writes template to disk

For SAP objects:

- Validates clip existence

- Creates template using specified parameters

- Updates SAP object with template information

## See also

[`create_audio_clip`](https://lxiao06.github.io/ASAP/reference/create_audio_clip.md)
for creating audio clips

## Examples

``` r
if (FALSE) { # \dontrun{
# Create template from WAV file
template <- create_template("path/to/song.wav",
                           template_name = "template1",
                           start_time = 1.0,
                           end_time = 2.0,
                           freq_min = 2,
                           freq_max = 8)

# Create and save template
template <- create_template("song.wav",
                           template_name = "template2",
                           start_time = 1.0,
                           end_time = 2.0,
                           write_template = TRUE)

# Create template from SAP object
sap_obj <- create_template(sap_object,
                          template_name = "template1",
                          clip_name = "clip1",
                          freq_min = 2,
                          freq_max = 8)
} # }
```
