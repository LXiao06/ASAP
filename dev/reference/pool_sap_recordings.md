# Pool Selected Recording Folders into a Standard SAP Object

Links or copies selected recording folders into one conventional SAP
folder tree, then creates a regular SAP object. This lets the existing
pipeline run without multi-animal path handling.

## Usage

``` r
pool_sap_recordings(
  selections,
  output_dir,
  method = c("link", "copy"),
  overwrite = FALSE
)
```

## Arguments

- selections:

  Data frame with character columns `animal_id`, `base_path`,
  `subfolder`, and `label`; one row per selected folder.

- output_dir:

  Directory for the pooled recording tree.

- method:

  Either `"link"` (default) to create symbolic links or `"copy"` to copy
  WAV files.

- overwrite:

  Logical; replace existing files in the pooled tree.

## Value

A regular SAP object. Its day folders are named
`<animal_id>__<subfolder>` to keep source files distinct.

## Examples

``` r
if (FALSE) { # \dontrun{
sap <- pool_sap_recordings(data.frame(
  animal_id = c("bird_1", "bird_2"),
  base_path = c("/data/bird_1", "/data/bird_2"),
  subfolder = c("60", "65"),
  label = c("60 dph", "65 dph")
), output_dir = "/data/pooled_recordings")
} # }
```
