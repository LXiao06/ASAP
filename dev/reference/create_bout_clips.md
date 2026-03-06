# Create Bout Audio Clips

Exports bout-level clips detected by
[`find_bout()`](https://lxiao06.github.io/ASAP/dev/reference/find_bout.md)
either as WAV files in a directory tree or as datasets in a single HDF5
file. Bouts typically span multiple motifs and are longer than
individual motif clips.

## Usage

``` r
create_bout_clips(x, ...)

# Default S3 method
create_bout_clips(
  x,
  wav_dir,
  indices = NULL,
  n_bouts = NULL,
  seed = 222,
  output_format = c("wav", "hdf5"),
  output_dir = NULL,
  hdf5_filename = "bouts.h5",
  metadata_filename = "metadata.csv",
  name_prefix = "bout",
  amp_normalize = c("none", "peak", "rms"),
  cores = NULL,
  overwrite = TRUE,
  write_metadata = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
create_bout_clips(
  x,
  indices = NULL,
  n_bouts = NULL,
  seed = 222,
  output_format = c("wav", "hdf5"),
  output_dir = NULL,
  hdf5_filename = "bouts.h5",
  metadata_filename = "metadata.csv",
  name_prefix = "bout",
  amp_normalize = c("none", "peak", "rms"),
  cores = NULL,
  overwrite = TRUE,
  write_metadata = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A SAP object or a bouts data frame.

- ...:

  Additional arguments passed to methods.

- wav_dir:

  For data-frame method: base directory containing source WAV files.

- indices:

  Optional numeric vector of bout row indices to export. If NULL, all
  bouts are exported.

- n_bouts:

  Optional integer. Number of bouts to randomly sample when `indices` is
  NULL. Sampling is applied per day/subdirectory (`day_post_hatch`) when
  available. If NULL, all bouts are processed.

- seed:

  Integer random seed used when `n_bouts` sampling is applied. Default
  is `222`.

- output_format:

  Output type: `"wav"` or `"hdf5"`.

- output_dir:

  Directory where output files are written.

- hdf5_filename:

  File name used when `output_format = "hdf5"`.

- metadata_filename:

  Name of metadata CSV written to `output_dir`.

- name_prefix:

  Prefix used for generated bout clip names.

- amp_normalize:

  Waveform amplitude normalization applied to exported clips: one of
  "none", "peak", or "rms" (default: "none").

- cores:

  Number of CPU cores used for clip processing.

- overwrite:

  Logical. Overwrite existing output file(s) when TRUE. Default is TRUE.

- write_metadata:

  Logical. Write metadata CSV when TRUE.

- verbose:

  Logical. Print one export summary per day and overall totals.

## Value

For SAP input: updated SAP object with export summary in
`x$misc$bout_clip_exports`. For data-frame input: metadata data frame.

## Details

Output layout for `output_format = "wav"`:

    output_dir/bouts/{bird_id}/{day_post_hatch}/{name_prefix}_001.wav

Output layout for `output_format = "hdf5"`:

    output_dir/{hdf5_filename}
      /{bird_id}/{day_post_hatch}/{name_prefix}_001

If `bird_id` is not available in bout rows, the function attempts to
infer it from metadata. Missing values are stored as `"unknown_bird"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Export up to 50 bouts per day as WAV files
sap <- sap |>
  find_bout(min_duration = 0.4, summary = TRUE) |>
  create_bout_clips(
    n_bouts = 50,
    output_format = "wav",
    output_dir = "exported_bouts"
  )

# Export all bouts to a single HDF5 file
sap <- create_bout_clips(
  sap,
  output_format = "hdf5",
  output_dir    = "exported_bouts",
  hdf5_filename = "bouts.h5"
)
} # }
```
