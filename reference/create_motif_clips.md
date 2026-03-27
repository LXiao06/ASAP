# Create Motif Audio Clips

Exports motif-level clips detected by
[`find_motif()`](https://lxiao06.github.io/ASAP/reference/find_motif.md)
either as WAV files in a directory tree or as datasets in a single HDF5
file.

## Usage

``` r
create_motif_clips(x, ...)

# Default S3 method
create_motif_clips(
  x,
  wav_dir,
  indices = NULL,
  n_motifs = NULL,
  seed = 222,
  output_format = c("wav", "hdf5"),
  output_dir = NULL,
  hdf5_filename = "motifs.h5",
  metadata_filename = "metadata.csv",
  name_prefix = NULL,
  keep_source_file_name = FALSE,
  amp_normalize = c("none", "peak", "rms"),
  noise_reduction = FALSE,
  cores = NULL,
  overwrite = TRUE,
  write_metadata = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
create_motif_clips(
  x,
  indices = NULL,
  n_motifs = NULL,
  seed = 222,
  output_format = c("wav", "hdf5"),
  output_dir = NULL,
  hdf5_filename = "motifs.h5",
  metadata_filename = "metadata.csv",
  name_prefix = NULL,
  keep_source_file_name = FALSE,
  amp_normalize = c("none", "peak", "rms"),
  noise_reduction = FALSE,
  cores = NULL,
  overwrite = TRUE,
  write_metadata = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A SAP object or a motif data frame.

- ...:

  Additional arguments passed to methods.

- wav_dir:

  For data-frame method: base directory containing source WAV files.

- indices:

  Optional numeric vector of motif row indices to export. If NULL, all
  motifs are exported.

- n_motifs:

  Optional integer. Number of motifs to randomly sample when `indices`
  is NULL. Sampling is applied per day/subdirectory (`day_post_hatch`)
  when available. If NULL, all motifs are processed.

- seed:

  Integer random seed used when `n_motifs` sampling is applied. In this
  repository, sampling is fixed to `222`.

- output_format:

  Output type: `"wav"` or `"hdf5"`.

- output_dir:

  Directory where output files are written.

- hdf5_filename:

  File name used when `output_format = "hdf5"`.

- metadata_filename:

  Name of metadata CSV written to `output_dir`.

- name_prefix:

  Prefix used for generated clip file names. Can be any string (e.g.
  `"motif"` or `"bout"`). If `NULL`, defaults to `"motif"`.

- keep_source_file_name:

  Logical. If `TRUE`, uses the stem of the originating WAV file combined
  with the `selec` column (or sequential index) as the clip identifier
  (e.g. `S237_42685_001.wav`). This option is especially useful for
  Scenario A exports where tracing a clip back to its original recording
  is important. Overrides `name_prefix`.

- amp_normalize:

  Waveform amplitude normalization applied to exported clips: one of
  "none", "peak", or "rms" (default: "none")

- noise_reduction:

  Logical. If TRUE, denoise source WAV files before extracting motif
  clips. Denoised sources are written under
  `output_dir/motifs/denoised_sources` and used for clip extraction.

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
`x$misc$motif_clip_exports`. For data-frame input: metadata data frame.

## Details

Output layout for `output_format = "wav"`:

    output_dir/motifs/{bird_id}/{day_post_hatch}/{name_prefix}_001.wav

Output layout for `output_format = "hdf5"`:

    output_dir/{hdf5_filename}
      /{bird_id}/{day_post_hatch}/{name_prefix}_001

If `bird_id` is not available in motif rows, the function attempts to
infer it from metadata. Missing values are stored as `"unknown_bird"`.

## Examples

``` r
if (FALSE) { # \dontrun{
sap <- sap |>
  find_motif(template_name = "syllable_d", pre_time = 0.7, lag_time = 0.5) |>
  create_motif_clips(indices = 1:50, output_format = "wav")

sap <- create_motif_clips(
  sap,
  output_format = "hdf5",
  hdf5_filename = "motifs.h5"
)

# Export RMS-normalized motif clips
sap <- create_motif_clips(
  sap,
  output_format = "wav",
  amp_normalize = "rms"
)
} # }
```
