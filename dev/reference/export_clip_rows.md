# Core Clip Extraction and Writing Engine

Iterates over a prepared clip data frame, extracts each clip from its
source WAV file, applies optional amplitude normalisation or denoising,
and writes the result to WAV files or a single HDF5 archive. Returns a
metadata data frame describing every successfully written clip.

## Usage

``` r
export_clip_rows(
  clips,
  wav_dir,
  output_format,
  output_dir,
  hdf5_filename,
  metadata_filename,
  name_prefix,
  keep_source_file_name,
  amp_normalize,
  noise_reduction = FALSE,
  cores,
  overwrite,
  write_metadata,
  verbose
)
```

## Arguments

- clips:

  Prepared clip data frame (output of `prepare_clip_export_rows`).

- wav_dir:

  Root directory containing source WAV files.

- output_format:

  One of `"wav"` or `"hdf5"`.

- output_dir:

  Directory where clips are written.

- hdf5_filename:

  File name for the HDF5 archive.

- metadata_filename:

  File name for the metadata CSV.

- name_prefix:

  Prefix for generated clip file names.

- keep_source_file_name:

  Logical; use WAV stem + index as clip name.

- amp_normalize:

  Amplitude normalisation method (`"none"`, `"peak"`, or `"rms"`).

- noise_reduction:

  Logical; denoise source files before extraction.

- cores:

  Number of parallel cores.

- overwrite:

  Logical; overwrite existing files.

- write_metadata:

  Logical; write a metadata CSV.

- verbose:

  Logical; print per-day progress messages.

## Value

Data frame of export metadata (one row per successfully written clip).
