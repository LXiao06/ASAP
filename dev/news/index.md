# Changelog

## ASAP 0.3.5 (2026-03-27)

### New Features

- [`create_bout_clips()`](https://lxiao06.github.io/ASAP/dev/reference/create_bout_clips.md),
  [`create_motif_clips()`](https://lxiao06.github.io/ASAP/dev/reference/create_motif_clips.md):
  Added new workflows for exporting curated bout- and motif-level audio
  clips as WAV or HDF5 files, with support for metadata export,
  amplitude normalization, sampling, source-aware clip naming, and a new
  `margin` parameter for expanding bout boundaries before export.
- [`denoise()`](https://lxiao06.github.io/ASAP/dev/reference/denoise.md):
  Added audio denoising for single WAV files and SAP objects using
  spectral median subtraction or spectral gating.
- Added a comprehensive song clip exporting vignette with end-to-end
  examples for compression, filtered-bout export, and motif-level clip
  extraction.

### Bug Fixes & Improvements

- [`parallel_apply()`](https://lxiao06.github.io/ASAP/dev/reference/parallel_apply.md):
  Improved parallel processing efficiency on Linux and Windows by
  reusing PSOCK workers and explicitly loading the `ASAP` namespace on
  workers to avoid repeated startup overhead.
- [`detect_template()`](https://lxiao06.github.io/ASAP/dev/reference/detect_template.md),
  [`denoise()`](https://lxiao06.github.io/ASAP/dev/reference/denoise.md),
  [`find_bout()`](https://lxiao06.github.io/ASAP/dev/reference/find_bout.md),
  [`segment()`](https://lxiao06.github.io/ASAP/dev/reference/segment.md):
  Fixed multi-core processing issues by adopting the updated PSOCK
  workflow, and suppressed C-level `fontconfig` warnings during template
  matching on Linux.
- Refined clip export validation and internal documentation, including
  motif filtering, time-boundary checks, and CRAN portability cleanups.

## ASAP 0.3.4 (2026-02-27)

### New Features

- `template_clips()`: Added feature for automated acoustic template
  matching.
- [`plot_traces()`](https://lxiao06.github.io/ASAP/dev/reference/plot_traces.md),
  [`plot_clusters()`](https://lxiao06.github.io/ASAP/dev/reference/plot_clusters.md):
  Expanded visualization capabilities with new plot types and arguments.
- Added comprehensive longitudinal analysis vignette tutorials for bout
  detection, motif detection, syllable segmentation, and syllable
  labeling.
- Built a new `pkgdown` website with improved navigation and flowchart
  logic.

### Bug Fixes & Improvements

- [`auto_label()`](https://lxiao06.github.io/ASAP/dev/reference/auto_label.md),
  [`manual_label()`](https://lxiao06.github.io/ASAP/dev/reference/manual_label.md):
  Refined workflows, fixed `dplyr` variable scoping errors, and improved
  argument naming (`label_type = "auto"` and `"manual"`).
- Cleaned up argument typos and reorganized default `NAMESPACE` imports.
