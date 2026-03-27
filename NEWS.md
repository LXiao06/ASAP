# ASAP 0.3.5 (2026-03-27)

## New Features
- `create_bout_clips()`, `create_motif_clips()`: Added new workflows for exporting curated bout- and motif-level audio clips as WAV or HDF5 files, with support for metadata export, amplitude normalization, sampling, source-aware clip naming, and a new `margin` parameter for expanding bout boundaries before export.
- `denoise()`: Added audio denoising for single WAV files and SAP objects using spectral median subtraction or spectral gating.
- Added a comprehensive song clip exporting vignette with end-to-end examples for compression, filtered-bout export, and motif-level clip extraction.

## Bug Fixes & Improvements
- `parallel_apply()`: Improved parallel processing efficiency on Linux and Windows by reusing PSOCK workers and explicitly loading the `ASAP` namespace on workers to avoid repeated startup overhead.
- `detect_template()`, `denoise()`, `find_bout()`, `segment()`: Fixed multi-core processing issues by adopting the updated PSOCK workflow, and suppressed C-level `fontconfig` warnings during template matching on Linux.
- Refined clip export validation and internal documentation, including motif filtering, time-boundary checks, and CRAN portability cleanups.

# ASAP 0.3.4 (2026-02-27)

## New Features
- `template_clips()`: Added feature for automated acoustic template matching.
- `plot_traces()`, `plot_clusters()`: Expanded visualization capabilities with new plot types and arguments.
- Added comprehensive longitudinal analysis vignette tutorials for bout detection, motif detection, syllable segmentation, and syllable labeling.
- Built a new `pkgdown` website with improved navigation and flowchart logic.

## Bug Fixes & Improvements
- `auto_label()`, `manual_label()`: Refined workflows, fixed `dplyr` variable scoping errors, and improved argument naming (`label_type = "auto"` and `"manual"`).
- Cleaned up argument typos and reorganized default `NAMESPACE` imports.
