# Changelog

## ASAP 0.3.4

### New Features

- `template_clips()`: Added feature for automated acoustic template
  matching.
- [`plot_traces()`](https://lxiao06.github.io/ASAP/reference/plot_traces.md),
  [`plot_clusters()`](https://lxiao06.github.io/ASAP/reference/plot_clusters.md):
  Expanded visualization capabilities with new plot types and arguments.
- Added comprehensive longitudinal analysis vignette tutorials for bout
  detection, motif detection, syllable segmentation, and syllable
  labeling.
- Built a new `pkgdown` website with improved navigation and flowchart
  logic.

### Bug Fixes & Improvements

- [`auto_label()`](https://lxiao06.github.io/ASAP/reference/auto_label.md),
  [`manual_label()`](https://lxiao06.github.io/ASAP/reference/manual_label.md):
  Refined workflows, fixed `dplyr` variable scoping errors, and improved
  argument naming (`label_type = "auto"` and `"manual"`).
- Cleaned up argument typos and reorganized default `NAMESPACE` imports.
