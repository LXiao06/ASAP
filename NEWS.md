# ASAP 0.3.4 (2026-02-27)

## New Features
- `template_clips()`: Added feature for automated acoustic template matching.
- `plot_traces()`, `plot_clusters()`: Expanded visualization capabilities with new plot types and arguments.
- Added comprehensive longitudinal analysis vignette tutorials for bout detection, motif detection, syllable segmentation, and syllable labeling.
- Built a new `pkgdown` website with improved navigation and flowchart logic.

## Bug Fixes & Improvements
- `auto_label()`, `manual_label()`: Refined workflows, fixed `dplyr` variable scoping errors, and improved argument naming (`label_type = "auto"` and `"manual"`).
- Cleaned up argument typos and reorganized default `NAMESPACE` imports.
