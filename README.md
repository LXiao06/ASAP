---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# ASAP <a href="https://lxiao06.github.io/ASAP/"><img src="man/figures/logo.png" align="right" height="90" alt="ASAP" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/LXiao06/ASAP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LXiao06/ASAP/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Overview

ASAP (Automated Sound Analysis Pipeline) is an R toolkit designed for the longitudinal analysis of birdsong development, specifically optimized for tracking and studying the long-term vocalization patterns of zebra finches.

The pipeline introduces a dedicated SAP object, which streamlines the processing of recordings from the SAP2011 system or other sound data organized with similar structures. Additionally, ASAP is built with flexibility in mind, allowing for the integration of other objects to ensure compatibility with various recording platforms.

### Key Features

-   **Bout Detection**: Automatically identifies periods of singing within recordings.

-   **Motif Extraction**: Extracts recurring song motifs for in-depth analysis.

-   **Syllable Segmentation**: Breaks down songs into individual syllables for detailed study.

-   **Template Matching**: Compares and matches song patterns against predefined templates.

-   **Feature Extraction**: Extracts detailed temporal and spectral features from vocalizations.

-   **Standardized Analytical Workflows**: Utilizes predefined workflows for consistent and efficient analysis.

## Installation

You can install the latest stable release of ASAP from GitHub:


``` r
# Install the latest stable release (recommended)
remotes::install_github("LXiao06/ASAP@*release")
```

Or install the development version with the latest features:


``` r
# Install the development version
remotes::install_github("LXiao06/ASAP")
```

## Tutorials

Explore our guides to get up and running with ASAP. You can find all documentation on the [ASAP website](https://lxiao06.github.io/ASAP/).

### Quick Start
- [ASAP 101](https://lxiao06.github.io/ASAP/articles/single_wav_analysis.html): Learn the single-file workflow with spectrograms, bout detection and export, and syllable segmentation.
- [Motif Detection](https://lxiao06.github.io/ASAP/articles/motif_detection.html): Introduction to automated motif extraction using template matching.
- [Acoustic Feature Analysis](https://lxiao06.github.io/ASAP/articles/acoustic_feature_analysis.html): Measure spectral entropy, pitch contours, and amplitude envelopes.

### Longitudinal Analysis Pipeline
Step-by-step guides for processing large longitudinal datasets spanning multiple developmental time points:

1. [Constructing the SAP Object](https://lxiao06.github.io/ASAP/articles/construct_sap_object.html)
2. [Longitudinal Bout Detection](https://lxiao06.github.io/ASAP/articles/longitudinal_bout_detection.html)
3. [Longitudinal Motif Detection](https://lxiao06.github.io/ASAP/articles/longitudinal_motif_detection.html)
4. [Longitudinal Syllable Segmentation](https://lxiao06.github.io/ASAP/articles/longitudinal_syllable_segmentation.html)
5. [Syllable Labeling](https://lxiao06.github.io/ASAP/articles/syllable_labeling.html)
6. [Exporting Curated Song Clips](https://lxiao06.github.io/ASAP/articles/exporting_song_clips.html)

## Acknowledgements

ASAP relies on the following R packages:
[seewave](https://cran.r-project.org/web/packages/seewave/index.html),
[tuneR](https://cran.r-project.org/web/packages/tuneR/index.html),
[warbleR](https://cran.r-project.org/web/packages/warbleR/index.html),
[monitoR](https://cran.r-project.org/web/packages/monitoR/index.html), and
[av](https://cran.r-project.org/web/packages/av/index.html).

This work was inspired by the pioneering contributions of
**Ofer Tchernichovski et al.** ([Sound Analysis Pro, SAP 2011](http://soundanalysispro.com/)),
**Tim Sainburg** ([timsainb](https://github.com/timsainb)),
**Therese Koch** ([AVN](https://avn.readthedocs.io/en/latest/index.html)), and
**Seurat** ([Butler et al.](https://cran.r-project.org/web/packages/Seurat/index.html)).
