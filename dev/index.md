# ASAP

## Overview

ASAP (Automated Sound Analysis Pipeline) is an R toolkit designed for
the longitudinal analysis of birdsong development, specifically
optimized for tracking and studying the long-term vocalization patterns
of zebra finches.

The pipeline introduces a dedicated SAP object, which streamlines the
processing of recordings from the SAP2011 system or other sound data
organized with similar structures. Additionally, ASAP is built with
flexibility in mind, allowing for the integration of other objects to
ensure compatibility with various recording platforms.

### Key Features

- **Bout Detection**: Automatically identifies periods of singing within
  recordings.

- **Motif Extraction**: Extracts recurring song motifs for in-depth
  analysis.

- **Syllable Segmentation**: Breaks down songs into individual syllables
  for detailed study.

- **Template Matching**: Compares and matches song patterns against
  predefined templates.

- **Feature Extraction**: Extracts detailed temporal and spectral
  features from vocalizations.

- **Standardized Analytical Workflows**: Utilizes predefined workflows
  for consistent and efficient analysis.

## Installation

To install ASAP, I recommend using remotes:

``` r
require(remotes)
remotes::install_github("LXiao06/ASAP")
```

## Tutorials

Explore our comprehensive guides to get up and running with ASAP. You
can find all documentation on the [ASAP
website](https://lxiao06.github.io/ASAP/).

### Quick Start

- [Single WAV File
  Analysis](https://lxiao06.github.io/ASAP/articles/single_wav_analysis.html):
  Learn ASAP’s core functions using individual audio recordings.
- [Motif
  Detection](https://lxiao06.github.io/ASAP/articles/motif_detection.html):
  Introduction to automated motif extraction using template matching.

### Longitudinal Analysis Pipeline

Step-by-step guides for processing large longitudinal datasets spanning
multiple developmental time points:

1.  [Constructing the SAP
    Object](https://lxiao06.github.io/ASAP/articles/construct_sap_object.html)
2.  [Longitudinal Bout
    Detection](https://lxiao06.github.io/ASAP/articles/longitudinal_bout_detection.html)
3.  [Longitudinal Motif
    Detection](https://lxiao06.github.io/ASAP/articles/longitudinal_motif_detection.html)
4.  [Longitudinal Syllable
    Segmentation](https://lxiao06.github.io/ASAP/articles/longitudinal_syllable_segmentation.html)
5.  [Syllable
    Labeling](https://lxiao06.github.io/ASAP/articles/syllable_labeling.html)
