# ASAP <a href="https://lxiao06.github.io/ASAP/"><img src="reference/figures/logo.png" align="right" height="90" alt="ASAP" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/LXiao06/ASAP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LXiao06/ASAP/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Overview

ASAP (Automated Sound Analysis Pipeline) is an R toolkit designed for the longitudinal analysis of birdsong development, specifically optimized for tracking and studying the long-term vocalization patterns of zebra finches.

The pipeline introduces a dedicated SAP object, which streamlines the processing of recordings from the SAP2011 system or other sound data organized with similar structures. Additionally, ASAP is built with flexibility in mind, allowing for the integration of other objects to ensure compatibility with various recording platforms.

---

## Installation

You can install the development version of ASAP from GitHub with:

``` r
require(remotes)
remotes::install_github("LXiao06/ASAP")
```

---

## Tutorials

Explore our comprehensive guides to get up and running with ASAP:

- [Overview: ASAP 101](articles/single_wav_analysis.html)
- [Motif Detection](articles/motif_detection.html)
- [Acoustic Feature Analysis](articles/acoustic_feature_analysis.html)
- [Constructing the SAP Object](articles/construct_sap_object.html)
- [Longitudinal Bout Detection](articles/longitudinal_bout_detection.html)
- [Longitudinal Motif Detection](articles/longitudinal_motif_detection.html)
- [Longitudinal Syllable Segmentation](articles/longitudinal_syllable_segmentation.html)
- [Syllable Labeling](articles/syllable_labeling.html)
- [Exporting Curated Song Clips](articles/exporting_song_clips.html)
