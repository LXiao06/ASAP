# Longitudinal Motif Detection

## Introduction

This vignette demonstrates how to detect and analyze song motifs across
**longitudinal recordings** using SAP objects and optimized template
parameters.

**Prerequisites**: Before reading this vignette, we recommend
completing:

- [Overview: Basic Audio
  Analysis](https://lxiao06.github.io/ASAP/articles/single_wav_analysis.md) -
  Core ASAP functions
- [Motif
  Detection](https://lxiao06.github.io/ASAP/articles/motif_detection.md) -
  Template optimization workflow
- [Constructing SAP
  Object](https://lxiao06.github.io/ASAP/articles/construct_sap_object.md) -
  SAP object creation

## Overview

The longitudinal motif detection workflow applies the template
parameters you optimized on a single recording (see [Motif
Detection](https://lxiao06.github.io/ASAP/articles/motif_detection.md))
across all recordings in a SAP object.

![Longitudinal Motif Detection
Pipeline](figures/longitudinal_pipeline.png)

## Complete Pipeline

``` r
library(ASAP)

# Load or create SAP object
sap <- create_sap_object(
  base_path = "/path/to/recordings",
  subfolders_to_include = c("190", "201", "203"),
  labels = c("BL", "Post", "Rec")
)

# Run the complete motif detection pipeline
sap <- sap |>
  create_audio_clip(indices = 1, 
                    start_time = 1, 
                    end_time = 2.5, 
                    clip_names = "motif_ref") |>
  create_template(template_name = "syllable_d", 
                  clip_name = "motif_ref",
                  start_time = 0.72, 
                  end_time = 0.84,
                  freq_min = 1, 
                  freq_max = 10, 
                  threshold = 0.5, 
                  write_template = TRUE) |>
  detect_template(template_name = "syllable_d",
                  threshold = 0.5,
                  proximity_window = 1) |>
  find_motif(template_name = "syllable_d", 
             pre_time = 0.7, 
             lag_time = 0.5)
```

### Understanding SAP Object Behavior

**1. Metadata-Based Lazy Loading**

After creating the SAP object, you have a metadata index (stored in
`sap$metadata`) that references all audio files by their paths and
timestamps. The actual WAV file data is **not loaded into memory**—files
are read on-demand during analysis. This means:

- The SAP object remains lightweight even with thousands of recordings
- **All subsequent pipeline steps require the original WAV files to be
  accessible** at their original paths
- If you move or delete the audio files, the pipeline will fail

**2. Centralized Result Storage**

The SAP object stores all analysis results in a structured format,
making it easy to access detection outcomes:

``` r
# Access template detection results
sap$templates$template_matches[["syllable_d"]]

# Access extracted motif boundaries (onset/offset timestamps)
sap$motifs  # Data frame with start_time, end_time for each detected motif

# Access spectral features (if analyze_spectral() was run)
sap$features$motif$spectral_feature
```

**3. Optional Template Detection Visualizations**

If you set `save_plot = TRUE` in
[`detect_template()`](https://lxiao06.github.io/ASAP/reference/detect_template.md),
spectrogram images of template detections are saved to a local directory
(typically `templates/` within your base path). These images are **not
stored in the SAP object** but can be useful for quality control and
manual inspection of detection accuracy.

**Example output:**

    === Starting Template Detection ===

    Processing 312 files for day 190 using 7 cores.
    Processed files in day 190. Total detections: 1247

    Processing 285 files for day 201 using 7 cores.
    Processed files in day 201. Total detections: 1089

    Processing 250 files for day 203 using 7 cores.
    Processed files in day 203. Total detections: 956

    Total detections across all days: 3292
    Access detection results via: sap$templates$template_matches[["syllable_d"]]

## Visualizing Results

``` r
# Visualize sample motifs from each time point
visualize_segments(sap, 
                   segment_type = "motifs", 
                   n_samples = 3)
```

![Sample motif spectrograms across developmental time
points.](figures/longitudinal_segments.png)

Sample motif spectrograms across developmental time points.

### Amplitude Envelope Heatmap

Amplitude envelope heatmaps visualize the temporal structure of detected
motifs. **The `balanced` argument** ensures equal representation across
time points (e.g., same number of motifs from each developmental stage),
which improves statistical comparisons.

However, **motif durations can vary** across renditions—some motifs may
be shorter than others due to natural variability or developmental
changes. This duration variability can make heatmaps appear misaligned
or noisy, with vertical bands appearing blurred or offset. To address
this issue, we can extract acoustic features from individual motifs,
cluster them based on these features within each developmental stage,
and order them by their latent feature structure (see “Ordered Heatmap”
section below).

``` r
# Create amplitude envelope heatmap with balanced sampling
plot_heatmap(sap, balanced = TRUE)
```

![Amplitude envelope heatmap showing temporal structure across time
points.](figures/longitudinal_heatmap.png)

Amplitude envelope heatmap showing temporal structure across time
points.

## Feature Extraction and Analysis

Extracting spectral features allows you to quantify motif acoustic
properties and identify clusters of similar motifs. This is particularly
useful for:

- **Grouping motifs by acoustic similarity** rather than just temporal
  alignment
- **Identifying developmental changes** in song structure
- **Reducing noise** in downstream visualizations

``` r
# Extract spectral features (frequency, entropy, duration, etc.)
sap <- sap |>
  analyze_spectral(balanced = TRUE) |>
  find_clusters() |>              # Group acoustically similar motifs
  run_umap()                       # Dimensionality reduction for visualization

# Visualize UMAP by time point
plot_umap(sap, split.by = "label")
```

**Results interpretation:**

- [`analyze_spectral()`](https://lxiao06.github.io/ASAP/reference/analyze_spectral.md)
  extracts ~30 acoustic features per motif (mean frequency, entropy,
  duration, spectral slope, etc.) and stores them in
  `sap$features$motif$spectral_feature`
- [`find_clusters()`](https://lxiao06.github.io/ASAP/reference/find_clusters.md)
  uses hierarchical clustering to group motifs by acoustic similarity,
  adding a `cluster` column to the feature data
- [`run_umap()`](https://lxiao06.github.io/ASAP/reference/run_umap.md)
  reduces high-dimensional feature space to 2D for visualization,
  storing coordinates in `sap$umap$motif`

![UMAP visualization of motif features colored by developmental time
point.](figures/longitudinal_umap.png)

UMAP visualization of motif features colored by developmental time
point.

### Ordered Heatmap by Acoustic Similarity

After extracting features and clustering, motifs are ordered by their
cluster membership and latent feature structure rather than
chronological order. Typically, ordering motifs by acoustic similarity
creates a **cleaner, more organized heatmap** that addresses the visual
noise from duration variability.

``` r
# Create heatmap ordered by cluster membership
plot_heatmap(sap, balanced = TRUE, ordered = TRUE)
```

**Why this helps:**

- Motifs with similar acoustic properties are grouped together
- Vertical bands align better because acoustically similar motifs tend
  to have similar durations
- Developmental patterns become more apparent when organized by acoustic
  structure
- Reduces visual noise from duration variability

![Ordered amplitude envelope heatmap grouped by acoustic
similarity.](figures/longitudinal_heatmap_ordered.png)

Ordered amplitude envelope heatmap grouped by acoustic similarity.

## Key Parameters for Longitudinal Analysis

| Parameter          | Location                                                                                                                                                                | Description                                                                                                                                                    |
|--------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `threshold`        | [`create_template()`](https://lxiao06.github.io/ASAP/reference/create_template.md) / [`detect_template()`](https://lxiao06.github.io/ASAP/reference/detect_template.md) | Minimum correlation score. Adjust in [`detect_template()`](https://lxiao06.github.io/ASAP/reference/detect_template.md) to refine without recreating template. |
| `proximity_window` | [`detect_template()`](https://lxiao06.github.io/ASAP/reference/detect_template.md)                                                                                      | Filter duplicate detections within this time window (seconds).                                                                                                 |
| `balanced`         | [`analyze_spectral()`](https://lxiao06.github.io/ASAP/reference/analyze_spectral.md) / [`plot_heatmap()`](https://lxiao06.github.io/ASAP/reference/plot_heatmap.md)     | Balance samples across time points.                                                                                                                            |

## Tips for Longitudinal Analysis

1.  **Optimize parameters first**: Use the single-file workflow in
    [Motif
    Detection](https://lxiao06.github.io/ASAP/articles/motif_detection.md)
    before bulk processing.

2.  **Check detection quality**: Visualize sample detections from each
    time point to verify template works across developmental stages.

3.  **Adjust threshold if needed**: If detection rates vary
    significantly across time points, consider adjusting the threshold
    in
    [`detect_template()`](https://lxiao06.github.io/ASAP/reference/detect_template.md).

4.  **Use balanced sampling**: When comparing across time points, use
    `balanced = TRUE` to ensure equal representation.

## Session Info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.56         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.30    lifecycle_1.0.5   cli_3.6.5         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1
#> [17] compiler_4.5.2    tools_4.5.2       ragg_1.5.0        evaluate_1.0.5   
#> [21] bslib_0.10.0      yaml_2.3.12       jsonlite_2.0.0    rlang_1.1.7      
#> [25] fs_1.6.6
```
