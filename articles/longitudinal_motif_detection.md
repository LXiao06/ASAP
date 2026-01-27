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
  Objects](https://lxiao06.github.io/ASAP/articles/construct_sap_object.md) -
  SAP object creation

## Overview

The longitudinal motif detection workflow applies the template
parameters you optimized on a single recording (see [Motif
Detection](https://lxiao06.github.io/ASAP/articles/motif_detection.md))
across all recordings in a SAP object.

## Complete Pipeline

``` r
library(ASAP)

# Load or create SAP object
sap <- create_sap_object(
  base_path = "/path/to/recordings",
  subfolders_to_include = c("190", "201", "203"),
  labels = c("Baseline", "Post", "Recovery")
)

# Run the complete motif detection pipeline
sap <- sap |>
  # Step 1: Create audio clip from reference recording
  create_audio_clip(indices = 1, 
                    start_time = 1, 
                    end_time = 2.5, 
                    clip_names = "motif_ref") |>
  
  # Step 2: Create template with optimized parameters
  create_template(template_name = "syllable_d", 
                  clip_name = "motif_ref",
                  start_time = 0.72, 
                  end_time = 0.84,
                  freq_min = 1, 
                  freq_max = 10, 
                  threshold = 0.5, 
                  write_template = TRUE) |>
  
  # Step 3: Detect template across all recordings
  detect_template(template_name = "syllable_d",
                  threshold = 0.5,
                  proximity_window = 1) |>
  
  # Step 4: Extract motif boundaries
  find_motif(template_name = "syllable_d", 
             pre_time = 0.7, 
             lag_time = 0.5)
```

## Visualizing Results

``` r
# View detection summary
summary(sap$motifs)

# Visualize sample motifs from each time point
visualize_segments(sap, 
                   segment_type = "motifs", 
                   n_samples = 3)

# Create amplitude envelope heatmap
sap |> plot_heatmap(balanced = TRUE)
```

## Feature Extraction and Analysis

``` r
# Extract spectral features
sap <- sap |>
  analyze_spectral(balanced = TRUE) |>
  find_clusters() |>
  run_umap()

# Visualize UMAP by time point
sap |> plot_umap(split.by = "label")
```

## Key Parameters for Bulk Processing

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
