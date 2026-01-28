# Constructing SAP Objects

## Introduction

This vignette demonstrates how to construct and organize **SAP (Song
Analysis Pipeline) objects** for longitudinal analysis of zebra finch
vocalizations. SAP objects serve as the central data structure in ASAP,
managing recordings across multiple developmental time points.

**Prerequisites**: Before reading this vignette, we recommend
completing:

- [Overview: Basic Audio
  Analysis](https://lxiao06.github.io/ASAP/articles/single_wav_analysis.md) -
  Core ASAP functions
- [Motif
  Detection](https://lxiao06.github.io/ASAP/articles/motif_detection.md) -
  Template optimization workflow

## Organizing Your Recording Data

ASAP expects recordings to be organized in a specific folder structure:

    base_path/
    ├── day_post_hatch_1/     # e.g., "190" for 190 dph
    │   ├── song_001.wav
    │   ├── song_002.wav
    │   └── ...
    ├── day_post_hatch_2/     # e.g., "201" for 201 dph
    │   ├── song_001.wav
    │   └── ...
    └── day_post_hatch_3/     # e.g., "203" for 203 dph
        ├── song_001.wav
        └── ...

Each subfolder represents a developmental time point (days post hatch).

## Creating a SAP Object

``` r
library(ASAP)

# Create SAP object from organized recording folders
sap <- create_sap_object(
  base_path = "/path/to/recordings",
  subfolders_to_include = c("190", "201", "203"),
  labels = c("Baseline", "Post", "Recovery")
)

# View the structure
print(sap)
summary(sap)
```

**Example output:**

    SAP Object
    ===========
    Base path: /path/to/recordings
    Time points: 3 (190, 201, 203)
    Labels: Baseline, Post, Recovery
    Total files: 847
      - 190 (Baseline): 312 files
      - 201 (Post): 285 files
      - 203 (Recovery): 250 files

### Key Arguments

| Argument                | Description                                        |
|-------------------------|----------------------------------------------------|
| `base_path`             | Root directory containing recording folders        |
| `subfolders_to_include` | Vector of folder names (typically days post hatch) |
| `labels`                | Human-readable labels for each time point          |

## Exploring SAP Object Contents

``` r
# View metadata
head(sap$metadata)

# Check number of files per time point
table(sap$metadata$label)

# Visualize sample recordings
visualize_song(sap, n_samples = 4, random = TRUE)
```

![Sample spectrograms from SAP object across time
points.](figures/sap_visualize_song.png)

Sample spectrograms from SAP object across time points.

## SAP Object Structure

A SAP object contains:

| Component    | Description                                                                                                     |
|--------------|-----------------------------------------------------------------------------------------------------------------|
| `$base_path` | Root directory path                                                                                             |
| `$metadata`  | Data frame with file information                                                                                |
| `$templates` | Template storage (after [`create_template()`](https://lxiao06.github.io/ASAP/reference/create_template.md))     |
| `$motifs`    | Detected motifs (after [`find_motif()`](https://lxiao06.github.io/ASAP/reference/find_motif.md))                |
| `$bouts`     | Detected bouts (after [`find_bout()`](https://lxiao06.github.io/ASAP/reference/find_bout.md))                   |
| `$features`  | Extracted features (after [`analyze_spectral()`](https://lxiao06.github.io/ASAP/reference/analyze_spectral.md)) |

## Next Steps

Once you have created a SAP object, proceed to:

- [Longitudinal Motif
  Detection](https://lxiao06.github.io/ASAP/articles/longitudinal_motif_detection.md) -
  Apply templates across all recordings

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
