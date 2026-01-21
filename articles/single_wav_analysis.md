# Single WAV File Analysis with ASAP

## Introduction

This vignette demonstrates how to use ASAP functions for analyzing
single WAV files of zebra finch vocalizations. While ASAP is designed
for large-scale longitudinal studies, all core analysis functions work
directly with individual audio files, making it easy to explore and
understand your data.

We’ll cover:

1.  **Audio Visualization** - View spectrograms
2.  **Bout Detection** - Find singing periods
3.  **Syllable Segmentation** - Detect individual syllables
4.  **Spectral Entropy** - Measure spectral structure
5.  **Fundamental Frequency** - Extract pitch contours
6.  **Amplitude Envelope** - Analyze temporal dynamics

## Setup

``` r
library(ASAP)

# Get path to example WAV file included with the package
wav_file <- system.file("extdata", "zf_example.wav", package = "ASAP")
```

## 1. Audio Visualization

The
[`visualize_song()`](https://lxiao06.github.io/ASAP/reference/visualize_song.md)
function creates spectrogram visualizations of audio recordings. This is
typically the first step in exploring your data.

### Full recording

``` r
# Visualize the entire recording
visualize_song(wav_file)
```

![](single_wav_analysis_files/figure-html/visualize-full-1.png)

    #> Song visualization completed for: zf_example.wav

### Specific time window

You can focus on a specific time range using `start_time_in_second` and
`end_time_in_second`:

``` r
# Visualize a 3-second segment
visualize_song(wav_file, 
               start_time_in_second = 1, 
               end_time_in_second = 4)
```

![](single_wav_analysis_files/figure-html/visualize-segment-1.png)

    #> Song visualization completed for: zf_example.wav

## 2. Bout Detection

A “bout” is a continuous period of singing. The
[`find_bout()`](https://lxiao06.github.io/ASAP/reference/find_bout.md)
function automatically detects bout boundaries using RMS (root mean
square) amplitude thresholding with bandpass filtering.

``` r
# Detect bouts in the recording
bouts <- find_bout(wav_file, 
                   rms_threshold = 0.1,    # Amplitude threshold (0-1)
                   min_duration = 0.7,      # Minimum bout length in seconds
                   plot = TRUE)             # Show detection plot
```

![](single_wav_analysis_files/figure-html/find-bout-1.png)

``` r

# View detected bouts
bouts
#>         filename selec start_time end_time
#> 1 zf_example.wav     1   1.056508 3.552653
#> 2 zf_example.wav     2   4.156372 4.945850
```

### Key parameters

- `rms_threshold`: Higher values require louder sounds; lower values
  detect quieter vocalizations
- `min_duration`: Minimum bout length (filters out short sounds/noise)
- `freq_range`: Bandpass filter range (default: 1-8 kHz for zebra finch)

## 3. Syllable Segmentation

The [`segment()`](https://lxiao06.github.io/ASAP/reference/segment.md)
function detects individual syllables within a specified time window
using dynamic spectral thresholding.

``` r
# Segment syllables in a time window
syllables <- segment(wav_file, 
                     start_time = 1,          # Start time (seconds)
                     end_time = 5,            # End time (seconds)
                     flim = c(1, 8),          # Frequency limits (kHz)
                     silence_threshold = 0.01,
                     min_syllable_ms = 20,    # Minimum syllable duration
                     max_syllable_ms = 240,   # Maximum syllable duration
                     min_level_db = 10,       # Starting threshold (dB)
                     verbose = FALSE)
```

![](single_wav_analysis_files/figure-html/segment-1.png)

``` r

# View detected syllables
syllables
#>          filename selec threshold     .start      .end start_time end_time
#> 1  zf_example.wav     1        10 0.06146572 0.1229314   1.061466 1.122931
#> 2  zf_example.wav     2        10 0.15130024 0.1985816   1.151300 1.198582
#> 3  zf_example.wav     3        10 0.26004728 0.3120567   1.260047 1.312057
#> 4  zf_example.wav     4        10 0.35933806 0.4113475   1.359338 1.411348
#> 5  zf_example.wav     5        10 0.44444444 0.5200946   1.444444 1.520095
#> 6  zf_example.wav     6        10 0.54846336 0.6099291   1.548463 1.609929
#> 7  zf_example.wav     7        10 0.62411348 0.8368794   1.624113 1.836879
#> 8  zf_example.wav     8        10 0.87470449 0.9834515   1.874704 1.983452
#> 9  zf_example.wav     9        10 0.99290780 1.0496454   1.992908 2.049645
#> 10 zf_example.wav    10        10 1.09219858 1.2151300   2.092199 2.215130
#> 11 zf_example.wav    11        10 1.33333333 1.5366430   2.333333 2.536643
#> 12 zf_example.wav    12        10 1.66903073 1.7352246   2.669031 2.735225
#> 13 zf_example.wav    13        10 1.75413712 1.8392435   2.754137 2.839243
#> 14 zf_example.wav    14        10 1.85815603 1.9196217   2.858156 2.919622
#> 15 zf_example.wav    15        10 1.93380615 2.1560284   2.933806 3.156028
#> 16 zf_example.wav    16        10 2.18439716 2.2931442   3.184397 3.293144
#> 17 zf_example.wav    17        10 2.30732861 2.3593381   3.307329 3.359338
#> 18 zf_example.wav    18        10 2.40189125 2.5437352   3.401891 3.543735
#> 19 zf_example.wav    19        10 3.17730496 3.2387707   4.177305 4.238771
#> 20 zf_example.wav    20        10 3.27186761 3.3427896   4.271868 4.342790
#> 21 zf_example.wav    21        10 3.35224586 3.4373522   4.352246 4.437352
#> 22 zf_example.wav    22        10 3.45626478 3.5130024   4.456265 4.513002
#> 23 zf_example.wav    23        10 3.52718676 3.7541371   4.527187 4.754137
#> 24 zf_example.wav    24        10 3.77777778 3.8817967   4.777778 4.881797
#> 25 zf_example.wav    25        10 3.90070922 3.9479905   4.900709 4.947991
#>      duration silence_gap
#> 1  0.06146572          NA
#> 2  0.04728132 0.028368794
#> 3  0.05200946 0.061465721
#> 4  0.05200946 0.047281324
#> 5  0.07565012 0.033096927
#> 6  0.06146572 0.028368794
#> 7  0.21276596 0.014184397
#> 8  0.10874704 0.037825059
#> 9  0.05673759 0.009456265
#> 10 0.12293144 0.042553191
#> 11 0.20330969 0.118203310
#> 12 0.06619385 0.132387707
#> 13 0.08510638 0.018912530
#> 14 0.06146572 0.018912530
#> 15 0.22222222 0.014184397
#> 16 0.10874704 0.028368794
#> 17 0.05200946 0.014184397
#> 18 0.14184397 0.042553191
#> 19 0.06146572 0.633569740
#> 20 0.07092199 0.033096927
#> 21 0.08510638 0.009456265
#> 22 0.05673759 0.018912530
#> 23 0.22695035 0.014184397
#> 24 0.10401891 0.023640662
#> 25 0.04728132 0.018912530
```

### Understanding the output

The returned data frame contains:

- `start_time`/`end_time`: Syllable boundaries (seconds)
- `duration`: Syllable length
- `silence_gap`: Gap before the next syllable (useful for syntax
  analysis)
- `selec`: Selection number for tracking

## 4. Spectral Entropy Analysis

Spectral entropy measures the “randomness” or structure in the frequency
distribution. Harmonic sounds (like syllables) have low entropy, while
noisy sounds have high entropy.

``` r
# Calculate spectral entropy for a segment
entropy_result <- spectral_entropy(wav_file,
                                   start_time = 1.5,
                                   end_time = 2.5,
                                   method = "wiener",  # or "shannon"
                                   normalize = TRUE,
                                   plot = TRUE)
```

![](single_wav_analysis_files/figure-html/spectral-entropy-1.png)

### Interpreting entropy values

- **Low entropy (near 0)**: Structured, harmonic sounds (e.g., tonal
  syllables)
- **High entropy (near 1)**: Noisy, unstructured sounds (e.g., calls,
  noise)

### Available methods

- `"wiener"`: Wiener entropy (ratio of geometric to arithmetic mean)
- `"shannon"`: Shannon entropy (information-theoretic measure)

## 5. Fundamental Frequency (Pitch) Analysis

The
[`FF()`](https://lxiao06.github.io/ASAP/reference/Fundamental_Frequency.md)
function extracts the fundamental frequency (F0) contour, which
represents the perceived pitch of the vocalization over time.

``` r
# Extract fundamental frequency
pitch_result <- FF(wav_file,
                   start_time = 1.5,
                   end_time = 2.5,
                   method = "cepstrum",  # Cepstral analysis
                   fmax = 1400,          # Maximum F0 to detect (Hz)
                   threshold = 10,       # Confidence threshold
                   plot = TRUE)
```

![](single_wav_analysis_files/figure-html/fundamental-frequency-1.png)

### Key parameters

- `fmax`: Maximum fundamental frequency to detect (1400 Hz typical for
  zebra finch)
- `threshold`: Higher values filter out uncertain estimates
- `method`: `"cepstrum"` (default) or `"yin"` (requires Python librosa)

### The result contains:

- `f0`: Fundamental frequency values over time
- `time`: Corresponding time stamps

## 6. Amplitude Envelope

The amplitude envelope represents the temporal dynamics of sound
intensity. Use
[`amp_env()`](https://lxiao06.github.io/ASAP/reference/amp_env.md) on
segmented data to extract envelope profiles.

``` r
# First, we need a segment data frame with proper structure
# Using the first syllable from our segmentation
if (!is.null(syllables) && nrow(syllables) > 0) {
  # Extract amplitude envelope for the first syllable
  env <- amp_env(syllables[1, ], 
                 wav_dir = dirname(wav_file),
                 msmooth = c(256, 50),  # Smoothing parameters
                 norm = TRUE,           # Normalize to 0-1
                 plot = TRUE)
}
```

![](single_wav_analysis_files/figure-html/amplitude-envelope-1.png)

### Smoothing parameters

The `msmooth` argument controls envelope smoothing:

- First value: Window length (samples)
- Second value: Overlap percentage

## Putting It All Together

Here’s a typical workflow for single-file analysis:

``` r
library(ASAP)

# 1. Load and visualize
wav_file <- "path/to/your/recording.wav"
visualize_song(wav_file, start_time_in_second = 0, end_time_in_second = 10)

# 2. Detect bouts
bouts <- find_bout(wav_file, rms_threshold = 0.1, min_duration = 0.7)

# 3. Segment syllables (focusing on first bout)
if (!is.null(bouts) && nrow(bouts) > 0) {
  syllables <- segment(wav_file,
                       start_time = bouts$start_time[1],
                       end_time = bouts$end_time[1],
                       min_level_db = 10)
}

# 4. Analyze acoustic features
entropy <- spectral_entropy(wav_file, 
                            start_time = bouts$start_time[1],
                            end_time = bouts$end_time[1])

pitch <- FF(wav_file,
            start_time = bouts$start_time[1],
            end_time = bouts$end_time[1])
```

## Next Steps

For analyzing multiple recordings across developmental time points, see
the main ASAP workflow using SAP objects:

``` r
# Create a SAP object for batch analysis
sap <- create_sap_object(
  base_path = "/path/to/recordings",
  subfolders_to_include = c("day1", "day2", "day3"),
  labels = c("Pre", "During", "Post")
)

# Then use the pipeline functions:
sap <- sap |>
  find_bout() |>
  segment() |>
  analyze_spectral() |>
  find_clusters() |>
  run_umap()
```

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
#> other attached packages:
#> [1] ASAP_0.3.3
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10        generics_0.1.4     tidyr_1.3.2        lattice_0.22-7    
#>  [5] digest_0.6.39      magrittr_2.0.4     evaluate_1.0.5     grid_4.5.2        
#>  [9] RColorBrewer_1.1-3 fastmap_1.2.0      jsonlite_2.0.0     Matrix_1.7-4      
#> [13] tuneR_1.4.7        purrr_1.2.1        scales_1.4.0       pbapply_1.7-4     
#> [17] textshaping_1.0.4  jquerylib_0.1.4    cli_3.6.5          rlang_1.1.7       
#> [21] pbmcapply_1.5.1    fftw_1.0-9         seewave_2.2.4      cachem_1.1.0      
#> [25] yaml_2.3.12        av_0.9.6           tools_4.5.2        parallel_4.5.2    
#> [29] dplyr_1.1.4        ggplot2_4.0.1      reticulate_1.44.1  vctrs_0.7.0       
#> [33] R6_2.6.1           png_0.1-8          lifecycle_1.0.5    fs_1.6.6          
#> [37] MASS_7.3-65        ragg_1.5.0         pkgconfig_2.0.3    desc_1.4.3        
#> [41] pkgdown_2.2.0      pillar_1.11.1      bslib_0.9.0        gtable_0.3.6      
#> [45] glue_1.8.0         Rcpp_1.1.1         systemfonts_1.3.1  xfun_0.56         
#> [49] tibble_3.3.1       tidyselect_1.2.1   knitr_1.51         farver_2.1.2      
#> [53] htmltools_0.5.9    patchwork_1.3.2    rmarkdown_2.30     signal_1.8-1      
#> [57] compiler_4.5.2     S7_0.2.1
```
