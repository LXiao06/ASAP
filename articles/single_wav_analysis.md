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

# View detected bouts as a table
knitr::kable(bouts, digits = 3)
```

| filename       | selec | start_time | end_time |
|:---------------|------:|-----------:|---------:|
| zf_example.wav |     1 |      1.057 |    3.553 |
| zf_example.wav |     2 |      4.156 |    4.946 |

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

# View detected syllables as a table
knitr::kable(syllables, digits = 2)
```

| filename       | selec | threshold | .start | .end | start_time | end_time | duration | silence_gap |
|:---------------|------:|----------:|-------:|-----:|-----------:|---------:|---------:|------------:|
| zf_example.wav |     1 |        10 |   0.06 | 0.12 |       1.06 |     1.12 |     0.06 |          NA |
| zf_example.wav |     2 |        10 |   0.15 | 0.20 |       1.15 |     1.20 |     0.05 |        0.03 |
| zf_example.wav |     3 |        10 |   0.26 | 0.31 |       1.26 |     1.31 |     0.05 |        0.06 |
| zf_example.wav |     4 |        10 |   0.36 | 0.41 |       1.36 |     1.41 |     0.05 |        0.05 |
| zf_example.wav |     5 |        10 |   0.44 | 0.52 |       1.44 |     1.52 |     0.08 |        0.03 |
| zf_example.wav |     6 |        10 |   0.55 | 0.61 |       1.55 |     1.61 |     0.06 |        0.03 |
| zf_example.wav |     7 |        10 |   0.62 | 0.84 |       1.62 |     1.84 |     0.21 |        0.01 |
| zf_example.wav |     8 |        10 |   0.87 | 0.98 |       1.87 |     1.98 |     0.11 |        0.04 |
| zf_example.wav |     9 |        10 |   0.99 | 1.05 |       1.99 |     2.05 |     0.06 |        0.01 |
| zf_example.wav |    10 |        10 |   1.09 | 1.22 |       2.09 |     2.22 |     0.12 |        0.04 |
| zf_example.wav |    11 |        10 |   1.33 | 1.54 |       2.33 |     2.54 |     0.20 |        0.12 |
| zf_example.wav |    12 |        10 |   1.67 | 1.74 |       2.67 |     2.74 |     0.07 |        0.13 |
| zf_example.wav |    13 |        10 |   1.75 | 1.84 |       2.75 |     2.84 |     0.09 |        0.02 |
| zf_example.wav |    14 |        10 |   1.86 | 1.92 |       2.86 |     2.92 |     0.06 |        0.02 |
| zf_example.wav |    15 |        10 |   1.93 | 2.16 |       2.93 |     3.16 |     0.22 |        0.01 |
| zf_example.wav |    16 |        10 |   2.18 | 2.29 |       3.18 |     3.29 |     0.11 |        0.03 |
| zf_example.wav |    17 |        10 |   2.31 | 2.36 |       3.31 |     3.36 |     0.05 |        0.01 |
| zf_example.wav |    18 |        10 |   2.40 | 2.54 |       3.40 |     3.54 |     0.14 |        0.04 |
| zf_example.wav |    19 |        10 |   3.18 | 3.24 |       4.18 |     4.24 |     0.06 |        0.63 |
| zf_example.wav |    20 |        10 |   3.27 | 3.34 |       4.27 |     4.34 |     0.07 |        0.03 |
| zf_example.wav |    21 |        10 |   3.35 | 3.44 |       4.35 |     4.44 |     0.09 |        0.01 |
| zf_example.wav |    22 |        10 |   3.46 | 3.51 |       4.46 |     4.51 |     0.06 |        0.02 |
| zf_example.wav |    23 |        10 |   3.53 | 3.75 |       4.53 |     4.75 |     0.23 |        0.01 |
| zf_example.wav |    24 |        10 |   3.78 | 3.88 |       4.78 |     4.88 |     0.10 |        0.02 |
| zf_example.wav |    25 |        10 |   3.90 | 3.95 |       4.90 |     4.95 |     0.05 |        0.02 |

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
# Extract amplitude envelope for the second bout
if (!is.null(bouts) && nrow(bouts) >= 2) {
  env_bout <- amp_env(bouts[2, ], 
                      wav_dir = dirname(wav_file),
                      msmooth = c(256, 50),
                      norm = TRUE,
                      plot = TRUE)
}
```

![](single_wav_analysis_files/figure-html/amplitude-envelope-bout-1.png)

We can also extract the envelope for individual syllables:

``` r
# Extract amplitude envelope for the 7th syllable
if (!is.null(syllables) && nrow(syllables) >= 7) {
  env_syl <- amp_env(syllables[7, ], 
                     wav_dir = dirname(wav_file),
                     msmooth = c(256, 50),
                     norm = TRUE,
                     plot = TRUE)
}
```

![](single_wav_analysis_files/figure-html/amplitude-envelope-syllable-1.png)

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

# 3. Segment syllables
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

# 5. Extract amplitude envelope
env <- amp_env(bouts[1, ], 
               wav_dir = dirname(wav_file),
               msmooth = c(256, 50),
               norm = TRUE,
               plot = TRUE)
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
#> [29] dplyr_1.1.4        ggplot2_4.0.1      reticulate_1.44.1  vctrs_0.7.1       
#> [33] R6_2.6.1           png_0.1-8          lifecycle_1.0.5    fs_1.6.6          
#> [37] MASS_7.3-65        ragg_1.5.0         pkgconfig_2.0.3    desc_1.4.3        
#> [41] pkgdown_2.2.0      pillar_1.11.1      bslib_0.9.0        gtable_0.3.6      
#> [45] glue_1.8.0         Rcpp_1.1.1         systemfonts_1.3.1  xfun_0.56         
#> [49] tibble_3.3.1       tidyselect_1.2.1   knitr_1.51         farver_2.1.2      
#> [53] htmltools_0.5.9    patchwork_1.3.2    rmarkdown_2.30     signal_1.8-1      
#> [57] compiler_4.5.2     S7_0.2.1
```
