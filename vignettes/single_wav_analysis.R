## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10,
  fig.height = 5,
  out.width = "100%"
)


## ----setup--------------------------------------------------------------------
library(ASAP)

# Get path to example WAV file included with the package
wav_file <- system.file("extdata", "zf_example.wav", package = "ASAP")


## ----visualize-full, fig.height=4---------------------------------------------
# Visualize the entire recording
visualize_song(wav_file)


## ----visualize-segment, fig.height=4------------------------------------------
# Visualize a 3-second segment
visualize_song(wav_file, 
               start_time_in_second = 1, 
               end_time_in_second = 4)


## ----find-bout----------------------------------------------------------------
# Detect bouts in the recording
bouts <- find_bout(wav_file, 
                   rms_threshold = 0.1,    # Amplitude threshold (0-1)
                   min_duration = 0.7,      # Minimum bout length in seconds
                   plot = TRUE)             # Show detection plot

# View detected bouts
bouts


## ----segment, fig.height=6----------------------------------------------------
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

# View detected syllables
syllables


## ----spectral-entropy, fig.height=5-------------------------------------------
# Calculate spectral entropy for a segment
entropy_result <- spectral_entropy(wav_file,
                                   start_time = 1.5,
                                   end_time = 2.5,
                                   method = "wiener",  # or "shannon"
                                   normalize = TRUE,
                                   plot = TRUE)


## ----fundamental-frequency, fig.height=5--------------------------------------
# Extract fundamental frequency
pitch_result <- FF(wav_file,
                   start_time = 1.5,
                   end_time = 2.5,
                   method = "cepstrum",  # Cepstral analysis
                   fmax = 1400,          # Maximum F0 to detect (Hz)
                   threshold = 10,       # Confidence threshold
                   plot = TRUE)


## ----amplitude-envelope, fig.height=4-----------------------------------------
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


## ----workflow, eval=FALSE-----------------------------------------------------
## library(ASAP)
## 
## # 1. Load and visualize
## wav_file <- "path/to/your/recording.wav"
## visualize_song(wav_file, start_time_in_second = 0, end_time_in_second = 10)
## 
## # 2. Detect bouts
## bouts <- find_bout(wav_file, rms_threshold = 0.1, min_duration = 0.7)
## 
## # 3. Segment syllables (focusing on first bout)
## if (!is.null(bouts) && nrow(bouts) > 0) {
##   syllables <- segment(wav_file,
##                        start_time = bouts$start_time[1],
##                        end_time = bouts$end_time[1],
##                        min_level_db = 10)
## }
## 
## # 4. Analyze acoustic features
## entropy <- spectral_entropy(wav_file,
##                             start_time = bouts$start_time[1],
##                             end_time = bouts$end_time[1])
## 
## pitch <- FF(wav_file,
##             start_time = bouts$start_time[1],
##             end_time = bouts$end_time[1])


## ----next-steps, eval=FALSE---------------------------------------------------
## # Create a SAP object for batch analysis
## sap <- create_sap_object(
##   base_path = "/path/to/recordings",
##   subfolders_to_include = c("day1", "day2", "day3"),
##   labels = c("Pre", "During", "Post")
## )
## 
## # Then use the pipeline functions:
## sap <- sap |>
##   find_bout() |>
##   segment() |>
##   analyze_spectral() |>
##   find_clusters() |>
##   run_umap()


## ----session-info-------------------------------------------------------------
sessionInfo()

