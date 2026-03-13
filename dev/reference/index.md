# Package index

## Object Construction

Functions for creating and writing the core SAP object

- [`create_sap_object()`](https://lxiao06.github.io/ASAP/dev/reference/create_sap_object.md)
  : Create a Sound Analysis Pro (SAP) Object from Audio Recordings
- [`compute_wav_durations()`](https://lxiao06.github.io/ASAP/dev/reference/compute_wav_durations.md)
  : Compute WAV File Durations
- [`create_sap_metadata()`](https://lxiao06.github.io/ASAP/dev/reference/create_sap_metadata.md)
  : Create metadata for audio files recorded by SAP2011 (Sound Analysis
  Pro)

## Curated Song Export

Functions for exporting audio clips of detected bouts and motifs.

- [`create_audio_clip()`](https://lxiao06.github.io/ASAP/dev/reference/create_audio_clip.md)
  : Create Audio Clips from Sound Files
- [`create_motif_clips()`](https://lxiao06.github.io/ASAP/dev/reference/create_motif_clips.md)
  : Create Motif Audio Clips
- [`create_bout_clips()`](https://lxiao06.github.io/ASAP/dev/reference/create_bout_clips.md)
  : Create Bout Audio Clips

## Detection & Segmentation

Functions for finding bouts, motives and segmenting syllables

- [`find_bout()`](https://lxiao06.github.io/ASAP/dev/reference/find_bout.md)
  : Detect Song Bouts in Audio Recordings
- [`find_motif()`](https://lxiao06.github.io/ASAP/dev/reference/find_motif.md)
  : Find Motifs in Song Data
- [`refine_motif_boundaries()`](https://lxiao06.github.io/ASAP/dev/reference/refine_motif_boundaries.md)
  : Refine Motif Boundaries Using Segment Alignment
- [`segment()`](https://lxiao06.github.io/ASAP/dev/reference/segment.md)
  : Segment Audio into Syllables
- [`as_segment()`](https://lxiao06.github.io/ASAP/dev/reference/as_segment.md)
  : Convert data to a segment object

## Template Matching

Functions for establishing and detecting specific acoustic templates

- [`create_template()`](https://lxiao06.github.io/ASAP/dev/reference/create_template.md)
  : Create Correlation Templates for Song Analysis
- [`detect_template()`](https://lxiao06.github.io/ASAP/dev/reference/detect_template.md)
  : Detect Templates in Song Data

## Feature Extraction & Clustering

Functions to extract acoustic features (FF, goodness, etc) and cluster
segments

- [`extract_spec()`](https://lxiao06.github.io/ASAP/dev/reference/extract_spec.md)
  : Extract and pad spectrograms from audio segments
- [`analyze_spectral()`](https://lxiao06.github.io/ASAP/dev/reference/analyze_spectral.md)
  : Analyze Spectral Features of Audio Segments
- [`FF()`](https://lxiao06.github.io/ASAP/dev/reference/Fundamental_Frequency.md)
  : Fundamental Frequency Analysis and Visualization
- [`goodness()`](https://lxiao06.github.io/ASAP/dev/reference/goodness.md)
  : Pitch Goodness Analysis
- [`spectral_entropy()`](https://lxiao06.github.io/ASAP/dev/reference/spectral_entropy.md)
  : Calculate Spectral Entropy for Audio Segments
- [`find_clusters()`](https://lxiao06.github.io/ASAP/dev/reference/find_clusters.md)
  : Find Clusters in Feature Data
- [`amp_env()`](https://lxiao06.github.io/ASAP/dev/reference/amp_env.md)
  : Calculate Amplitude Envelope for Audio Segment
- [`refine_FF()`](https://lxiao06.github.io/ASAP/dev/reference/refine_FF.md)
  : Refine Fundamental Frequency Detection
- [`refine_sh()`](https://lxiao06.github.io/ASAP/dev/reference/refine_sh.md)
  : Refine Spectral Entropy Using Temporal Template
- [`anova_analysis()`](https://lxiao06.github.io/ASAP/dev/reference/anova_analysis.md)
  : Perform ANOVA and Multiple Comparisons Analysis

## Dimensionality Reduction

Functions dedicated to reducing acoustic parameter dimensionality

- [`run_pca()`](https://lxiao06.github.io/ASAP/dev/reference/run_pca.md)
  : Run Principal Component Analysis
- [`run_umap()`](https://lxiao06.github.io/ASAP/dev/reference/run_umap.md)
  : Run UMAP Dimensionality Reduction
- [`create_trajectory_matrix()`](https://lxiao06.github.io/ASAP/dev/reference/create_trajectory_matrix.md)
  : Create Spectrogram Matrices for Song Trajectory Analysis

## Plotting & Visualization

Plotting functionality for exploring ASAP objects

- [`plot_traces()`](https://lxiao06.github.io/ASAP/dev/reference/plot_traces.md)
  : Plot Traces of Song Features
- [`plot_clusters()`](https://lxiao06.github.io/ASAP/dev/reference/plot_clusters.md)
  : Plot Clusters in Time Series Data
- [`plot_heatmap()`](https://lxiao06.github.io/ASAP/dev/reference/plot_heatmap.md)
  : Plot Heatmap of Amplitude Envelopes
- [`plot_umap()`](https://lxiao06.github.io/ASAP/dev/reference/plot_umap.md)
  : Plot UMAP Visualization
- [`plot_umap2()`](https://lxiao06.github.io/ASAP/dev/reference/plot_umap2.md)
  : Plot UMAP Visualization for Trajectory Analysis
- [`plot_motif_boundaries()`](https://lxiao06.github.io/ASAP/dev/reference/plot_motif_boundaries.md)
  : Visualize Motif Boundaries in Audio Data
- [`visualize_song()`](https://lxiao06.github.io/ASAP/dev/reference/visualize_song.md)
  : Visualize Song Data
- [`visualize_segments()`](https://lxiao06.github.io/ASAP/dev/reference/visualize_segments.md)
  : Visualize Song Segments

## Labeling

Functions for annotating clustered syllables

- [`auto_label()`](https://lxiao06.github.io/ASAP/dev/reference/auto_label.md)
  : Automatic Syllable Labeling in Song Motif
- [`manual_label()`](https://lxiao06.github.io/ASAP/dev/reference/manual_label.md)
  : Manual Syllable Labeling for Bird Song Analysis

## Segment Methods

S3 methods for dealing with segment objects

- [`print(`*`<segment>`*`)`](https://lxiao06.github.io/ASAP/dev/reference/print.segment.md)
  : Print method for segment objects
- [`summary(`*`<segment>`*`)`](https://lxiao06.github.io/ASAP/dev/reference/summary.segment.md)
  : Summary method for segment objects
