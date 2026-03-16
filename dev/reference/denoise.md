# Denoise Audio Files

Removes stationary background noise from audio recordings using
spectral-domain techniques. Two methods are available:

- **Spectral Median Subtraction** (`"spectral_median"`): Estimates the
  per-frequency noise floor as the quantile magnitude across all time
  frames, subtracts it from every frame (half-wave rectified), then
  reconstructs the waveform via inverse STFT. Fast and reliable gold
  standard for flat, stationary broadband noise.

- **Spectral Gating** (`"spectral_gate"`): Audacity-style noise
  reduction. Builds a per-frequency noise profile (floor + spread), then
  applies a smooth sigmoid gate so that time-frequency bins
  significantly above the noise floor pass through nearly unmodified,
  while bins close to the noise floor are gently attenuated. A
  configurable floor (`gate_floor`) prevents any bin from being
  completely silenced, avoiding the spectral "twist" artefacts that hard
  gating can introduce in bird-song recordings.

## Usage

``` r
denoise(x, ...)

# Default S3 method
denoise(
  x,
  method = c("spectral_median", "spectral_gate"),
  output_dir = NULL,
  overwrite = FALSE,
  wl = 256L,
  ovlp = 50L,
  wn = "hanning",
  plot = TRUE,
  view_window = NULL,
  freq_range = c(0, 10),
  noise_quantile = NULL,
  gain = 1,
  gate_threshold = 1.5,
  gate_smoothing = 3L,
  gate_floor = 0.1,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
denoise(
  x,
  method = c("spectral_median", "spectral_gate"),
  output_dir = NULL,
  overwrite = FALSE,
  wl = 256L,
  ovlp = 50L,
  wn = "hanning",
  plot = FALSE,
  view_window = NULL,
  freq_range = c(0, 10),
  noise_quantile = NULL,
  gain = 1,
  gate_threshold = 1.5,
  gate_smoothing = 3L,
  gate_floor = 0.1,
  day = NULL,
  indices = NULL,
  cores = NULL,
  update_base_path = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to process: a WAV file path (character) or a SAP object.

- ...:

  Additional arguments (currently unused).

- method:

  Character, the denoising method to use: `"spectral_median"` (default)
  or `"spectral_gate"`.

- output_dir:

  Character, directory for denoised WAV files. If `NULL` (default), a
  `denoised/` subdirectory is created beside each source file (default
  method) or under `x$base_path` (SAP method).

- overwrite:

  Logical, whether to overwrite existing denoised files (default:
  `FALSE`).

- wl:

  Integer, STFT window length in samples (default: 256).

- ovlp:

  Integer, STFT overlap in percent (default: 50).

- wn:

  Character, window function passed to
  [`seewave::ftwindow`](https://rdrr.io/pkg/seewave/man/ftwindow.html)
  (default: `"hanning"`).

- plot:

  Logical, for default method only. When `TRUE` (default), plots the
  spectrogram of the denoised output audio using the same visualization
  style as
  [`segment`](https://lxiao06.github.io/ASAP/dev/reference/segment.md).

- view_window:

  Numeric vector of length 2 (seconds) for the plot window in the
  denoised spectrogram. `NULL` (default) shows the full file. Example:
  `c(1, 4)` plots 1s to 4s.

- freq_range:

  Numeric vector of length 2 giving the frequency range (in kHz) to
  denoise, e.g. `c(0, 10)` (default). Only frequency bins within this
  range are passed through the denoising algorithm; bins outside the
  range are left untouched. Set to `NULL` to denoise the entire
  spectrum. Restricting the range speeds up processing and avoids
  altering frequency bands that contain no noise.

- noise_quantile:

  Numeric in (0, 1\] or `NULL`. The role differs by method
  (auto-selected when `NULL`):

  - **spectral_median**: quantile of the per-frequency magnitude
    distribution across *all* time frames used as the noise floor
    estimate. Auto-default: `0.5` (median).

  - **spectral_gate**: fraction of time frames (ranked by total frame
    energy, lowest first) selected as *quiet/noise-only* frames. The
    noise floor and spread are computed exclusively from those frames,
    so song syllables never contaminate the noise profile. Auto-default:
    `0.25` (lowest-energy 25% of frames).

  Lower values are more conservative; higher values capture more of the
  noise distribution.

- gain:

  Numeric \\\ge 0\\, over-subtraction factor for **spectral_median**
  (default: 1.0). Values \> 1 increase noise removal at the risk of
  artefacts.

- gate_threshold:

  Numeric \\\ge 0\\, for **spectral_gate**: number of noise-spread
  standard deviations above the noise floor at which the sigmoid gate
  reaches 50% transmission (default: 1.5). Lower = less aggressive (more
  signal preserved); higher = more aggressive (more noise removed).

- gate_smoothing:

  Integer, for **spectral_gate**: half-width (in frequency bins) of the
  box-car smoothing applied to the gate mask (default: 3). Smoothing
  reduces sharp spectral edges that can distort song structure. Set to 0
  to disable.

- gate_floor:

  Numeric in \[0, 1), for **spectral_gate**: minimum gate value applied
  to every bin (default: 0.1). A non-zero floor ensures no frequency bin
  is completely silenced, preserving tonal continuity and preventing the
  "hollow" or "twisted" sound artefacts of hard gating.

- verbose:

  Logical, print progress messages (default: `TRUE`).

- day:

  For SAP objects: Numeric vector of days post-hatch to process. `NULL`
  (default) = all days.

- indices:

  For SAP objects: Numeric vector of row indices in `x$metadata`. `NULL`
  = all rows within selected days.

- cores:

  For SAP objects: Number of parallel cores (default:
  `parallel::detectCores() - 1`).

- update_base_path:

  Logical, for SAP objects only. When `TRUE`, `x$base_path` is replaced
  with the denoised output directory after processing completes
  (default: `FALSE`). This makes all subsequent pipeline steps
  (`detect_template`, `export_clips`, `extract_envelope`, etc.)
  automatically read from the denoised files without any further
  configuration, since the denoised directory mirrors the original
  `day_post_hatch/filename` structure exactly.

## Value

- Default method:

  Character path to the denoised WAV file (invisibly).

- SAP method:

  Updated SAP object with `x$denoised_path` set to the output directory,
  mirroring the `day_post_hatch/filename` structure of the originals.

## Details

**Spectral Median Subtraction** steps:

1.  Compute STFT magnitude and phase.

2.  Estimate the per-frequency noise floor as the `noise_quantile`
    quantile of magnitude across all time frames.

3.  Subtract `gain * noise_floor` from every frame; zero out negatives
    (half-wave rectification).

4.  Reconstruct the waveform from cleaned magnitude + original phase via
    inverse STFT (overlap-add).

**Spectral Gating** steps:

1.  Compute STFT magnitude and phase.

2.  Rank all time frames by their total energy (\\E_t = \sum_f
    \|S(t,f)\|^2\\). Select the lowest-energy `noise_quantile` fraction
    as *quiet frames*. Because bird-song syllables are brief and
    high-energy, they are excluded from this set, so the noise profile
    is estimated from background-only frames only.

3.  Compute the per-frequency noise floor (\\\mu_f\\) as the mean and
    spread (\\\sigma_f\\) as the SD of those quiet frames.

4.  Build a soft sigmoid gate for every time-frequency bin: \$\$G(t,f) =
    \max\\\left(g\_{\min},\\ \sigma\\\left(\frac{\|S(t,f)\| - \mu_f -
    \theta\sigma_f}{\sigma_f/4}\right)\right)\$\$ where \\\theta\\ =
    `gate_threshold` and \\g\_{\min}\\ = `gate_floor`.

5.  Optionally smooth the gate mask along the frequency axis
    (`gate_smoothing`) to reduce spectral edge artefacts.

6.  Multiply the gate mask by the original magnitude and reconstruct the
    waveform.

## See also

[`detect_template`](https://lxiao06.github.io/ASAP/dev/reference/detect_template.md)
for template matching after denoising.

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Single WAV file -------------------------------------------------------
# Spectral median (default — fast, reliable)
clean_path <- denoise("path/to/recording.wav")

# Spectral gate — cleaner background with auto noise_quantile = 0.25
clean_path <- denoise("path/to/recording.wav",
  method = "spectral_gate"
)

# Override gate aggressiveness
clean_path <- denoise("path/to/recording.wav",
  method         = "spectral_gate",
  gate_threshold = 2.0, # more aggressive
  gate_floor     = 0.05, # tighter floor
  gate_smoothing = 5L # extra smoothing
)

# Spectral median — auto noise_quantile = 0.5
clean_path <- denoise("path/to/recording.wav",
  method = "spectral_median",
  gain   = 1.3
)

# --- SAP object (batch) ---------------------------------------------------
sap_obj <- denoise(sap_obj, method = "spectral_median", cores = 4)

sap_obj <- denoise(sap_obj,
  method         = "spectral_gate",
  day            = c(70, 75, 80),
  gate_threshold = 1.0,
  gate_floor     = 0.1,
  cores          = 8
)
} # }
```
