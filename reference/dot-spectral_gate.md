# Spectral Gating (soft sigmoid gate with configurable floor)

Builds a per-frequency noise profile (floor \\\mu_f\\ and spread
\\\sigma_f\\) then computes a gate value for every time-frequency bin:
\$\$G(t,f) = \max\\\left(g\_{\min},\\ \frac{1}{1 + e^{-(\|S(t,f)\| -
\mu_f - \theta\sigma_f)\\/\\(\sigma_f/4)}} \right)\$\$ The floor
\\g\_{\min}\\ (`floor`) prevents complete silencing and is the key
parameter for preserving tonal continuity of bird song.

## Usage

``` r
.spectral_gate(
  mag,
  noise_quantile = 0.25,
  threshold = 1.5,
  smoothing = 3L,
  floor = 0.1
)
```

## Arguments

- mag:

  Numeric matrix (freq x time) - magnitude spectrogram.

- noise_quantile:

  Quantile for noise-floor estimation.

- threshold:

  Number of \\\sigma\\ above noise floor where gate reaches 50%
  transmission (default: 1.5).

- smoothing:

  Half-width (bins) of frequency-axis box-car mask smoother; 0 =
  disabled (default: 3).

- floor:

  Minimum gate value in \[0, 1).

## Value

Cleaned magnitude matrix.
