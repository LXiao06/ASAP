# Normalize Spectrogram

Internal function to normalize spectrogram using different methods.

## Usage

``` r
normalize_spec(spec, max_level_db = NULL, ref_level_db = NULL, method = "db")
```

## Arguments

- spec:

  Input spectrogram matrix

- max_level_db:

  Maximum threshold level

- ref_level_db:

  Reference level

- method:

  Normalization method

## Value

Normalized spectrogram matrix

## Details

Supports multiple normalization methods:

- db: dB-scale normalization

- minmax: Min-max scaling

- zscore: Z-score normalization
