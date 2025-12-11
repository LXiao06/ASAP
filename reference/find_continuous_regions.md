# Find Continuous Regions

Internal function to identify continuous regions in a binary signal.

## Usage

``` r
find_continuous_regions(signal)
```

## Arguments

- signal:

  Binary signal vector

## Value

Matrix of start and end indices

## Details

Processes binary signal to:

- Find transitions between states

- Identify continuous regions

- Handle edge cases
