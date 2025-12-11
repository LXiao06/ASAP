# Split Template at Local Minima in Goodness

Split Template at Local Minima in Goodness

## Usage

``` r
split_at_dips(template, goodness, min_duration_samples)
```

## Arguments

- template:

  Logical vector representing the template

- goodness:

  Numeric vector of goodness values

- min_duration_samples:

  Minimum number of samples for a valid segment

## Value

Logical vector with segments split at significant dips
