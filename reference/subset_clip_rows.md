# Subset Clip Rows by Index or Random Sample

Subsets a clip data frame using explicit row indices or by randomly
sampling up to `n_clips` rows per day. Attaches a `.source_index` column
recording each row's original position in the input.

## Usage

``` r
subset_clip_rows(rows, indices = NULL, n_clips = NULL, seed = 222)
```

## Arguments

- rows:

  Data frame of clip rows (e.g. motifs or bouts).

- indices:

  Optional integer vector of explicit row indices.

- n_clips:

  Optional integer. Maximum clips to sample per day when `indices` is
  NULL.

- seed:

  Integer random seed (fixed to 222 in this repository).

## Value

Subset of `rows` with an added `.source_index` column.
