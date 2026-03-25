# Restore C-level stderr output

Restores file descriptor 2 from a previously saved descriptor returned
by
[`suppress_stderr()`](https://lxiao06.github.io/ASAP/dev/reference/suppress_stderr.md).

## Usage

``` r
restore_stderr(saved_fd)
```

## Arguments

- saved_fd:

  Integer file descriptor from
  [`suppress_stderr()`](https://lxiao06.github.io/ASAP/dev/reference/suppress_stderr.md).
