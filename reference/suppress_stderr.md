# Suppress C-level stderr output

Redirects file descriptor 2 (C stderr) to /dev/null and returns the
saved file descriptor so it can be restored later. This suppresses
C-library warnings (e.g. fontconfig's "using without calling FcInit()")
that R's sink() cannot intercept.

## Usage

``` r
suppress_stderr()
```

## Value

Integer file descriptor to pass to
[`restore_stderr()`](https://lxiao06.github.io/ASAP/reference/restore_stderr.md),
or -1 on Windows where this is a no-op.
