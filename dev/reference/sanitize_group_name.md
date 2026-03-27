# Sanitize a String for Use as a Group Name

Replaces any character that is not alphanumeric, a dot, underscore, or
hyphen with an underscore. Returns `"unknown"` for empty or NA input.

## Usage

``` r
sanitize_group_name(x)
```

## Arguments

- x:

  Character scalar to sanitize.

## Value

Sanitized character string.
