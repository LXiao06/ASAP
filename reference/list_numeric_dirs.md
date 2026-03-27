# List Numeric Subdirectory Names

Detects all subfolders within a specified directory that have purely
numeric names (e.g., "1", "100", "042") and returns them as a character
vector.

## Usage

``` r
list_numeric_dirs(directory)
```

## Arguments

- directory:

  Path to the parent directory to search within

## Value

A character vector of the numeric subfolder names. Returns an empty
character vector if none are found.
