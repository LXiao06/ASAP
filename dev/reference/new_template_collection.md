# Internal Constructor for Template Collection

An internal constructor function for creating a new template collection
with a standardized structure for storing and managing sound analysis
templates.

## Usage

``` r
new_template_collection()
```

## Value

A template_collection object with an empty but structured template
collection

## Details

Creates a template collection object with three key components:

- template_info: A data frame containing metadata about templates

- template_list: A list to store corTemplate objects

- template_matches: A list to store template detection results

The template_info data frame includes the following columns:

- template_name: Name of the template

- start_time: Start time of the template

- end_time: End time of the template

- duration: Duration of the template

- freq_min: Minimum frequency

- freq_max: Maximum frequency

- threshold: Detection threshold

- clip_name: Name of the audio clip

- clip_path: Path to the audio clip

- source_file: Original source file name

- source_file_path: Path to the source file

- creation_date: Date and time of template creation

## Note

This is an internal function not intended to be called directly by
users. It is used within the package for creating template collections.
