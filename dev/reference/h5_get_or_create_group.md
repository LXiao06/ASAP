# Get or Create an HDF5 Group

Returns an existing group from an HDF5 parent object, or creates it if
it does not already exist.

## Usage

``` r
h5_get_or_create_group(parent, group_name)
```

## Arguments

- parent:

  An `H5File` or `H5Group` object.

- group_name:

  Character name of the group to retrieve or create.

## Value

The requested HDF5 group object.
