# Element metadata from a ModelArray

Reads element metadata (e.g., greyordinates for cifti data) from the h5
file if present. Returns NULL if no element metadata is found.

## Usage

``` r
elementMetadata(x)

# S4 method for class 'ModelArray'
elementMetadata(x)
```

## Arguments

- x:

  A ModelArray object

## Value

A matrix or data.frame of element metadata, or NULL
