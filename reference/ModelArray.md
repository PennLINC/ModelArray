# Construct a ModelArray object

Load element-wise data from an .h5 file as a \`ModelArray\` object.

## Usage

``` r
ModelArray(filepath, scalar_types = c("FD"), analysis_names = character(0))
```

## Arguments

- filepath:

  Path to an .h5 file

- scalar_types:

  Expected scalars

- analysis_names:

  The subfolder names for results in the .h5 file. If empty (default),
  results are not read.

## Value

A \`ModelArray\` object
