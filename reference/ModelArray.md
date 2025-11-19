# An S4 class to represent element-wise scalar data and statistics.

An S4 class to represent element-wise scalar data and statistics.

Load element-wise data from .h5 file as an ModelArray object

## Usage

``` r
ModelArray(filepath, scalar_types = c("FD"), analysis_names = character(0))

ModelArray(filepath, scalar_types = c("FD"), analysis_names = character(0))
```

## Arguments

- filepath:

  file

- scalar_types:

  expected scalars

- analysis_names:

  the subfolder names for results in .h5 file. If empty (default),
  results are not read.

## Value

ModelArray object

## Details

: Tips for debugging: if you run into this error: "Error in
h(simpleError(msg, call)) : error in evaluating the argument 'seed' in
selecting a method for function 'DelayedArray': HDF5. Symbol table.
Can't open object." Then please check if you give correct
"scalar_types" - check via rhdf5::h5ls(filename_for_h5)

## Slots

- `sources`:

  A list of source filenames

- `scalars`:

  A list of element-wise scalar matrix

- `results`:

  A list of statistical result matrix

- `path`:

  Path to the h5 file on disk
