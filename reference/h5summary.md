# Summarize an h5 file without loading the full ModelArray

Reads the h5 file structure and returns a summary of available scalars,
their dimensions, and any saved analyses. Useful for inspecting large
files without constructing a full ModelArray object.

## Usage

``` r
h5summary(filepath)
```

## Arguments

- filepath:

  Path to an h5 file

## Value

A list with components:

- scalars:

  A data.frame with columns: name, nElements, nInputFiles

- analyses:

  Character vector of analysis names

- filepath:

  The input filepath
