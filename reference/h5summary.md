# Summarize an HDF5 file without loading a full ModelArray

Reads the HDF5 file structure and returns a summary of available
scalars, their dimensions, and any saved analyses. Useful for inspecting
large files without constructing a full
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
object.

## Usage

``` r
h5summary(filepath)

# S3 method for class 'h5summary'
print(x, ...)
```

## Arguments

- filepath:

  Character. Path to an HDF5 (`.h5`) file.

- x:

  An `h5summary` object as returned by `h5summary`.

- ...:

  Additional arguments (currently ignored).

## Value

An object of class `"h5summary"`, which is a list with components:

- scalars:

  A data.frame with columns `name`, `nElements`, and `nInputFiles`.

- analyses:

  Character vector of analysis names found under `/results/`.

- filepath:

  The input filepath.

Invisible `x`. Called for its side effect of printing a human-readable
summary to the console.

## Details

This function opens the HDF5 file read-only via
[`h5ls`](https://huber-group-embl.github.io/rhdf5/reference/h5ls.html),
inspects the group structure under `/scalars/` and `/results/`, and
closes the file. It does not load any data into memory. The returned
object has a `print` method that displays a formatted summary.

## See also

[`ModelArray`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html)
for loading the full object,
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
for the class definition.

## Examples

``` r
if (FALSE) { # \dontrun{
h5summary("path/to/data.h5")

# Inspect before deciding which scalars to load
info <- h5summary("path/to/data.h5")
info$scalars$name
ma <- ModelArray("path/to/data.h5", scalar_types = info$scalars$name)
} # }
```
