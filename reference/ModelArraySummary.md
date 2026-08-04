# Summarize ModelArray storage without loading a full ModelArray

Reads the storage structure and returns a summary of available scalars,
their dimensions, and any saved analyses. Useful for inspecting large
files without constructing a full
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
object.

## Usage

``` r
ModelArraySummary(filepath, backend = c("auto", "hdf5", "tiledb"))

h5summary(filepath, backend = c("auto", "hdf5", "tiledb"))

# S3 method for class 'ModelArraySummary'
print(x, ...)

# S3 method for class 'h5summary'
print(x, ...)
```

## Arguments

- filepath:

  Character. Path to an HDF5 (`.h5`) file or TileDB (`.tdb`) store.

- backend:

  Character. Storage backend: `"auto"` (default), `"hdf5"`, or
  `"tiledb"`.

- x:

  A `ModelArraySummary` or `h5summary` object.

- ...:

  Additional arguments (currently ignored).

## Value

`ModelArraySummary()` returns an object of class
`c("ModelArraySummary", "h5summary")`. `h5summary()` returns an object
of class `"h5summary"` for backward compatibility. Both are lists with
components:

- scalars:

  A data.frame with columns `name`, `nElements`, and `nInputFiles`.

- analyses:

  Character vector of analysis names found under `/results/`.

- filepath:

  The input filepath.

- backend:

  The resolved storage backend.

Invisible `x`. Called for its side effect of printing a human-readable
summary to the console.

## Details

For HDF5, this function opens the file read-only via
[`h5ls`](https://huber-group-embl.github.io/rhdf5/reference/h5ls.html).
For TileDB, it inspects arrays under `/scalars/` and `/results/`. It
does not load full data matrices into memory. The returned object has a
`print` method that displays a formatted summary.

`h5summary()` is a backward-compatible alias for `ModelArraySummary()`.
Despite the historical name, it supports both HDF5 files and TileDB
stores.

## See also

[`ModelArray`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html)
for loading the full object,
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
for the class definition.

## Examples

``` r
if (FALSE) { # \dontrun{
ModelArraySummary("path/to/data.h5")
ModelArraySummary("path/to/store.tdb")

# Inspect before deciding which scalars to load
info <- ModelArraySummary("path/to/data.h5")
info$scalars$name
ma <- ModelArray("path/to/data.h5", scalar_types = info$scalars$name)

# Historical alias, still supported
h5summary("path/to/data.h5")
} # }
```
