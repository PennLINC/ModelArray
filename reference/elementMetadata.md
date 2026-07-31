# Element metadata from a ModelArray

Reads element metadata (e.g., greyordinates for CIFTI data, or
fixel/voxel coordinate information) from the HDF5 file or TileDB store
if present. The function searches for known metadata dataset names
(`"greyordinates"`, `"fixels"`, `"voxels"`) at the top level of the
storage backend. It uses the backend retained when the `ModelArray` was
constructed, so an explicit constructor choice is preserved.

## Usage

``` r
elementMetadata(x)

# S4 method for class 'ModelArray'
elementMetadata(x)
```

## Arguments

- x:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

## Value

A matrix or data.frame of element metadata if found, or `NULL` if no
known metadata dataset exists.

## See also

[`nElements`](https://pennlinc.github.io/ModelArray/reference/nElements.md),
[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5", scalar_types = c("FD"))
meta <- elementMetadata(ma)
if (!is.null(meta)) head(meta)
} # }
```
