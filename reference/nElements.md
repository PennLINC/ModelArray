# Number of elements in a ModelArray

Returns the number of elements (e.g., fixels or voxels) for a given
scalar in a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).
This is the row count of the scalar matrix.

## Usage

``` r
nElements(x, scalar = NULL)

# S4 method for class 'ModelArray'
nElements(x, scalar = NULL)
```

## Arguments

- x:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

- scalar:

  Optional character string. Name of the scalar to query. Defaults to
  the first scalar in `names(scalars(x))`.

## Value

Integer. The number of elements (rows) in the scalar matrix.

## See also

[`nInputFiles`](https://pennlinc.github.io/ModelArray/reference/nInputFiles.md),
[`numElementsTotal`](https://pennlinc.github.io/ModelArray/reference/numElementsTotal.md),
[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md),
[`scalarNames`](https://pennlinc.github.io/ModelArray/reference/scalarNames.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5", scalar_types = c("FD"))
nElements(ma)
nElements(ma, "FD")
} # }
```
