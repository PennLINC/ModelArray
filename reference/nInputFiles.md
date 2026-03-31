# Number of input files in a ModelArray

Returns the number of input files (i.e., subjects or source files) for a
given scalar in a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).
This is the column count of the scalar matrix.

## Usage

``` r
nInputFiles(x, scalar = NULL)

# S4 method for class 'ModelArray'
nInputFiles(x, scalar = NULL)
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

Integer. The number of input files (columns) in the scalar matrix.

## See also

[`nElements`](https://pennlinc.github.io/ModelArray/reference/nElements.md),
[`sources`](https://pennlinc.github.io/ModelArray/reference/sources.md),
[`scalarNames`](https://pennlinc.github.io/ModelArray/reference/scalarNames.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5", scalar_types = c("FD"))
nInputFiles(ma)
} # }
```
