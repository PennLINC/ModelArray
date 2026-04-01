# Element-wise scalar data of a ModelArray object

Retrieve scalar matrices from a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).
When called with no additional arguments, returns the full named list of
all scalar matrices. When called with a single scalar name, returns that
one matrix.

## Usage

``` r
scalars(x, ...)

# S4 method for class 'ModelArray'
scalars(x, ...)
```

## Arguments

- x:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

- ...:

  Optional: a single character string giving the scalar name to extract.
  If omitted, the entire named list is returned.

## Value

If called with no extra arguments, a named list of DelayedArray matrices
(elements as rows, source files as columns). If called with a scalar
name, the corresponding single DelayedArray matrix.

## See also

[`sources`](https://pennlinc.github.io/ModelArray/reference/sources.md),
[`results`](https://pennlinc.github.io/ModelArray/reference/results.md),
[`scalarNames`](https://pennlinc.github.io/ModelArray/reference/scalarNames.md),
[`nElements`](https://pennlinc.github.io/ModelArray/reference/nElements.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5", scalar_types = c("FD", "FC"))
scalars(ma) # named list of all scalars
scalars(ma, "FD") # single DelayedArray matrix
} # }
```
