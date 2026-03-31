# Names of scalars in a ModelArray

Returns the names of all scalar datasets loaded in a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
(e.g. `"FD"`, `"FC"`, `"log_FC"`).

## Usage

``` r
scalarNames(x)

# S4 method for class 'ModelArray'
scalarNames(x)
```

## Arguments

- x:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

## Value

Character vector of scalar names.

## See also

[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md),
[`analysisNames`](https://pennlinc.github.io/ModelArray/reference/analysisNames.md),
[`nElements`](https://pennlinc.github.io/ModelArray/reference/nElements.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5", scalar_types = c("FD", "FC"))
scalarNames(ma)   # c("FD", "FC")
} # }
```
