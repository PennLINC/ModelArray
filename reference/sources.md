# Source filenames of a ModelArray object

Retrieve the named list of source filename vectors from a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).
Each element of the list corresponds to one scalar and contains a
character vector of filenames, one per input file/subject.

## Usage

``` r
sources(x)

# S4 method for class 'ModelArray'
sources(x)
```

## Arguments

- x:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

## Value

A named list of character vectors. Names correspond to scalar names
(e.g. `"FD"`, `"FC"`).

## See also

[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md),
[`results`](https://pennlinc.github.io/ModelArray/reference/results.md),
[`scalarNames`](https://pennlinc.github.io/ModelArray/reference/scalarNames.md),
[`nInputFiles`](https://pennlinc.github.io/ModelArray/reference/nInputFiles.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5", scalar_types = c("FD"))
sources(ma) # named list
sources(ma)[["FD"]] # character vector of filenames
} # }
```
