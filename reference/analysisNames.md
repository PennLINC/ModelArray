# Names of analyses in a ModelArray

Returns the names of all analysis result sets currently loaded in a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).
These correspond to subfolder names under `/results/` in the HDF5 file.

## Usage

``` r
analysisNames(x)

# S4 method for class 'ModelArray'
analysisNames(x)
```

## Arguments

- x:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

## Value

Character vector of analysis names. Returns `character(0)` if no
analyses have been loaded or saved.

## See also

[`results`](https://pennlinc.github.io/ModelArray/reference/results.md),
[`scalarNames`](https://pennlinc.github.io/ModelArray/reference/scalarNames.md),
[`writeResults`](https://pennlinc.github.io/ModelArray/reference/writeResults.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5",
  scalar_types = c("FD"),
  analysis_names = c("lm_age")
)
analysisNames(ma) # "lm_age"
} # }
```
