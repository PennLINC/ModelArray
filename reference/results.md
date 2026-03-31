# Statistical results of a ModelArray object

Retrieve previously saved analysis results from a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).
When called with no additional arguments, returns the full named list of
all result sets. When called with a single analysis name, returns that
one result set.

## Usage

``` r
results(x, ...)

# S4 method for class 'ModelArray'
results(x, ...)
```

## Arguments

- x:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

- ...:

  Optional: a single character string giving the analysis name to
  extract. If omitted, the entire named list is returned.

## Value

If called with no extra arguments, a named list of result lists. If
called with an analysis name, the corresponding result list (containing
at minimum `results_matrix`).

## Details

Each result set is itself a list containing at minimum `results_matrix`
(a DelayedArray with elements as rows and statistics as columns). Column
names are stored alongside the matrix in the HDF5 file.

Results are only available if `analysis_names` was supplied when the
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
was constructed, or if results have been written back with
[`writeResults`](https://pennlinc.github.io/ModelArray/reference/writeResults.md).

## See also

[`sources`](https://pennlinc.github.io/ModelArray/reference/sources.md),
[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md),
[`analysisNames`](https://pennlinc.github.io/ModelArray/reference/analysisNames.md),
[`writeResults`](https://pennlinc.github.io/ModelArray/reference/writeResults.md),
[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5",
  scalar_types = c("FD"),
  analysis_names = c("lm_age")
)
results(ma) # named list of all results
results(ma, "lm_age") # single result set
results(ma, "lm_age")$results_matrix
} # }
```
