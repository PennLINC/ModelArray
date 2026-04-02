# Load element-wise data from an HDF5 file

Reads scalar matrices and (optionally) saved analysis results from an
HDF5 file and returns a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
object.

## Usage

``` r
ModelArray(filepath, scalar_types = c("FD"), analysis_names = character(0))
```

## Arguments

- filepath:

  Character. Path to an existing HDF5 (`.h5`) file containing
  element-wise scalar data.

- scalar_types:

  Character vector. Names of scalar groups to read from `/scalars/` in
  the HDF5 file. Default is `c("FD")`. Must match group names in the
  file.

- analysis_names:

  Character vector. Subfolder names under `/results/` to load. Default
  is `character(0)` (none).

## Value

A
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
object.

## Details

The constructor reads each scalar listed in `scalar_types` from
`/scalars/<scalar_type>/values`, wrapping them as
[DelayedArray::DelayedArray](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
objects. Source filenames are extracted from HDF5 attributes or
companion datasets.

If `analysis_names` is non-empty, saved results are loaded from
`/results/<name>/results_matrix`.

**Debugging tip:** If you encounter
`"error in evaluating the argument 'seed'..."`, check that
`scalar_types` matches groups in the file. Inspect with
`rhdf5::h5ls(filepath)`.

## See also

[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
for the class definition,
[`h5summary`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)
for inspecting an HDF5 file.

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("path/to/data.h5", scalar_types = c("FD"))
ma
} # }
```
