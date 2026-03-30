# Load element-wise data from an HDF5 file as a ModelArray object

Reads scalar matrices and (optionally) saved analysis results from an
HDF5 file and returns a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
object.

## Usage

``` r
modelarray(filepath, scalar_types = c("FD"), analysis_names = character(0))
```

## Arguments

- filepath:

  Character. Path to an existing HDF5 (`.h5`) file containing
  element-wise scalar data.

- scalar_types:

  Character vector. Names of scalar groups to read from `/scalars/` in
  the HDF5 file. Default is `c("FD")`. Must match group names in the
  file; verify with `rhdf5::h5ls(filepath)`.

- analysis_names:

  Character vector. Subfolder names under `/results/` to load. Default
  is `character(0)` (no results loaded).

## Value

A
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
object.

## Details

The constructor reads each scalar listed in `scalar_types` from
`/scalars/<scalar_type>/values` in the HDF5 file, wrapping them as
DelayedArray objects for lazy access. Source filenames are extracted
from HDF5 attributes or companion datasets.

If `analysis_names` is non-empty, saved results are read from
`/results/<name>/results_matrix` along with column name metadata.

**Debugging tip:** If you encounter
`"Error in h(simpleError(msg, call)) : error in evaluating the argument 'seed'..."`,
check that `scalar_types` matches the groups present in the file.
Inspect with `rhdf5::h5ls(filepath)`.

## See also

[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
for the class definition,
[`h5summary`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)
for inspecting an HDF5 file without loading,
[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md),
[`sources`](https://pennlinc.github.io/ModelArray/reference/sources.md)
for accessing loaded data.

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("path/to/data.h5", scalar_types = c("FD", "FC"))
ma
scalars(ma)
} # }
```
