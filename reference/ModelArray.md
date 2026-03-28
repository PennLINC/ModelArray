# ModelArray class

A ModelArray wraps one or more element-wise scalar matrices (e.g., FD,
FC, log_FC for fixel data) read lazily via DelayedArray, along with any
previously saved analysis results. The object holds references to the
underlying HDF5 file and reads data on demand, making it suitable for
large-scale neuroimaging datasets.

Reads scalar matrices and (optionally) saved analysis results from an
HDF5 file and returns a ModelArray object.

## Usage

``` r
ModelArray(filepath, scalar_types = c("FD"), analysis_names = character(0))

ModelArray(filepath, scalar_types = c("FD"), analysis_names = character(0))
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

A ModelArray object.

## Details

ModelArray is an S4 class that represents element-wise scalar data and
associated statistical results backed by an HDF5 file on disk.

Each scalar in the HDF5 file is stored at `/scalars/<name>/values` as a
matrix of elements (rows) by source files (columns). Source filenames
are read from HDF5 attributes or companion datasets. Analysis results,
if present, live under `/results/<analysis_name>/results_matrix`.

ModelArray objects are typically created with the `ModelArray`
constructor function. Element-wise models are fit with
[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
or
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md).

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

## Slots

- `sources`:

  A named list of character vectors. Each element corresponds to a
  scalar and contains the source filenames (one per input file/subject).

- `scalars`:

  A named list of DelayedArray matrices. Each matrix has elements as
  rows and source files as columns.

- `results`:

  A named list of analysis results. Each element is itself a list
  containing at minimum `results_matrix` (a DelayedArray).

- `path`:

  Character. Path(s) to the HDF5 file(s) on disk.

## See also

`ModelArray` for the constructor,
[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
for analysis,
[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md),
[`sources`](https://pennlinc.github.io/ModelArray/reference/sources.md),
[`results`](https://pennlinc.github.io/ModelArray/reference/results.md)
for accessors.

ModelArray for the class definition,
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
