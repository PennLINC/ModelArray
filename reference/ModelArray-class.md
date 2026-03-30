# ModelArray class

A ModelArray wraps one or more element-wise scalar matrices (e.g., FD,
FC, log_FC for fixel data) read lazily via DelayedArray, along with any
previously saved analysis results. The object holds references to the
underlying HDF5 file and reads data on demand, making it suitable for
large-scale neuroimaging datasets.

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
