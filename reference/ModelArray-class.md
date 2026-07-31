# ModelArray class

A ModelArray wraps one or more element-wise scalar matrices (e.g., FD,
FC, log_FC for fixel data) read lazily via DelayedArray, along with any
previously saved analysis results. The object holds references to the
underlying storage and reads data on demand, making it suitable for
large-scale neuroimaging datasets.

Prints a summary of the ModelArray including file path, source count,
and any saved analysis names.

## Usage

``` r
ModelArray(
  filepath,
  scalar_types = c("FD"),
  analysis_names = character(0),
  backend = c("auto", "hdf5", "tiledb")
)

# S4 method for class 'ModelArray'
show(object)
```

## Arguments

- object:

  A ModelArray object.

## Value

Invisible `NULL`. Called for its side effect of printing to the console.

## Details

ModelArray is an S4 class that represents element-wise scalar data and
associated statistical results backed by an HDF5 file or TileDB store on
disk.

Each scalar is stored at `/scalars/<name>/values` as a matrix of
elements (rows) by source files (columns). Source filenames are read
from storage metadata or companion datasets. Analysis results, if
present, live under `/results/<analysis_name>/results_matrix`.

ModelArray objects are typically created with the
[`ModelArray`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html)
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

  A named list of
  [DelayedArray::DelayedArray](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
  matrices. Each matrix has elements as rows and source files as
  columns.

- `results`:

  A named list of analysis results. Each element is itself a list
  containing at minimum `results_matrix` (a
  [DelayedArray::DelayedArray](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)).

- `path`:

  Character. Path(s) to the HDF5 file(s) or TileDB store(s) on disk.

- `backend`:

  Character. Resolved storage backend(s) for `path`: `"hdf5"` or
  `"tiledb"`.

## See also

[`ModelArray`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html)
for the constructor,
[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
for analysis,
[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md),
[`sources`](https://pennlinc.github.io/ModelArray/reference/sources.md),
[`results`](https://pennlinc.github.io/ModelArray/reference/results.md)
for accessors.
