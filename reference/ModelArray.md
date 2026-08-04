# Load element-wise data from an HDF5 file or TileDB store

Reads scalar matrices and (optionally) saved analysis results from an
HDF5 file or TileDB store and returns a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
object.

## Usage

``` r
ModelArray(
  filepath,
  scalar_types = c("FD"),
  analysis_names = character(0),
  backend = c("auto", "hdf5", "tiledb")
)
```

## Arguments

- filepath:

  Character. Path to an existing HDF5 (`.h5`) file or TileDB (`.tdb`)
  store containing element-wise scalar data.

- scalar_types:

  Character vector. Names of scalar groups to read from `/scalars/`.
  Default is `c("FD")`. Must match group names in the file.

- analysis_names:

  Character vector. Subfolder names under `/results/` to load. Default
  is `character(0)` (none).

- backend:

  Character. Storage backend: `"auto"` (default), `"hdf5"`, or
  `"tiledb"`. Auto-detection resolves to TileDB for a `.tdb` path, or
  for a directory that already contains a `scalars/` or `results/`
  subdirectory; everything else is treated as HDF5. Pass `"tiledb"`
  explicitly for a TileDB store that is neither.

## Value

A
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
object.

## Details

The constructor reads each scalar listed in `scalar_types` from
`/scalars/<scalar_type>/values`, wrapping them as
[DelayedArray::DelayedArray](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
objects. Source filenames are extracted from storage metadata or
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
[`ModelArraySummary`](https://pennlinc.github.io/ModelArray/reference/ModelArraySummary.md)
for inspecting storage.

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("path/to/data.h5", scalar_types = c("FD"))
ma
} # }
```
