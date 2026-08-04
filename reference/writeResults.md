# Write outputs from element-wise statistical analysis to storage

Creates a group named `analysis_name` under `/results/` in the storage
backend, then writes the statistical results data.frame (i.e. for one
analysis) into it as `results_matrix` along with column names.

## Usage

``` r
writeResults(
  fn.output,
  df.output,
  analysis_name = "myAnalysis",
  overwrite = TRUE,
  backend = c("auto", "hdf5", "tiledb")
)
```

## Arguments

- fn.output:

  Character. The HDF5 (`.h5`) filename or TileDB (`.tdb`) store for the
  output. Use an absolute path if you encounter file-not-found errors.

- df.output:

  A data.frame of element-wise statistical results, as returned by
  [`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
  [`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
  or
  [`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md).
  Must inherit from `data.frame`.

- analysis_name:

  Character. The name for this set of results. Used as the group name
  under `/results/` in the HDF5 file. Default is `"myAnalysis"`.

- overwrite:

  Logical. If a group with the same `analysis_name` already exists in
  the storage backend, whether to overwrite it (`TRUE`) or skip with a
  warning (`FALSE`). Default is `TRUE`.

- backend:

  Character. Storage backend: `"auto"` (default), `"hdf5"`, or
  `"tiledb"`. Auto-detection resolves to TileDB for a `.tdb` path, or
  for a directory that already contains a `scalars/` or `results/`
  subdirectory; everything else is treated as HDF5. To create a new
  TileDB store at a path without a `.tdb` suffix, pass `"tiledb"`
  explicitly.

## Value

Invisible `NULL`. Called for its side effect of writing results.

## Details

The results are stored at `/results/<analysis_name>/results_matrix` with
column names saved as a separate HDF5 dataset or TileDB metadata.

If any column of `df.output` is not numeric or integer, it is coerced to
numeric via [`factor()`](https://rdrr.io/r/base/factor.html) and the
factor levels are saved as a look-up table at
`/results/<analysis_name>/lut_forcol<i>`.

**Debugging tip:** If you encounter
`"Error in H5File.open(filename, mode, file_create_pl, file_access_pl)"`,
check if the message mentions "No such file or directory". Try using an
absolute path for the `fn.output` argument.

## See also

[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md),
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
which produce the `df.output`,
[`results`](https://pennlinc.github.io/ModelArray/reference/results.md)
for reading results back from a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md),
[`ModelArraySummary`](https://pennlinc.github.io/ModelArray/reference/ModelArraySummary.md)
for inspecting what has been written.

## Examples

``` r
if (FALSE) { # \dontrun{
ma <- ModelArray("data.h5", scalar_types = c("FD"))
phenotypes <- read.csv("cohort.csv")

results <- ModelArray.lm(
  FD ~ age + sex,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FD"
)

writeResults(
  fn.output = "data.h5",
  df.output = results,
  analysis_name = "lm_age_sex",
  overwrite = TRUE
)

# Verify
ModelArraySummary("data.h5")
} # }
```
