# Write outputs from element-wise statistical analysis to an HDF5 file

Creates a group named `analysis_name` under `/results/` in the HDF5
file, then writes the statistical results data.frame (i.e. for one
analysis) into it as `results_matrix` along with column names.

## Usage

``` r
writeResults(
  fn.output,
  df.output,
  analysis_name = "myAnalysis",
  overwrite = TRUE
)
```

## Arguments

- fn.output:

  Character. The HDF5 (`.h5`) filename for the output. The file must
  already exist; use an absolute path if you encounter file-not-found
  errors.

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
  the HDF5 file, whether to overwrite it (`TRUE`) or skip with a warning
  (`FALSE`). Default is `TRUE`.

## Value

Invisible `NULL`. Called for its side effect of writing results to the
HDF5 file.

## Details

The results are stored at `/results/<analysis_name>/results_matrix` with
column names saved as a separate dataset at
`/results/<analysis_name>/column_names`.

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
[`h5summary`](https://pennlinc.github.io/ModelArray/reference/h5summary.md)
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
h5summary("data.h5")
} # }
```
