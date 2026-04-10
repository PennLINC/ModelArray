# Run a user-supplied function for element-wise data

`ModelArray.wrap` runs a user-supplied function `FUN` at each requested
element and returns a tibble of results combined across elements.

## Usage

``` r
ModelArray.wrap(
  FUN,
  data,
  phenotypes,
  scalar,
  element.subset = NULL,
  num.subj.lthr.abs = 10,
  num.subj.lthr.rel = 0.2,
  verbose = TRUE,
  pbar = TRUE,
  n_cores = 1,
  on_error = "stop",
  write_scalar_name = NULL,
  write_scalar_file = NULL,
  write_scalar_columns = NULL,
  write_scalar_column_names = NULL,
  write_scalar_flush_every = 1000L,
  write_scalar_storage_mode = "double",
  write_scalar_compression_level = 4L,
  write_results_name = NULL,
  write_results_file = NULL,
  write_results_flush_every = 1000L,
  write_results_storage_mode = "double",
  write_results_compression_level = 4L,
  return_output = TRUE,
  ...
)
```

## Arguments

- FUN:

  A function that accepts at least `data` (a data.frame) and returns a
  one-row data.frame/tibble, a named list, or a named vector of results
  for that element.

- data:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

- phenotypes:

  A data.frame of the cohort with columns of independent variables and
  covariates to be added to the model. It must contain a column called
  `"source_file"` whose entries match those in
  `sources(data)[[scalar]]`.

- scalar:

  Character. The name of the element-wise scalar to analyse. Must be one
  of `names(scalars(data))`.

- element.subset:

  Integer vector of element indices (1-based) to run. Default is `NULL`,
  i.e. all elements in `data`.

- num.subj.lthr.abs:

  Integer. Lower threshold for the absolute number of subjects with
  finite scalar values (not `NaN`, `NA`, or `Inf`) required per element.
  Elements below this threshold are skipped (outputs set to `NaN`).
  Default is 10.

- num.subj.lthr.rel:

  Numeric between 0 and 1. Lower threshold for the proportion of
  subjects with finite values. Used together with `num.subj.lthr.abs`
  (the effective threshold is the maximum of the two). Default is 0.2.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

- pbar:

  Logical. Show progress bar. Default `TRUE`.

- n_cores:

  Positive integer. Number of CPU cores for parallel processing via
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html). Default is 1
  (serial).

- on_error:

  Character: one of `"stop"`, `"skip"`, or `"debug"`. When an error
  occurs fitting one element: `"stop"` halts execution; `"skip"` returns
  all-`NaN` for that element; `"debug"` drops into
  [`browser`](https://rdrr.io/r/base/browser.html) (if interactive) then
  skips. Default: `"stop"`.

- write_scalar_name:

  Optional character. If provided, selected output columns are written
  into `scalars/<write_scalar_name>/values` in the HDF5 file specified
  by `write_scalar_file`.

- write_scalar_file:

  Optional character. HDF5 output file path. Required when
  `write_scalar_name` is provided.

- write_scalar_columns:

  Optional character or integer vector selecting which output columns to
  save as scalar values. If `NULL` (default), uses all output columns
  except `element_id`.

- write_scalar_column_names:

  Optional character vector of source names saved as scalar
  `column_names`. If `NULL`, uses `phenotypes$source_file`.

- write_scalar_flush_every:

  Positive integer. Elements per write block for scalar writes. Default
  1000.

- write_scalar_storage_mode:

  Character. Storage mode for scalar writes (e.g. `"double"`). Default
  `"double"`.

- write_scalar_compression_level:

  Integer 0–9. Gzip compression level for scalar writes. Default 4.

- write_results_name:

  Optional character. If provided, results are incrementally written to
  `results/<write_results_name>/results_matrix` in the HDF5 file
  specified by `write_results_file`.

- write_results_file:

  Optional character. HDF5 file path for incremental result writes.
  Required when `write_results_name` is provided.

- write_results_flush_every:

  Positive integer. Number of elements per write block. Default 1000.

- write_results_storage_mode:

  Character. Storage mode for HDF5 writes (e.g. `"double"`). Default
  `"double"`.

- write_results_compression_level:

  Integer 0–9. Gzip compression level for HDF5 writes. Default 4.

- return_output:

  Logical. If `TRUE` (default), return the combined data.frame. If
  `FALSE`, return `invisible(NULL)`; useful when writing large outputs
  directly to HDF5.

- ...:

  Additional arguments forwarded to `FUN`.

## Value

If `flag_initiate = TRUE`, a list with one component:

- column_names:

  Character vector. The column names derived from the return value of
  `user_fun`, with `"element_id"` prepended. For unnamed list or atomic
  returns, columns are named `v1`, `v2`, etc. Set to `NaN` if the
  element was skipped or errored.

If `flag_initiate = FALSE`, a numeric vector of length `num.stat.output`
with `element_id` (0-based) first and the coerced output of `user_fun`
in subsequent positions. All-`NaN` (except `element_id`) if the element
was skipped or if an error occurred with `on_error = "skip"`.

## Details

This provides a generic framework reusing ModelArray's per-element
looping, alignment, subject-thresholding, and parallelization. The user
function is called as `FUN(data = dat, ...)` where `dat` is `phenotypes`
with all scalar columns appended for the current element. The return
value from `FUN` for a single element must be one of:

- a one-row `data.frame` or `tibble`

- a named list

- a named atomic vector

The column names from the first successful element determine the final
schema.

Note: `ModelArray.wrap` never performs any p-value corrections or
modifications. If you need adjusted p-values (e.g. FDR), implement them
inside `FUN`.

Use
[`exampleElementData`](https://pennlinc.github.io/ModelArray/reference/exampleElementData.md)
to construct a sample per-element data.frame for testing your function
before committing to a full run.

## See also

[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
for linear models,
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
for GAMs,
[`exampleElementData`](https://pennlinc.github.io/ModelArray/reference/exampleElementData.md)
for building test data,
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
for the input class.

## Examples

``` r
{
if (FALSE) { # \dontrun{
ma <- ModelArray("path/to/data.h5", scalar_types = c("FD"))
phenotypes <- read.csv("cohort.csv")

# Simple custom function
my_fun <- function(data, ...) {
  mod <- lm(FD ~ age + sex, data = data)
  tidy_out <- broom::tidy(mod)
  # Return a one-row tibble
  tibble::tibble(
    age_estimate = tidy_out$estimate[tidy_out$term == "age"],
    age_pvalue   = tidy_out$p.value[tidy_out$term == "age"]
  )
}


# Test on one element first
test_df <- exampleElementData(ma, scalar = "FD",
                               i_element = 1,
                               phenotypes = phenotypes)
my_fun(data = test_df)

# Run across all elements
results <- ModelArray.wrap(
  FUN = my_fun,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FD"
)
} # }
}
```
