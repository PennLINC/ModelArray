# Fit element-wise linear models

`ModelArray.lm` fits a linear model at each requested element in a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
and returns a tibble of requested model statistics.

## Usage

``` r
ModelArray.lm(
  formula,
  data,
  phenotypes,
  scalar = NULL,
  element.subset = NULL,
  full.outputs = FALSE,
  var.terms = c("estimate", "statistic", "p.value"),
  var.model = c("adj.r.squared", "p.value"),
  correct.p.value.terms = c("fdr"),
  correct.p.value.model = c("fdr"),
  num.subj.lthr.abs = 10,
  num.subj.lthr.rel = 0.2,
  verbose = TRUE,
  pbar = TRUE,
  n_cores = 1,
  on_error = "stop",
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

- formula:

  Formula (passed to [`lm`](https://rdrr.io/r/stats/lm.html)).

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

- full.outputs:

  Logical. If `TRUE`, return the full set of statistics (ignoring
  `var.*` arguments). If `FALSE` (default), only return those requested
  in `var.*` and `correct.p.value.*`.

- var.terms:

  Character vector. Statistics to save per term, from
  [`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html). See
  Details.

- var.model:

  Character vector. Statistics to save for the overall model, from
  [`broom::glance()`](https://generics.r-lib.org/reference/glance.html).
  See Details.

- correct.p.value.terms:

  Character vector. P-value correction method(s) for each term. Default:
  `"fdr"`. See Details.

- correct.p.value.model:

  Character vector. P-value correction method(s) for the model-level
  p-value. Default: `"fdr"`. See Details.

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

  Additional arguments passed to
  [`lm`](https://rdrr.io/r/stats/lm.html).

## Value

A tibble with one row per element. The first column is `element_id`
(0-based). Remaining columns contain the requested statistics, named as
`<term>.<statistic>` for per-term statistics and `model.<statistic>` for
model-level statistics. If p-value corrections were requested,
additional columns are appended with the correction method as suffix
(e.g. `<term>.p.value.fdr`).

## Details

You may request returning specific statistical variables by setting
`var.*`, or you can get all by setting `full.outputs = TRUE`. Note that
statistics covered by `full.outputs` or `var.*` are the ones from
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`broom::glance()`](https://generics.r-lib.org/reference/glance.html)
only, and do not include corrected p-values. However FDR-corrected
p-values (`"fdr"`) are generated by default.

List of acceptable statistic names for each of `var.*`:

- `var.terms`: `c("estimate", "std.error", "statistic", "p.value")`; For
  interpretation please see
  [tidy.lm](https://broom.tidymodels.org/reference/tidy.lm.html).

- `var.model`:
  `c("r.squared", "adj.r.squared", "sigma", "statistic", "p.value", "df", "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs")`;
  For interpretation please see
  [glance.lm](https://broom.tidymodels.org/reference/glance.lm.html).

For p-value corrections (arguments `correct.p.value.*`), supported
methods include all methods in `p.adjust.methods` except `"none"`. You
can request more than one method. FDR-corrected p-values (`"fdr"`) are
calculated by default. Turn it off by setting to `"none"`.

Arguments `num.subj.lthr.abs` and `num.subj.lthr.rel` are mainly for
input data with subject-specific masks, i.e. currently only for volume
data. For fixel-wise data, you may ignore these arguments.

## See also

[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
for generalized additive models,
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
for user-supplied functions,
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
for the input class,
[`ModelArray`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html)
for the constructor,
[`exampleElementData`](https://pennlinc.github.io/ModelArray/reference/exampleElementData.md)
for testing formulas on a single element.

## Examples

``` r
if (FALSE) { # interactive()
ma <- ModelArray("path/to/data.h5", scalar_types = c("FD"))
phenotypes <- read.csv("cohort.csv")

# Fit linear model with default outputs
results <- ModelArray.lm(
  FD ~ age + sex,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FD"
)
head(results)

# Full outputs, no p-value correction
results_full <- ModelArray.lm(
  FD ~ age + sex,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FD",
  full.outputs = TRUE,
  correct.p.value.terms = "none",
  correct.p.value.model = "none"
)
}
```
