# Fit element-wise generalized additive models no model-level p-value for GAMs, so there is no `correct.p.value.model` argument.

Fit element-wise generalized additive models no model-level p-value for
GAMs, so there is no `correct.p.value.model` argument.

## Usage

``` r
ModelArray.gam(
  formula,
  data,
  phenotypes,
  scalar,
  element.subset = NULL,
  full.outputs = FALSE,
  var.smoothTerms = c("statistic", "p.value"),
  var.parametricTerms = c("estimate", "statistic", "p.value"),
  var.model = c("dev.expl"),
  changed.rsq.term.index = NULL,
  correct.p.value.smoothTerms = c("fdr"),
  correct.p.value.parametricTerms = c("fdr"),
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

  Formula (passed to [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)).

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

- var.smoothTerms:

  Character vector. Statistics to save for smooth terms, from
  `broom::tidy(parametric = FALSE)`. See Details.

- var.parametricTerms:

  Character vector. Statistics to save for parametric terms, from
  `broom::tidy(parametric = TRUE)`. See Details.

- var.model:

  Character vector. Statistics to save for the overall model, from
  [`broom::glance()`](https://generics.r-lib.org/reference/glance.html)
  and [`summary()`](https://rdrr.io/r/base/summary.html). See Details.

- changed.rsq.term.index:

  A list of positive integers. Each value is the index of a term on the
  right-hand side of `formula` for which delta adjusted R-squared and
  partial R-squared should be computed. Usually the term of interest is
  a smooth term or interaction term. Default `NULL` (not computed). See
  Details for warnings.

- correct.p.value.smoothTerms:

  Character vector. P-value correction method(s) for each smooth term.
  Default: `"fdr"`.

- correct.p.value.parametricTerms:

  Character vector. P-value correction method(s) for each parametric
  term. Default: `"fdr"`.

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
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

## Value

A tibble with one row per element. The first column is `element_id`
(0-based). Remaining columns contain the requested statistics, named as
`<term>.<statistic>`. If `changed.rsq.term.index` was requested,
additional columns `<term>.delta.adj.rsq` and `<term>.partial.rsq` are
appended.

## See also

[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
for linear models,
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
for user-supplied functions,
[`gen_gamFormula_fxSmooth`](https://pennlinc.github.io/ModelArray/reference/gen_gamFormula_fxSmooth.md)
and
[`gen_gamFormula_contIx`](https://pennlinc.github.io/ModelArray/reference/gen_gamFormula_contIx.md)
for formula helpers,
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
for the input class.

## Examples

``` r
{
if (FALSE) { # \dontrun{
ma <- ModelArray("path/to/data.h5", scalar_types = c("FD"))
phenotypes <- read.csv("cohort.csv")

results <- ModelArray.gam(
  FD ~ s(age, fx = TRUE) + sex,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FD"
)
head(results)

# With changed R-squared for the smooth term (term index 1)
results_rsq <- ModelArray.gam(
  FD ~ s(age, fx = TRUE) + sex,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FD",
  changed.rsq.term.index = list(1)
)
} # }
}
```
