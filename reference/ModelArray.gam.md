# Fit element-wise generalized additive models

`ModelArray.gam` fits a generalized additive model at each requested
element in a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
and returns a tibble of requested model statistics. There is no
model-level p-value for GAMs, so there is no `correct.p.value.model`
argument.

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
  and [`summary.gam`](https://rdrr.io/pkg/mgcv/man/summary.gam.html).
  See Details.

- changed.rsq.term.index:

  A list of positive integers. Each value is the index of a term on the
  right-hand side of `formula` for which delta adjusted R-squared and
  partial R-squared should be computed. Usually the term of interest is
  a smooth term or interaction term. Default `NULL` (not computed). See
  Details.

- correct.p.value.smoothTerms:

  Character vector. P-value correction method(s) for each smooth term.
  Default: `"fdr"`. See Details.

- correct.p.value.parametricTerms:

  Character vector. P-value correction method(s) for each parametric
  term. Default: `"fdr"`. See Details.

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
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (e.g.
  `method = "REML"`).

## Value

A data.frame with one row per element. The first column is `element_id`
(0-based). Remaining columns contain the requested statistics, named as
`<term>.<statistic>` for per-term statistics and `model.<statistic>` for
model-level statistics. Smooth term names are normalized (e.g.
`s_age.statistic`). If p-value corrections were requested, additional
columns are appended with the correction method as suffix (e.g.
`s_age.p.value.fdr`). If `changed.rsq.term.index` was requested,
additional columns `<term>.delta.adj.rsq` and `<term>.partial.rsq` are
appended.

## Details

You may request returning specific statistical variables by setting
`var.*`, or you can get all by setting `full.outputs = TRUE`. Note that
statistics covered by `full.outputs` or `var.*` are the ones from
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`broom::glance()`](https://generics.r-lib.org/reference/glance.html),
and `summary.gam()` only, and do not include corrected p-values. However
FDR-corrected p-values (`"fdr"`) are generated by default.

List of acceptable statistic names for each of `var.*`:

- `var.smoothTerms`: `c("edf", "ref.df", "statistic", "p.value")`; From
  `broom::tidy(parametric = FALSE)`.

- `var.parametricTerms`:
  `c("estimate", "std.error", "statistic", "p.value")`; From
  `broom::tidy(parametric = TRUE)`.

- `var.model`:
  `c("adj.r.squared", "dev.expl", "sp.criterion", "scale", "df", "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs")`;
  From
  [`broom::glance()`](https://generics.r-lib.org/reference/glance.html)
  and [`summary.gam`](https://rdrr.io/pkg/mgcv/man/summary.gam.html).

Smooth term names in the output are normalized: `s(age)` becomes
`s_age`, `ti(x,z)` becomes `ti_x_z`, and `s(x):oFactor` becomes
`s_x_BYoFactor`.

For p-value corrections (arguments `correct.p.value.*`), supported
methods include all methods in `p.adjust.methods` except `"none"`. You
can request more than one method. FDR-corrected p-values (`"fdr"`) are
calculated by default. Turn it off by setting to `"none"`.

When `changed.rsq.term.index` is provided, a reduced model (dropping the
specified term) is fit at each element to compute delta adjusted
R-squared and partial R-squared. This approximately doubles execution
time per requested term. The term index refers to the position on the
right-hand side of `formula` (use `labels(terms(formula))` to see the
ordering).

Arguments `num.subj.lthr.abs` and `num.subj.lthr.rel` are mainly for
input data with subject-specific masks, i.e. currently only for volume
data. For fixel-wise data, you may ignore these arguments.

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
for the input class,
[`ModelArray`](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.html)
for the constructor,
[`exampleElementData`](https://pennlinc.github.io/ModelArray/reference/exampleElementData.md)
for testing formulas on a single element.

## Examples

``` r
{
if (FALSE) { # \dontrun{
# Fit GAM with default outputs
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

# Full outputs, no p-value correction
results_full <- ModelArray.gam(
  FD ~ s(age, fx = TRUE) + sex,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FD",
  full.outputs = TRUE,
  correct.p.value.smoothTerms = "none",
  correct.p.value.parametricTerms = "none"
)
} # }
}
```
