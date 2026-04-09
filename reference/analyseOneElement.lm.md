# Fit a linear model for a single element

Fits [`lm`](https://rdrr.io/r/stats/lm.html) on one element's data. When
a precomputed context (`ctx`) is provided, all loop-invariant work
(formula parsing, collision checks, source alignment) is skipped. When
`ctx` is `NULL`, the function falls back to computing everything inline
(legacy behaviour for direct calls / debugging).

If the number of subjects with finite scalar values (not `NaN`, `NA`, or
`Inf`) does not exceed `num.subj.lthr`, the element is skipped and all
statistics are set to `NaN`.

## Usage

``` r
analyseOneElement.lm(
  i_element,
  formula = NULL,
  modelarray = NULL,
  phenotypes = NULL,
  scalar = NULL,
  var.terms = c("estimate", "statistic", "p.value"),
  var.model = c("adj.r.squared", "p.value"),
  num.subj.lthr,
  num.stat.output = NULL,
  flag_initiate = FALSE,
  on_error = "stop",
  ctx = NULL,
  ...
)
```

## Arguments

- i_element:

  Integer. The 1-based index of the element to analyse.

- formula:

  A [`formula`](https://rdrr.io/r/stats/formula.html) passed to
  [`lm`](https://rdrr.io/r/stats/lm.html). Ignored when `ctx` is
  provided (the formula is taken from the context).

- modelarray:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object. Ignored when `ctx` is provided.

- phenotypes:

  A data.frame of the cohort with columns of independent variables and
  covariates. Must contain a `"source_file"` column matching
  `sources(modelarray)[[scalar]]`. Ignored when `ctx` is provided.

- scalar:

  Character. The name of the element-wise scalar to analyse. Must be one
  of `names(scalars(modelarray))`. Ignored when `ctx` is provided.

- var.terms:

  Character vector. Statistics to extract per term from
  [`tidy.lm`](https://broom.tidymodels.org/reference/tidy.lm.html) (e.g.
  `"estimate"`, `"statistic"`, `"p.value"`).

- var.model:

  Character vector. Statistics to extract for the overall model from
  [`glance.lm`](https://broom.tidymodels.org/reference/glance.lm.html)
  (e.g. `"adj.r.squared"`, `"p.value"`).

- num.subj.lthr:

  Numeric. The pre-computed minimum number of subjects with finite
  values required for this element to be analysed. Elements below this
  threshold are skipped. This value is typically computed by the parent
  function from `num.subj.lthr.abs` and `num.subj.lthr.rel`.

- num.stat.output:

  Integer or `NULL`. The total number of output columns (including
  `element_id`). Used when `flag_initiate = FALSE` to generate an
  all-`NaN` row for skipped elements. Must be `NULL` when
  `flag_initiate = TRUE`.

- flag_initiate:

  Logical. If `TRUE`, fit the model once and return metadata for
  initialising the output data.frame (column names and term names). If
  `FALSE`, return a numeric vector of results for this element.

- on_error:

  Character. One of `"stop"`, `"skip"`, or `"debug"`. When an error
  occurs fitting the model: `"stop"` halts execution; `"skip"` returns
  all-`NaN` for this element; `"debug"` drops into
  [`browser`](https://rdrr.io/r/base/browser.html) (if interactive) then
  skips. Default: `"stop"`.

- ctx:

  A precomputed context list from `.build_lm_context()`, or `NULL` for
  legacy inline computation.

- ...:

  Additional arguments passed to
  [`lm`](https://rdrr.io/r/stats/lm.html).

## Value

If `flag_initiate = TRUE`, a list with components:

- column_names:

  Character vector. The column names for the output data.frame, with
  `"element_id"` first.

- list.terms:

  Character vector. The names of the model terms (from
  [`tidy.lm`](https://broom.tidymodels.org/reference/tidy.lm.html)).

If `flag_initiate = FALSE`, a numeric vector of length `num.stat.output`
with `element_id` (0-based) as the first value and the requested
statistics in subsequent positions. All-`NaN` (except `element_id`) if
the element had insufficient valid subjects or if an error occurred with
`on_error = "skip"`.

## See also

[`ModelArray.lm`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md),
`.build_lm_context`,
[`analyseOneElement.gam`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.gam.md)
for the GAM equivalent,
[`analyseOneElement.wrap`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.wrap.md)
for user-supplied functions
