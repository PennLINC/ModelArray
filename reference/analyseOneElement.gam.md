# Fit a GAM for a single element Returns metadata (column names, smooth term names, parametric term names, and the smoothing parameter criterion attribute name) used by [`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md) to initialise the output data.frame. When `flag_initiate = FALSE`, it returns a numeric vector representing one row of the final results matrix.

If the number of subjects with finite scalar values does not exceed
`num.subj.lthr`, the element is skipped and all statistics are set to
`NaN`.

## Usage

``` r
analyseOneElement.gam(
  i_element,
  formula = NULL,
  modelarray = NULL,
  phenotypes = NULL,
  scalar = NULL,
  var.smoothTerms = c("statistic", "p.value"),
  var.parametricTerms = c("estimate", "statistic", "p.value"),
  var.model = c("dev.expl"),
  num.subj.lthr,
  num.stat.output = NULL,
  flag_initiate = FALSE,
  flag_sse = FALSE,
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
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

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

- var.smoothTerms:

  Character vector. Statistics to extract for smooth terms from
  [`tidy.gam`](https://broom.tidymodels.org/reference/tidy.gam.html)
  with `parametric = FALSE` (e.g. `"edf"`, `"ref.df"`, `"statistic"`,
  `"p.value"`).

- var.parametricTerms:

  Character vector. Statistics to extract for parametric terms from
  [`tidy.gam`](https://broom.tidymodels.org/reference/tidy.gam.html)
  with `parametric = TRUE` (e.g. `"estimate"`, `"std.error"`,
  `"statistic"`, `"p.value"`).

- var.model:

  Character vector. Statistics to extract for the overall model from
  [`glance.gam`](https://broom.tidymodels.org/reference/glance.gam.html)
  and [`summary.gam`](https://rdrr.io/pkg/mgcv/man/summary.gam.html)
  (e.g. `"adj.r.squared"`, `"dev.expl"`, `"sp.criterion"`).

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

- flag_sse:

  Logical. If `TRUE`, also compute the error sum of squares
  (`model.sse`) for the model, which is needed for partial R-squared
  calculations in
  [`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md).
  Default: `FALSE`.

- on_error:

  Character. One of `"stop"`, `"skip"`, or `"debug"`. When an error
  occurs fitting the model: `"stop"` halts execution; `"skip"` returns
  all-`NaN` for this element; `"debug"` drops into
  [`browser`](https://rdrr.io/r/base/browser.html) (if interactive) then
  skips. Default: `"stop"`.

- ctx:

  A precomputed context list from `.build_gam_context()`, or `NULL`.

- ...:

  Additional arguments passed to
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

## Value

If `flag_initiate = TRUE`, a list with components:

- column_names:

  Character vector of output column names.

- list.smoothTerms:

  Character vector of smooth term names.

- list.parametricTerms:

  Character vector of parametric term names.

- sp.criterion.attr.name:

  Character. The name attribute of the smoothing parameter selection
  criterion (e.g. `"REML"` or `"GCV.Cp"`).

If `flag_initiate = FALSE`, a numeric vector of length `num.stat.output`
with `element_id` (0-based) first and requested statistics in subsequent
positions. All-`NaN` (except `element_id`) if the element was skipped.

## See also

[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
which calls this function iteratively,
[`analyseOneElement.lm`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.lm.md)
for the linear model equivalent,
[`analyseOneElement.wrap`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.wrap.md)
for user-supplied functions.
