# Run a user-supplied function for a single element

Runs a user-supplied function on one element's data, preparing the
per-element data.frame by attaching all scalar values as new columns to
the provided `phenotypes`. This is the per-element workhorse called
iteratively by
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md).

## Usage

``` r
analyseOneElement.wrap(
  i_element,
  user_fun,
  modelarray = NULL,
  phenotypes = NULL,
  scalar = NULL,
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

- user_fun:

  A function that accepts at least an argument named `data` (a
  data.frame: `phenotypes` with scalar columns appended for the current
  element) and returns one of:

  - A one-row `data.frame` or `tibble`. Multi-row returns will error.

  - A named list. Unnamed lists are accepted and auto-named as `v1`,
    `v2`, etc.

  - A named atomic vector. Unnamed vectors are accepted and auto-named
    as `v1`, `v2`, etc.

  All return values are coerced to a numeric vector internally.
  Unsupported types (e.g., environments, S4 objects) will error.

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

  A precomputed context list from `.build_wrap_context()`, or `NULL`.

- ...:

  Additional arguments forwarded to `user_fun`.

## Value

If `flag_initiate = TRUE`, a list with one component:

- column_names:

  Character vector. The column names derived from the return value of
  `user_fun`, with `"element_id"` prepended.

If `flag_initiate = FALSE`, a numeric vector of length `num.stat.output`
with `element_id` (0-based) first and the coerced output of `user_fun`
in subsequent positions. All-`NaN` (except `element_id`) if the element
was skipped or if an error occurred with `on_error = "skip"`.

## Details

Most users should call
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
directly, which handles looping, parallelisation, and result assembly.
`analyseOneElement.wrap` is exported for advanced use cases such as
debugging a single element or building custom analysis loops.

The user-supplied `user_fun` is called as `user_fun(data = dat, ...)`
where `dat` is `phenotypes` with columns appended for **all** scalars in
the
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
(not just the one named by `scalar`). The function should return a
one-row data.frame/tibble, a named list, or a named atomic vector. The
result is coerced to a numeric vector for assembly into the final
results matrix.

Unlike
[`analyseOneElement.lm`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.lm.md)
and
[`analyseOneElement.gam`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.gam.md),
this function appends **all** scalar columns (not just the response) and
checks for column name collisions between scalar names and existing
`phenotypes` columns.

If the number of subjects with finite values across all scalars does not
exceed `num.subj.lthr`, the element is skipped and all statistics are
set to `NaN`.

## See also

[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
which calls this function iteratively,
[`exampleElementData`](https://pennlinc.github.io/ModelArray/reference/exampleElementData.md)
for building a test data.frame matching the format passed to `user_fun`,
[`analyseOneElement.lm`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.lm.md)
for the linear model equivalent,
[`analyseOneElement.gam`](https://pennlinc.github.io/ModelArray/reference/analyseOneElement.gam.md)
for the GAM equivalent.
