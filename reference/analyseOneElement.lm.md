# Fit linear model for one element.

\`analyseOneElement.lm\` fits a linear model for one element data, and
returns requested model statistics.

## Usage

``` r
analyseOneElement.lm(
  i_element,
  formula,
  modelarray,
  phenotypes,
  scalar,
  var.terms,
  var.model,
  num.subj.lthr,
  num.stat.output = NULL,
  flag_initiate = FALSE,
  on_error = "stop",
  ...
)
```

## Arguments

- i_element:

  An integer, the i_th element, starting from 1. For initiating
  (flag_initiate = TRUE), use i_element=1

- formula:

  Formula (passed to \`stats::lm()\`)

- modelarray:

  ModelArray class

- phenotypes:

  A data.frame of the cohort with columns of independent variables and
  covariates to be added to the model.

- scalar:

  A character. The name of the element-wise scalar to be analysed

- var.terms:

  A list of characters. The list of variables to save for terms (got
  from \`broom::tidy()\`).

- var.model:

  A list of characters. The list of variables to save for the model (got
  from \`broom::glance()\`).

- num.subj.lthr:

  The minimal number of subjects with valid value in input h5 file, i.e.
  number of subjects with finite values (defined by \`is.finite()\`,
  i.e. not NaN or NA or Inf) in h5 file \> `num.subj.lthr`, then this
  element will be run normally; otherwise, this element will be skipped
  and statistical outputs will be set as NaN.

- num.stat.output:

  The number of output stat metrics (for generating all NaN stat when \#
  subjects does not meet criteria). This includes column \`element_id\`.
  This is required when flag_initiate = TRUE.

- flag_initiate:

  TRUE or FALSE, Whether this is to initiate the new analysis. If TRUE,
  it will return column names etc to be used for initiating data.frame;
  if FALSE, it will return the list of requested statistic values.

- on_error:

  Character: one of "stop", "skip", or "debug". When an error occurs
  while fitting one element, choose whether to stop, skip returning
  all-NaN values for that element, or drop into \`browser()\` (if
  interactive) then skip. Default: "stop".

- ...:

  Additional arguments for \`stats::lm()\`

## Value

If flag_initiate==TRUE, returns column names, and list of term names of
final results; if flag_initiate==FALSE, it will return the list of
requested statistic values for a element.

## Details

\`ModelArray.lm\` iteratively calls this function to get statistics for
all requested elements.
