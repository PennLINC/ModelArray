# Run a user-supplied function for one element

\`analyseOneElement.wrap\` runs a user-supplied function \`FUN\` on a
single element's data. It prepares the per-element data by attaching the
element values as a new column named by \`scalar\` to the provided
\`phenotypes\` data.frame, then calls \`FUN(data = dat, ...)\`.

## Usage

``` r
analyseOneElement.wrap(
  i_element,
  user_fun,
  modelarray,
  phenotypes,
  scalar,
  num.subj.lthr,
  num.stat.output = NULL,
  flag_initiate = FALSE,
  on_error = "stop",
  ...
)
```

## Arguments

- i_element:

  An integer, the i_th element, starting from 1.

- user_fun:

  A function that accepts at least an argument named \`data\` (the
  per-element \`phenotypes\` with the response column appended) and
  returns a one-row data.frame/tibble, named list, or named vector.

- modelarray:

  ModelArray class

- phenotypes:

  A data.frame of the cohort with columns of independent variables and
  covariates

- scalar:

  A character. The name of the element-wise scalar to be analysed

- num.subj.lthr:

  The minimal number of subjects with valid value in input h5 file

- num.stat.output:

  The number of output stat metrics (including \`element_id\`). Required
  when \`flag_initiate = TRUE\`.

- flag_initiate:

  TRUE or FALSE, whether this is to initiate the new analysis to get
  column names

- on_error:

  Character: one of "stop", "skip", or "debug". When an error occurs
  while executing the user function, choose whether to stop, skip
  returning all-NaN values for this element, or drop into \`browser()\`
  (if interactive) then skip. Default: "stop".

- ...:

  Additional arguments forwarded to \`FUN\`

## Value

If \`flag_initiate==TRUE\`, returns a list with \`column_names\`. If
\`flag_initiate==FALSE\`, returns a numeric vector representing the
one-row result for this element with \`element_id\` as the first value.

## Details

The user-supplied \`FUN\` should return either a one-row
data.frame/tibble, a named list, or a named vector. The result will be
coerced to a one-row tibble and combined into the final results matrix
across elements.
