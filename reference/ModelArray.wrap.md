# Run a user-supplied function for element-wise data

\`ModelArray.wrap\` runs a user-supplied function \`FUN\` for each
requested element and returns a tibble/data.frame of results combined
across elements.

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
  ...
)
```

## Arguments

- FUN:

  A function that accepts at least \`data\` and returns a one-row
  data.frame/tibble, a named list, or a named vector of results for that
  element.

- data:

  ModelArray class

- phenotypes:

  A data.frame with cohort covariates; must contain column
  \`source_file\`. It must match \`sources(data)\[\[scalar\]\]\` order
  and contents (reordered if needed).

- scalar:

  Character, the scalar name to analyze

- element.subset:

  Optional integer vector of element indices (1-based); defaults to all

- num.subj.lthr.abs:

  Integer lower threshold for subjects with finite values (default 10)

- num.subj.lthr.rel:

  Relative lower threshold (0-1) (default 0.2)

- verbose:

  TRUE/FALSE to print messages

- pbar:

  TRUE/FALSE to show progress bar

- n_cores:

  Positive integer number of CPU cores

- on_error:

  Character: one of "stop", "skip", or "debug". When an error occurs in
  the user function for an element, choose whether to stop, skip
  returning all-NaN values for that element, or drop into \`browser()\`
  (if interactive) then skip. Default: "stop".

- ...:

  Additional arguments forwarded to \`FUN\`

## Value

Tibble/data.frame with one row per element and first column
\`element_id\`

## Details

This provides a generic framework reusing ModelArray's per-element
looping, alignment, subject-thresholding, and parallelization. The user
function will be called as \`FUN(data = dat, ...)\` where \`dat\` is
\`phenotypes\` with the response column named \`scalar\` appended for
the current element. The return value from \`FUN\` for a single element
must be either a one-row data.frame/tibble, a named list, or a named
atomic vector. The column names from the first successful element
determine the final schema.

Note: \`ModelArray.wrap\` never performs any p-value corrections or
modifications. If you need adjusted p-values (e.g., FDR), implement them
inside your provided \`FUN\`.
