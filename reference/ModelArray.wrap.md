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

- write_scalar_name:

  Optional character scalar name. If provided, \`ModelArray.wrap\`
  writes selected output columns into
  \`scalars/\<write_scalar_name\>/values\` while processing.

- write_scalar_file:

  Optional character HDF5 output filename used when
  \`write_scalar_name\` is provided.

- write_scalar_columns:

  Optional character/integer selector for output columns to save as
  scalar values. If NULL, uses all wrap output columns except
  \`element_id\`.

- write_scalar_column_names:

  Optional character vector of source names saved as scalar
  \`column_names\`. If NULL, uses \`phenotypes\$source_file\`.

- write_scalar_flush_every:

  Positive integer number of elements per write block.

- write_scalar_storage_mode:

  Storage mode for scalar writes (e.g., \`"double"\`).

- write_scalar_compression_level:

  Gzip compression level (0-9) for scalar writes.

- write_results_name:

  Optional analysis name for incremental writes to
  \`results/\<write_results_name\>/results_matrix\`.

- write_results_file:

  Optional HDF5 file path used when \`write_results_name\` is provided.

- write_results_flush_every:

  Positive integer number of elements per write block for results
  writes.

- write_results_storage_mode:

  Storage mode for results writes (e.g., \`"double"\`).

- write_results_compression_level:

  Gzip compression level (0-9) for results writes.

- return_output:

  If TRUE (default), return the combined data.frame. If FALSE, returns
  \`invisible(NULL)\`; useful when writing large outputs directly to
  HDF5.

- ...:

  Additional arguments forwarded to \`FUN\`

## Value

Tibble/data.frame with one row per element and first column
\`element_id\` when \`return_output = TRUE\`; otherwise
\`invisible(NULL)\`.

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
