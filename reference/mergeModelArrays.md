# Merge multiple ModelArrays from different h5 files

Combines scalars from multiple ModelArray objects into a single
ModelArray, aligning subjects via phenotype columns. Uses
\`DelayedArray::acbind()\` for virtual column-binding — no h5 rewriting
is needed.

## Usage

``` r
mergeModelArrays(modelarrays, phenotypes_list, merge_on)
```

## Arguments

- modelarrays:

  A list of ModelArray objects

- phenotypes_list:

  A list of data.frames, one per ModelArray. Each must contain a
  \`source_file\` column matching its corresponding ModelArray's
  sources.

- merge_on:

  Character vector of column names to join phenotypes on (e.g.,
  \`c("subject_id")\`). These columns must exist in all phenotypes
  data.frames.

## Value

A list with:

- data:

  A combined ModelArray with scalars from all inputs

- phenotypes:

  The inner-joined phenotypes data.frame. The original \`source_file\`
  columns are renamed to \`source_file.\<scalar_name\>\`.
