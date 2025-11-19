# Example per-element data.frame for user functions

Generic for constructing a per-element data.frame from a \`ModelArray\`.
See the \`ModelArray\` method for details.

Returns a copy of \`phenotypes\` with an extra column named by
\`scalar\` populated with the selected element's values from the
\`ModelArray\`. This mirrors the per-element data that
\`ModelArray.wrap\` passes to user functions (\`data = dat\`).

## Usage

``` r
exampleElementData(x, ...)

# S4 method for class 'ModelArray'
exampleElementData(x, scalar = "FD", i_element = 1L, phenotypes)
```

## Arguments

- x:

  An ModelArray object

- ...:

  Additional arguments (ignored)

- scalar:

  A character. The name of the element-wise scalar to append

- i_element:

  An integer, the i_th element (1-based)

- phenotypes:

  A data.frame of the cohort with independent variables/covariates

## Value

A data.frame with the additional response column named by \`scalar\`

## Details

Example per-element data.frame for user functions

## Examples

``` r
if (FALSE) { # \dontrun{
h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
ma <- ModelArray(h5_path, scalar_types = c("FD"))
phen <- read.csv(csv_path)
df <- exampleElementData(ma, scalar = "FD", i_element = 1, phenotypes = phen)
} # }
```
