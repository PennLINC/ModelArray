# Example per-element data.frame for user functions

Constructs a per-element data.frame from a
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
that mirrors the `data` argument passed to user functions by
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md).
This is useful for testing and debugging user-supplied functions outside
of the full element-wise analysis loop.

## Usage

``` r
exampleElementData(x, ...)

# S4 method for class 'ModelArray'
exampleElementData(x, scalar = "FD", i_element = 1L, phenotypes)
```

## Arguments

- x:

  A
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md)
  object.

- ...:

  Additional arguments passed to the method (currently unused at the
  generic level).

- scalar:

  Character. The name of the element-wise scalar to append as a column.
  Must be one of `names(scalars(x))`. Default is `"FD"`.

- i_element:

  Integer. The 1-based index of the element whose values should be
  extracted. Must be between 1 and the number of elements for the given
  scalar.

- phenotypes:

  A data.frame of the cohort with columns of independent variables and
  covariates. Must have the same number of rows as the number of source
  files in the
  [ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).

## Value

A data.frame: the input `phenotypes` with one additional column named by
`scalar` containing that element's values across all subjects.

## Details

Returns a copy of `phenotypes` with an extra column named by `scalar`
populated with the selected element's values from the
[ModelArray](https://pennlinc.github.io/ModelArray/reference/ModelArray-class.md).
This mirrors the per-element data that
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
passes to user functions (as `data = dat`).

Use this to verify that your custom function works correctly on a single
element before committing to a full
[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md)
run.

## See also

[`ModelArray.wrap`](https://pennlinc.github.io/ModelArray/reference/ModelArray.wrap.md),
[`scalars`](https://pennlinc.github.io/ModelArray/reference/scalars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
ma <- ModelArray(h5_path, scalar_types = c("FD"))
phen <- read.csv(csv_path)
df <- exampleElementData(ma, scalar = "FD", i_element = 1, phenotypes = phen)

# Now test your custom function on this single element:
my_fun <- function(data, ...) {
  mod <- lm(FD ~ age + sex, data = data)
  broom::tidy(mod)
}
my_fun(data = df)
} # }
```
