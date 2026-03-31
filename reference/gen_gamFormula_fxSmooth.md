# Generate GAM formula with factor-smooth interaction

Generates a formula in the format
`y ~ orderedFactor + s(x) + s(x, by = orderedFactor)`, where `y` is
`response.var`, `x` is `smooth.var`, and `orderedFactor` is
`factor.var`. The formula generated can be further modified, e.g. by
adding covariates.

## Usage

``` r
gen_gamFormula_fxSmooth(
  response.var,
  factor.var,
  smooth.var,
  phenotypes,
  reference.group = NULL,
  prefix.ordered.factor = "o",
  fx = TRUE,
  k = NULL
)
```

## Arguments

- response.var:

  Character. The variable name for the response (dependent variable),
  typically a scalar name like `"FD"`.

- factor.var:

  Character. The variable name for the factor. It should be an ordered
  factor in `phenotypes`. If not, an ordered factor will be generated as
  a new column, which requires `reference.group`.

- smooth.var:

  Character. The variable name for the smooth term main effect (e.g.
  `"age"`).

- phenotypes:

  A data.frame of the cohort with columns of independent variables,
  including `factor.var` and `smooth.var`.

- reference.group:

  Character. The reference (baseline) group for the ordered factor of
  `factor.var`. Required when `factor.var` in `phenotypes` is not
  already an ordered factor.

- prefix.ordered.factor:

  Character. Prefix for the ordered factor column name. Required when
  `factor.var` in `phenotypes` is not already an ordered factor. Default
  is `"o"`.

- fx:

  Logical. Passed to [`s`](https://rdrr.io/pkg/mgcv/man/s.html). If
  `TRUE` (recommended), the smooth is treated as fixed degrees of
  freedom. Default is `TRUE`.

- k:

  Integer or `NULL`. Basis dimension passed to
  [`s`](https://rdrr.io/pkg/mgcv/man/s.html) for both the main smooth
  and interaction terms. If `NULL` (default), uses the default from
  [`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html).

## Value

A list with two components:

- formula:

  The generated [`formula`](https://rdrr.io/r/stats/formula.html)
  object.

- phenotypes:

  The (possibly updated) data.frame. If `factor.var` was not already an
  ordered factor, a new column named
  `paste0(prefix.ordered.factor, factor.var)` is added. Otherwise
  identical to the input.

## Details

This helper exists because setting up factor-smooth interactions in
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) requires an ordered
factor and a specific formula structure. If `factor.var` in `phenotypes`
is not already an ordered factor, this function creates one using
`reference.group` as the baseline level and adds it as a new column
(named with `prefix.ordered.factor` prepended to `factor.var`).

The returned `phenotypes` data.frame must be used in the subsequent
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
call so that the ordered factor column is available to the model.

## See also

[`gen_gamFormula_contIx`](https://pennlinc.github.io/ModelArray/reference/gen_gamFormula_contIx.md)
for continuous-by-continuous interactions,
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
which accepts the generated formula.

## Examples

``` r
if (FALSE) { # \dontrun{
phenotypes <- read.csv("cohort.csv")

# factor.var is not yet ordered - function creates it
result <- gen_gamFormula_fxSmooth(
  response.var = "FD",
  factor.var = "sex",
  smooth.var = "age",
  phenotypes = phenotypes,
  reference.group = "female"
)
result$formula

# Use the updated phenotypes (contains the ordered factor column)
results <- ModelArray.gam(
  result$formula,
  data = ma,
  phenotypes = result$phenotypes,
  scalar = "FD"
)
} # }
```
