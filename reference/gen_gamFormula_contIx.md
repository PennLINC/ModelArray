# Generate GAM formula with continuous-by-continuous interaction

Generates a formula in the format `y ~ ti(x) + ti(z) + ti(x, z)`, where
`y` is `response.var`, `x` is `cont1.var`, and `z` is `cont2.var`. The
formula generated can be further modified, e.g. by adding covariates.

## Usage

``` r
gen_gamFormula_contIx(response.var, cont1.var, cont2.var, fx = TRUE, k = NULL)
```

## Arguments

- response.var:

  Character. The variable name for the response (dependent variable),
  typically a scalar name like `"FD"`.

- cont1.var:

  Character. The name of the first continuous variable.

- cont2.var:

  Character. The name of the second continuous variable.

- fx:

  Logical. Passed to [`ti`](https://rdrr.io/pkg/mgcv/man/te.html). If
  `TRUE` (recommended), the smooth is treated as fixed degrees of
  freedom. Default is `TRUE`.

- k:

  Integer or `NULL`. Basis dimension passed to
  [`ti`](https://rdrr.io/pkg/mgcv/man/te.html) for all three terms (both
  main effects and the interaction). If `NULL` (default), uses the
  default from [`mgcv::ti()`](https://rdrr.io/pkg/mgcv/man/te.html).

## Value

A [`formula`](https://rdrr.io/r/stats/formula.html) object.

## Details

This helper uses [`ti`](https://rdrr.io/pkg/mgcv/man/te.html) (tensor
product interaction) terms so that the interaction `ti(x, z)` captures
only the interaction effect, separate from the main effects `ti(x)` and
`ti(z)`. This decomposition is important for interpretability and for
requesting `changed.rsq.term.index` in
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md).

## See also

[`gen_gamFormula_fxSmooth`](https://pennlinc.github.io/ModelArray/reference/gen_gamFormula_fxSmooth.md)
for factor-smooth interactions,
[`ModelArray.gam`](https://pennlinc.github.io/ModelArray/reference/ModelArray.gam.md)
which accepts the generated formula.

## Examples

``` r
if (FALSE) { # \dontrun{
formula <- gen_gamFormula_contIx(
  response.var = "FD",
  cont1.var = "age",
  cont2.var = "cognition"
)
formula

# Use in ModelArray.gam with changed R-squared for the interaction
results <- ModelArray.gam(
  formula,
  data = ma,
  phenotypes = phenotypes,
  scalar = "FD",
  changed.rsq.term.index = list(3)
)
} # }
```
