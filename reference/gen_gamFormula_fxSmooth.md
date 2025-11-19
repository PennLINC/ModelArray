# Generate GAM formula with factor-smooth interaction

This function will generate a formula in the following format:
`y ~ orderedFactor + s(x) + s(x, by=orderedFactor)`, where `y` is
`response.var`, `x` is `smooth.var`, and `orderedFactor` is
`factor.var` - see `factor.var` for more. The formula generated could be
further modified, e.g. adding covariates.

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

  character class, the variable name for response

- factor.var:

  character class, the variable name for factor. It should be an ordered
  factor. If not, it will generate it as a new column in \`phenotypes\`,
  which requires \`reference.group\`.

- smooth.var:

  character class, the variable name in smooth term as main effect

- phenotypes:

  data.frame class, the cohort matrix with columns of independent
  variables (including `factor.var` and `smooth.var`) to be added to the
  model

- reference.group:

  character class, the reference group for ordered factor of
  \`factor.var\`; required when \`factor.var\` in \`phenotypes\` is not
  an ordered factor.

- prefix.ordered.factor:

  character class, the prefix for ordered factor; required when
  \`factor.var\` in \`phenotypes\` is not an ordered factor.

- fx:

  TRUE or FALSE, to be used in smooth term s(). Recommend TRUE.

- k:

  integer, to be used in smooth term including the interaction term. If
  NULL (no entry), will use default value as in mgcv::s()

## Value

a list, including: 1) formula generated; 2) data.frame phenotypes -
updated if argument factor.var is not an ordered factor
