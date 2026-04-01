# Generate GAM formula with continuous\*continuous interaction

This function will generate a formula in the following format:
`y ~ ti(x) + ti(z) + ti(x,z)`, where `y` is `response.var`, `x` is
`cont1.var`, and `z` is `cont2.var`. The formula generated could be
further modified, e.g. adding covariates.

## Usage

``` r
gen_gamFormula_contIx(response.var, cont1.var, cont2.var, fx = TRUE, k = NULL)
```

## Arguments

- response.var:

  character class, the variable name for response

- cont1.var:

  character class, the name of the first continuous variable

- cont2.var:

  character class, the name of the second continuous variable

- fx:

  TRUE or FALSE, to be used in smooth term s(). Recommend TRUE.

- k:

  integer, to be used in smooth term including the interaction term. If
  NULL (no entry), will use default value as in mgcv::s()

## Value

The formula generated
