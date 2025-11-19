# An intro to R

***UNDER DEVELOPMENT…***

This document is a quick introduction to basic R. Audience would be
users who are not familiar with R.

## R syntax

- Dot “.” is valid in a variable name or a function (e.g. variable
  `element.subset`, function
  [`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md))
- Syntax for formula: Below is an example formula:

``` r
formula <- FD ~ age + sex
```

As you can see even though we did not define variable “FD”,“age” or
“sex”, it’s a valid formula. Left hand side of “~” is dependent variable
(y, or response), and right hand side includes the independent variables
/ covariates.

## Prepare the phenotypes data.frame: manipulations

You may want to manipulate the data.frame `phenotypes` before passing it
to
[`ModelArray.lm()`](https://pennlinc.github.io/ModelArray/reference/ModelArray.lm.md)
etc functions for model fitting. Examples like de-mean and rescale
covariates. R package `dplyr` is useful and easy to use for
manipulation.

``` r
# install is first via: install.packages("dplyr")
library(dplyr) # load into memory
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union
