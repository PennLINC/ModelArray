
<!-- TODO README.md is generated from README.Rmd. Please edit that file -->

# ModelArray

<!-- badges: start -->

[![CircleCI build
status](https://circleci.com/gh/PennLINC/ModelArray.svg?style=svg)](https://circleci.com/gh/PennLINC/ModelArray)
<!-- badges: end -->

ModelArray is a generalizable, memory-efficient R package for
statistical analysis of fixel data. Its features include:

-   Easy to use: set up your analysis with just several lines of codes;
-   Low memory requirement, even for large datasets;
-   At present, ModelArray supports linear models as well as generalized
    additive models (GAM) with penalized splines, which are particularly
    useful for studying nonlinear effects in lifespan data.

## Installation

Before you install ModelArray R package, if you are using Linux system,
please check if libhdf5-dev has been installed in your system:

``` console
foo@bar:~$ ldconfig -p | grep libhdf5*
```

If you got more than one line of outputs, congrats, you have libhdf5-dev
installed. Otherwise, please install it. For Ubuntu user, you may
install via:

``` console
foo@bar:~$ sudo apt-get update -y
foo@bar:~$ sudo apt-get install -y libhdf5-dev
```

You can install the latest version from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("PennLINC/ModelArray")
```

## How to use

Load the package into R via:

``` r
library(ModelArray)
```

We provide a quick start guide
[here](https://pennlinc.github.io/ModelArray/articles/ModelArray_basics.html).  
For documentations of functions, you can either refer to [this
webpage](https://pennlinc.github.io/ModelArray/reference/index.html); Or
in R Studio console: `?<function_name>`. For example:
`?ModelArray.lm()`  
For users who are not familiar with R, you may check out [this
webpage](https://pennlinc.github.io/ModelArray/articles/basic_r_intro.html).
