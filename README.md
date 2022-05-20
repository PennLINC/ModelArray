
<!-- TODO README.md is generated from README.Rmd. Please edit that file -->

# ModelArray

<!-- badges: start -->

[![CircleCI build
status](https://circleci.com/gh/PennLINC/ModelArray.svg?style=svg)](https://circleci.com/gh/PennLINC/ModelArray)
<!-- badges: end -->

`ModelArray` is a generalizable, memory-efficient R package for
statistical analysis of fixel data. Its features include:

-   Easy to use: set up your analysis with just several lines of codes;
-   Low memory requirement, even for large datasets;
-   At present, `ModelArray` supports linear models as well as
    generalized additive models (GAM) with penalized splines, which are
    particularly useful for studying nonlinear effects in lifespan data.

## Installation

Please refer to webpage
[Installation](https://pennlinc.github.io/ModelArray/articles/a01_installations.html)
for a full guidance of installation of `ModelArray` and its companion
python package [ConFixel](https://github.com/PennLINC/ConFixel). The
most important steps for installing `ModelArray` are:

-   Make sure you have necessary libraries for HDF5 - see [this
    section](https://pennlinc.github.io/ModelArray/articles/a01_installations.html#install-hdf5-libraries-in-the-system)
-   Install `ModelArray` from GitHub - see [this
    section](https://pennlinc.github.io/ModelArray/articles/a01_installations.html#install-modelarray-r-package-from-github)

## How to use

Load the `ModelArray` package into R via:

``` r
library(ModelArray)
```

We provide a example walkthrough
[here](https://pennlinc.github.io/ModelArray/articles/a02_walkthrough.html).

For documentations of `ModelArray` functions, you can:

-   Either go to [this
    webpage](https://pennlinc.github.io/ModelArray/reference/index.html);
-   Or in R console, type: `help(<function_name>)`. For example:
    `help(ModelArray.lm)`
