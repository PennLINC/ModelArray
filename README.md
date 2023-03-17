
<!-- TODO README.md is generated from README.Rmd. Please edit that file -->

# ModelArray

<!-- badges: start -->

[![CircleCI build
status](https://circleci.com/gh/PennLINC/ModelArray.svg?style=svg)](https://circleci.com/gh/PennLINC/ModelArray)
[![GitHub
Clones](https://img.shields.io/badge/dynamic/json?color=success&label=Clone&query=count&url=https://gist.githubusercontent.com/zhao-cy/374b45552335a37d6bd613359eb9bf67/raw/clone.json&logo=github)](https://github.com/MShawon/github-clone-count-badge)
[![Docker
pulls](https://img.shields.io/docker/pulls/pennlinc/modelarray_confixel.svg)](https://hub.docker.com/r/pennlinc/modelarray_confixel)
<!-- badges: end -->

`ModelArray` is an R package for statistical analysis of fixel-wise data
and beyond. Its features include:

- Easy to use: set up your statistical analysis with just several lines
  of code;
- Supporting linear and nonlinear modeling, and extensible to more
  models:
  - At present, `ModelArray` supports linear models as well as
    generalized additive models (GAMs) with and without penalized
    splines, which are particularly useful for studying nonlinear
    effects in lifespan data. `ModelArray` is also extensible to diverse
    models available in R;
- Scalable for large-scale datasets;
- Compatible with fixel-wise data and voxel-wise data.

Please cite our [NeuroImage
paper](https://doi.org/10.1016/j.neuroimage.2023.120037) if you use
`ModelArray`:

> Zhao, C., Tapera, T. M., Bagautdinova, J., Bourque, J., Covitz, S.,
> Gur, R. E., Gur, R. C., Larsen, B., Mehta, K., Meisler, S. L., Murtha,
> K., Muschelli, J., Roalf, D. R., Sydnor, V. J., Valcarcel, A. M.,
> Shinohara, R. T., Cieslak, M. & Satterthwaite, T. D. (2023).
> ModelArray: an R package for statistical analysis of fixel-wise data,
> *NeuroImage*, In press.
> <https://doi.org/10.1016/j.neuroimage.2023.120037>

## Overview

<center>

![Overview](vignettes/overview_structure.png)

</center>

ModelArray is packaged with the companion software
[ConFixel](https://github.com/PennLINC/ConFixel) for converting
fixel-wise data and voxel-wise data to the expected file format that
ModelArray uses. Specifically,
[ConFixel](https://github.com/PennLINC/ConFixel) is Python-based
command-line interface software, and it converts between the original
image format (`.mif` for fixel-wise data, NIfTI for voxel-wise data) and
the HDF5 file format (`.h5`) used for ModelArray.

<!-- if there is any changes in this overview section, please also update ConFixel's frontpage! -->

## Installation

Please refer to webpage
[Installation](https://pennlinc.github.io/ModelArray/articles/installations.html)
for a full guidance of installation of `ModelArray` and its companion
python package [ConFixel](https://github.com/PennLINC/ConFixel). The
most important steps for installing `ModelArray` are:

- Make sure you have necessary libraries for HDF5 - see [this
  section](https://pennlinc.github.io/ModelArray/articles/installations.html#install-hdf5-libraries-in-the-system)
- Install `ModelArray` from GitHub - see [this
  section](https://pennlinc.github.io/ModelArray/articles/installations.html#install-modelarray-r-package-from-github)

Additionally, we also provide a [container
image](https://hub.docker.com/r/pennlinc/modelarray_confixel) that
includes `ModelArray` and `ConFixel`. With this container image, there
is no need for the user to install `ModelArray`, `ConFixel`, and
dependent R and Python packages. Please see [this
webpage](https://pennlinc.github.io/ModelArray/articles/container.html)
for how to use this container image.

<!-- check above links work, esp those with section titles!!! -->

## How to use

Load the `ModelArray` package into R via:

``` r
library(ModelArray)
```

We provide a walkthrough
[here](https://pennlinc.github.io/ModelArray/articles/walkthrough.html)
with example fixel-wise data. For additional notes on application to
voxel-wise data, please refer to
[here](https://pennlinc.github.io/ModelArray/articles/voxel-wise_data.html).

For documentation of `ModelArray` functions, you can:

- Either go to [this
  webpage](https://pennlinc.github.io/ModelArray/reference/index.html);
- Or in R console, type: `help(<function_name>)`. For example:
  `help(ModelArray.lm)`

Full documentation of `ModelArray` can be found
[here](https://pennlinc.github.io/ModelArray/).

Source code of `ModelArray` can be found
[here](https://github.com/PennLINC/ModelArray).
