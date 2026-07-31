# Get Started

This page walks through the installation of `ModelArray` and its
dependencies. Once installed, head to
[`vignette("walkthrough")`](https://pennlinc.github.io/ModelArray/articles/walkthrough.md)
for a complete end-to-end example.

> **Alternative**: If you have difficulty installing packages (e.g., on
> HPC clusters), you can use the container image of
> `ModelArray + ModelArrayIO` instead. See
> [`vignette("container")`](https://pennlinc.github.io/ModelArray/articles/container.md)
> for details.

We’ll directly install `ModelArray` from
[GitHub](https://github.com/PennLINC/ModelArray). Please notice that,
this is only fully tested and supported on Linux Ubuntu system and macOS
system. If you are a Windows user, you will face errors when running
with more than one CPU cores (`n_cores = 2 or above`). You may install a
Linux subsystem or virtual machine on your Windows computer; or using
the container image of `ModelArray` and running on High Performance
Computing (HPC) clusters (see
[`vignette("container")`](https://pennlinc.github.io/ModelArray/articles/container.md)
for more).

We will first set up the conda environment, then install some dependent
packages and libraries, and finally install `ModelArray`.

### Set up a conda environment

We first create a conda environment `modelarray` for installing the
companion software `ModelArrayIO` etc. We’ll install python version 3.9:

``` console
foo@bar:~$ conda create --name modelarray python=3.11
foo@bar:~$ conda activate modelarray
```

### Install MRtrix (Only required for fixel-wise data)

When working with fixel-wise data (`.mif`), install MRtrix so you can
inspect and visualize MIF inputs and outputs. It can either be installed
via `conda` in this conda environment we just created, or be compiled
from source. See [MRtrix’s webpage](https://www.mrtrix.org/download/)
for more. Type `mrview` in the terminal to check whether MRtrix
installation is successful.

If your input data is voxel-wise data, you can skip this step.

### Install HDF5 libraries in the system

ModelArray always supports data backed by the Hierarchical Data Format 5
(HDF5) file format, so we need to make sure the necessary HDF5 libraries
are installed in the system. TileDB support is optional and requires
additional R packages; see the optional TileDB section below if you plan
to use TileDB-backed stores.

#### On a Linux Ubuntu system

If you’re on Linux Ubuntu system: First, please check if `libhdf5-dev`
has been installed in the system:

``` console
foo@bar:~$ ldconfig -p | grep libhdf5*
```

If you got more than one line of outputs, congrats, you have
`libhdf5-dev` installed. Otherwise, please install it via:

``` console
foo@bar:~$ sudo apt-get update -y
foo@bar:~$ sudo apt-get install -y libhdf5-dev
```

#### On a macOS system

Use `Homebrew` to install:

``` console
foo@bar:~$ brew install hdf5
```

For details you may refer to the webpage
[here](https://formulae.brew.sh/formula/hdf5)

### Install ModelArrayIO python package from GitHub

[ModelArrayIO](https://github.com/PennLINC/ModelArrayIO) provides the
`modelarrayio` command-line tool for converting fixel-wise data
(`.mif`), voxel-wise data (NIfTI), and CIFTI data into ModelArray
storage, and for exporting saved ModelArray results back to neuroimaging
formats. Follow the commands below to install it from GitHub:

``` console
# We first activate the conda environment we just created:
foo@bar:~$ conda activate modelarray  

# Then install ModelArrayIO:
foo@bar:~$ cd ~/myProject
foo@bar:myProject$ git clone https://github.com/PennLINC/ModelArrayIO.git
foo@bar:myProject$ cd ModelArrayIO
foo@bar:myProject$ pip install .

# You may remove the original source code if you are an end user instead of a developer:
foo@bar:myProject$ cd ..
foo@bar:myProject$ rm -r ModelArrayIO
```

After installation, conversion is done with
`modelarrayio to-modelarray`. The command auto-detects the imaging
modality from the source files in the cohort CSV. For example:

``` console
foo@bar:myProject$ modelarrayio to-modelarray --help
```

### Install R

`ModelArray` requires R \>=4.1.2, and we currently only tests ModelArray
on R 4.1.2. If you have already installed it, you can skip this step. If
you don’t, you may download it from [CRAN](https://cran.r-project.org/).

### (Optional) Install RStudio

RStudio provides a good IDE for using R. However it’s optional to
install RStudio. If you haven’t got one but want to install it, you may
download it from
[here](https://www.rstudio.com/products/rstudio/download/).

### Install ModelArray R package from GitHub

ModelArray R package’s source code is available on GitHub. To install it
in R:

``` r

# First, load library "devtools":
library(devtools)

# if you got error, it means you don't have devtools installed, please install it by: install.packages("devtools")

# Then, install ModelArray:
devtools::install_github("PennLINC/ModelArray")
```

### (Optional) Install TileDB backend dependencies

ModelArray can also read and write TileDB-backed stores. This is
optional: you only need these packages if you plan to use
`backend = "tiledb"` or paths ending in `.tdb`.

``` r

install.packages(c("jsonlite", "tiledb"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("TileDBArray")
```

You can verify TileDB support with:

``` r

library(ModelArray)

ModelArraySummary("path/to/store.tdb", backend = "tiledb")
```

If installing the `tiledb` R package fails, install the TileDB system
library first. On macOS:

``` console
foo@bar:~$ brew install tiledb
```

Now, ModelArray is ready to use:

``` r

library(ModelArray)
```

### (Optional) Other potential packages only needed for downloading demo data in `vignette("walkthrough")` page

There are several potential packages only needed for downloading demo
data in
[`vignette("walkthrough")`](https://pennlinc.github.io/ModelArray/articles/walkthrough.md)
page. You may skip this step if you will not use the demo data.

- wget
- tar

If you don’t have them, please install them first. For macOS system, you
may try out `brew` to install them. For Linux Ubuntu system, you may try
out `sudo apt-get` to install them.
