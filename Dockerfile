# Building Docker file for ModelArray (R package) + ConVoxel (python package). Core is config.yml for ModelArray's circle.ci.
# Note: this is NOT final version for ModelArray + ConFixel (MRtrix mrconvert wasn't added)
# When update this Dockerfile, please update:
# 1. rocker/verse:<R_version>
# 2. commitSHA_confixel
# 3. commitSHA_modelarray - see towards the end of this file
# 4. ModelArray's dependent R packages - see DESCRIPTION file

# ConFixel requires `mrconvert` from MRtrix3:
FROM pennbbl/qsiprep-mrtrix3:22.1.0 as build_mrtrix3
    # ^^ pennbbl/qsiprep-mrtrix3:22.1.0: uses MRtrix SHA = 3498ff4, commited Jul 16, 2021
    # https://github.com/MRtrix3/mrtrix3/commit/3498ff469b843d5b023c3675f1d955ba4105c5d1

## Base image https://hub.docker.com/u/rocker/
FROM rocker/verse:4.1.2

## versions and parameters:
# specify the commit SHA:  # e.g. https://github.com/PennLINC/qsiprep/blob/master/Dockerfile#L174
    # should be the full SHA
ENV commitSHA_confixel="5d0e9c43ec26a29f3bd18c315e1bfb1429b872e0"

RUN mkdir /home/data
# RUN mkdir /home/ModelArray

## MRtrix3
COPY --from=build_mrtrix3 /opt/mrtrix3-latest /opt/mrtrix3-latest
ENV PATH="$PATH:/opt/mrtrix3-latest/bin" \
    MRTRIX3_DEPS="bzip2 ca-certificates curl libpng16-16 libblas3 liblapack3 libtiff5"

## Install libraries
# ref: .circleci/config.yml from ModelArray:
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    texlive-fonts-recommended \
    git \
    ${MRTRIX3_DEPS}

# Install python: # ref: https://github.com/PennLINC/flaudit/blob/master/Dockerfile#L23
RUN apt-get update && apt-get install -y python3-pip python3-dev


## Install ConFixel (python package)
RUN git clone -n https://github.com/PennLINC/ConFixel.git
WORKDIR ConFixel
RUN git checkout ${commitSHA_confixel}
RUN pip install .
WORKDIR /
# RUN rm -r ConFixel


## Install dependent R packages:
# from CRAN:   # removed base packages from the list (otherwise warning in docker build): methods and parallel
RUN install2.r --error --ncpus -4 \
    matrixStats \
    magrittr \
    dplyr \
    tidyr \
    tibble \
    stringr \
    glue \
    doParallel \
    hdf5r \
    mgcv \
    rlang \
    broom \
    pbmcapply \
    pbapply \
    crayon

# from Bioc: # first, install BiocManager from CRAN:
RUN install2.r --error BiocManager
RUN R -e 'BiocManager::install("HDF5Array")'
    # ^^ this step takes a long time - around 20-30min...
RUN R -e 'BiocManager::install("rhdf5")'
RUN R -e 'BiocManager::install("DelayedArray")'


## install ModelArray (R package)
COPY . /ModelArray
WORKDIR ModelArray
RUN R -e 'devtools::install()'
