# Building Docker file for ModelArray (R package) + ConFixel (python package). Core is config.yml for ModelArray's circle.ci.
# When update this Dockerfile, please update:
# 1. rocker/verse:<R_version>
# 2. commitSHA_confixel
# 3. commitSHA_modelarray - see towards the end of this file
# 4. ModelArray's dependent R packages - see DESCRIPTION file

## Base image https://hub.docker.com/u/rocker/
FROM rocker/verse:4.1.2

## versions and parameters:
# specify the commit SHA:  # e.g. https://github.com/PennLINC/qsiprep/blob/master/Dockerfile#L174
ENV commitSHA_confixel="eda00bc1c9837b5586a65f24b4c986845f41582c"


RUN mkdir /home/data
# RUN mkdir /home/ModelArray

## Install libraries 
# ref: .circleci/config.yml from ModelArray:
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    texlive-fonts-recommended \
    git

# # Install libraries
# RUN apt-get update && apt-get install -y --no-install-recommends libhdf5-dev
# # Install libraries for latex
# RUN apt-get update && apt-get install -y --no-install-recommends texlive-fonts-recommended
# Install git:   # ref a bit from: https://github.com/PennLINC/qsiprep_build/blob/main/Dockerfile_DSIStudio
# RUN apt-get update && apt-get install -y --no-install-recommends git

# Install python: # ref: https://github.com/PennLINC/flaudit/blob/master/Dockerfile#L23
RUN apt-get update && apt-get install -y python3-pip python3-dev


## Install ConFixel (python package)
# git clone xxxx.git${commitSHA}
# cd to directory
# pip install .

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
RUN R -e 'BiocManager::install("rhdf5")'
RUN R -e 'BiocManager::install("DelayedArray")'

# # copy necessary files and folders:
# COPY ./install_packages.R ./install_packages.R
# # install R-packages:
# RUN Rscript ./install_packages.R


## install ModelArray (R package)
# please update commit SHA for ModelArray if needed
RUN R -e 'devtools::install_github("PennLINC/ModelArray@9e735b93f2d6b756f8ef18aadc14e2d10c6cc191")'

