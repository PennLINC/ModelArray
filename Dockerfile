# Building Docker file for ModelArray (R package) + ConFixel (python package). Core is config.yml for ModelArray's circle.ci.
# When update this Dockerfile, please update:
# 1. tag of pre-built docker image `pennlinc/modelarray_build:<tag>`
    # Please make sure there is nothing to update in this pre-built docker image. To update that, see `ModelArray_tests` GitHub repo.
# 2. commitSHA_confixel

# Base image: using pre-built docker image:
FROM pennlinc/modelarray_build:0.0.1

## Versions and parameters:
# specify the commit SHA:  # e.g. https://github.com/PennLINC/qsiprep/blob/master/Dockerfile#L174
    # should be the full SHA
ENV commitSHA_confixel="5d0e9c43ec26a29f3bd18c315e1bfb1429b872e0"


## Install ConFixel (python package)
RUN git clone -n https://github.com/PennLINC/ConFixel.git
WORKDIR ConFixel
RUN git checkout ${commitSHA_confixel}
RUN pip install .
WORKDIR /
# RUN rm -r ConFixel


## Install ModelArray (R package)
COPY . /ModelArray
WORKDIR ModelArray
RUN R -e 'devtools::install()'
