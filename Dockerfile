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
ENV commitSHA_confixel="ec7ad92a51fa80a484b80be7b83e5a8067ed9418"


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

## Add metadata:
ARG BUILD_DATE
ARG VCS_REF
#ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="modelarray_confixel" \
      org.label-schema.description="ModelArray - an R package for statistical analysis of fixel-wise data" \
      org.label-schema.url="https://pennlinc.github.io/ModelArray/" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/PennLINC/ModelArray" \
      #org.label-schema.version=$VERSION \
        # ^^ I did not add this, because users should check out version by `packageVersion("ModelArray")`
        # in R when running this Docker image
        # also, it's a bit hard to get this version in circleci (as the base image of docker building does not have R...)
        # but someone says it is "git branch name"?? ref: https://guide.opencord.org/cord-5.0/build_images.html
      org.label-schema.schema-version="1.0"
# ^^these information can be viewed by:
    # docker inspect pennlinc/modelarray_confixel:<docker_tag>