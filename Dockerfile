FROM rocker/r2u:noble

# Install tricky bioconductor packages and minimal LaTeX for PDF generation
RUN sed -i 's#https://cloud.r-project.org#https://cran.r-project.org#' /etc/apt/sources.list.d/cran.sources \
        && apt-get update \
        && apt-get install -y --no-install-recommends \
            r-bioc-delayedarray \
            r-bioc-hdf5array \
            r-cran-broom \
            r-cran-crayon \
            r-cran-devtools \
            r-cran-doparallel \
            r-cran-dplyr \
            r-cran-glue \
            r-cran-gratia \
            r-cran-hdf5r \
            r-cran-hdf5r.extra \
            r-cran-magrittr \
            r-cran-mgcv \
            r-cran-pbapply \
            r-cran-pbmcapply \
            r-bioc-rhdf5 \
            r-cran-tibble \
            r-cran-tidyr \
            r-cran-tidyverse \
            git \
            python3-venv \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

## Install optional TileDB backend dependencies for ModelArray
RUN R -e 'install.packages(c("jsonlite", "tiledb", "BiocManager"))' \
        && R -e 'BiocManager::install("TileDBArray", ask = FALSE, update = FALSE)'

## Install ModelArrayIO (Python package and modelarrayio CLI)
RUN python3 -m venv /opt/modelarrayio \
        && /opt/modelarrayio/bin/pip install --no-cache-dir git+https://github.com/PennLINC/ModelArrayIO.git

ENV PATH="/opt/modelarrayio/bin:${PATH}"


## Install ModelArray (R package)
COPY . /ModelArray
WORKDIR /ModelArray
RUN R -e 'devtools::install()'

## Verify the R and Python TileDB paths are available
RUN R -e 'library(ModelArray); stopifnot(requireNamespace("tiledb", quietly = TRUE)); stopifnot(requireNamespace("TileDBArray", quietly = TRUE))' \
        && /opt/modelarrayio/bin/modelarrayio --version

## Add metadata:
ARG BUILD_DATE
ARG VCS_REF
#ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="modelarray_confixel" \
      org.label-schema.description="ModelArray - an R package for statistical analysis of fixel-wise data and beyond" \
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
