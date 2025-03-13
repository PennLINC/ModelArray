FROM rocker/r2u:jammy

# Install tricky bioconductor packages and minimal LaTeX for PDF generation
RUN apt update \
        && apt install -y --no-install-recommends \
 		 r-cran-devtools \
         r-bioc-rhdf5 \
         r-bioc-delayedarray


## Install ModelArray (R package)
COPY . /ModelArray
WORKDIR /ModelArray
RUN R -e 'devtools::install()'

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