# Binder image for vignette("rbc-s3-tiledb"): everything the example report needs.
# CRAN packages come from the Posit Package Manager snapshot named in runtime.txt.
install.packages(c("remotes", "BiocManager", "tiledb", "mgcv", "ggplot2", "gratia",
                   "gamm4", "lme4", "MASS", "rmarkdown", "knitr", "jsonlite", "cowplot", "png"))
BiocManager::install(c("TileDBArray", "HDF5Array", "rhdf5", "DelayedArray"),
                     update = FALSE, ask = FALSE)
options(repos = BiocManager::repositories())
remotes::install_github("andy1764/ComBatFamily", upgrade = "never")
# This repository is ModelArray itself, at the branch Binder was launched from.
remotes::install_local(".", upgrade = "never", force = TRUE)
