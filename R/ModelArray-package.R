#' ModelArray: Statistical Analysis of Element-Wise Neuroimaging Data
#'
#' @description
#' The ModelArray package provides an S4 class and associated methods for
#' performing massively univariate statistical analyses on element-wise
#' (fixel, voxel, or vertex) neuroimaging data stored in HDF5 files.
#'
#' @details
#' The core workflow is:
#' \enumerate{
#'   \item Inspect an HDF5 file with \code{\link{h5summary}()}
#'   \item Load data with \code{\link{ModelArray}()}
#'   \item Fit models with \code{\link{ModelArray.lm}},
#'     \code{\link{ModelArray.gam}}, or \code{\link{ModelArray.wrap}}
#'   \item Access results via \code{\link{scalars}}, \code{\link{results}},
#'     \code{\link{sources}}
#' }
#'
#' For multi-session or multi-modality analyses, combine ModelArrays with
#' \code{\link{mergeModelArrays}} before fitting cross-scalar models.
#'
#' @references
#' Zhao, C., et al. (2023). ModelArray--An R package for statistical
#' analysis of fixel-wise data. \emph{NeuroImage}, 271, 120037.
#' \doi{10.1016/j.neuroimage.2023.120037}
#'
#' @seealso
#' \linkS4class{ModelArray}, \code{\link{ModelArray}},
#' \code{\link{ModelArray.lm}}, \code{\link{ModelArray.gam}},
#' \code{\link{ModelArray.wrap}}, \code{\link{mergeModelArrays}},
#' \code{\link{h5summary}}
#'
#' @section Centralized imports:
#' The following imports are consolidated here because they cannot be
#' expressed as \code{pkg::fun()} qualified calls. All other add-on
#' package functions are called with explicit namespacing throughout the
#' package source code.
#'
## ---------------------------------------------------------------
## Infix operators (cannot use :: syntax)
## ---------------------------------------------------------------
#' @importFrom dplyr %>%
#'
## ---------------------------------------------------------------
## S4 infrastructure
## ---------------------------------------------------------------
#' @import methods
#' @importClassesFrom DelayedArray DelayedArray
#'
## ---------------------------------------------------------------
## S4 class constructor (needed for setClass slot type validation)
## ---------------------------------------------------------------
#' @importFrom DelayedArray DelayedArray
#'
## ---------------------------------------------------------------
## HDF5 seed constructor (no method dispatch, but called in
## ModelArraySeed which feeds into DelayedArray(); keeping here
## avoids a NOTE about undeclared dependency usage)
## ---------------------------------------------------------------
#' @importFrom HDF5Array HDF5ArraySeed
#'
#' @keywords internal
"_PACKAGE"
