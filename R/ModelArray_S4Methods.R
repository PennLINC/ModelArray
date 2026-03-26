### Methods of "ModelArray" ###

### Show ModelArray #####
#' Show ModelArray object
#'
#' @description
#' Print the basic information for an ModelArray object, including number of source files,
#' scalar names, and any analysis names.
#'
#' @param object An ModelArray object
#' @export
setMethod("show", "ModelArray", function(object) { # , group_name_results="results"

  # # check if there is a group of results:
  # flag_results_exist <- flagResultsExist(object, group_name_results)
  # if (flag_results_exist==TRUE) {
  #   str_results <- paste0("There is ", group_name_results, " in this ModelArray")
  # } else {
  #   str_results <- paste0("There is no ", group_name_results, " in this ModelArray")
  # }

  paths <- object@path
  if (length(paths) == 1) {
    path_str <- paths
  } else {
    path_str <- paste0("\n    ", paste(names(paths), paths, sep = " -> ", collapse = "\n    "))
  }
  cat(is(object)[[1]], " located at ", path_str, "\n\n", sep = "")

  scalar_names <- names(scalars(object))
  for (sn in scalar_names) {
    nr <- nrow(scalars(object)[[sn]])
    nc <- ncol(scalars(object)[[sn]])
    cat(format(paste0("  ", sn, ":"), justify = "left", width = 20),
      nr, " elements x ", nc, " input files\n",
      sep = ""
    )
  }
  analysis_names <- names(results(object))
  if (length(analysis_names) > 0) {
    cat(format("  Analyses:", justify = "left", width = 20),
      paste0(analysis_names, collapse = ", "), "\n",
      sep = ""
    )
  }
})

### Accessors for ModelArray #####


#' Source filenames of a ModelArray object
#'
#' @param x A ModelArray object
#' @return A list of source filenames
#' @name sources
#' @export
setGeneric("sources", function(x) standardGeneric("sources"))

#' @rdname sources
#' @export
setMethod("sources", "ModelArray", function(x) x@sources)


#' Element-wise scalar data of a ModelArray object
#'
#' @param x A ModelArray object
#' @param ... Additional arguments. Currently accepts a scalar name (character).
#' @return A matrix of element-wise scalar data: elements (row) by source files (column).
#' @name scalars
#' @export
setGeneric("scalars", function(x, ...) standardGeneric("scalars"))

#' @rdname scalars
#' @export
setMethod(
  "scalars",
  "ModelArray",
  function(x, ...) {
    dots <- list(...)

    if (length(dots) == 1) {
      scalar <- dots[[1]]
      x@scalars[[scalar]]
    } else {
      x@scalars
    }
  }
)

#' Statistical results of a ModelArray object
#'
#' @param x A ModelArray object
#' @param ... Additional arguments. Currently accepts an analysis name (character).
#' @return Statistical results in this ModelArray object
#' @name results
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))

#' @rdname results
#' @export
setMethod(
  "results", "ModelArray", function(x, ...) {
    dots <- list(...)

    if (length(dots) == 1) {
      analysis_name <- dots[[1]]
      x@results[[analysis_name]]
    } else {
      x@results
      # message: if the type of $results_matrix is character,
      # this may not be the case for all columns, but for columns that involve look-up table (lut)
    }
  }
)

### Trying function setGeneric.... #####
# ----------------below works:
# setGeneric("lm", function(x,...) standardGeneric("lm"))
# setMethod("lm",
#           "ModelArray",    # this can be multiple classes
#                            # e.g. signature(e1 = "foo", e2 = "numeric"),
#                            # or signature("A1", "A2") - from: http://adv-r.had.co.nz/S4.html
#
#           function(x, ...) {
#             message("run ModelArray.lm!")
#           }
#           # function(formula, data, phenotypes, scalar, verbose = TRUE,
#           # idx = NULL, pbar = TRUE, n_cores = 1, write = TRUE, ...){
#           #   print("xyz")
#           # }
#           # ModelArray.lm(object, ...)
#           # ModelArray.lm(formula, object, phenotypes, scalar, verbose = TRUE,
#           # idx = NULL, pbar = TRUE, n_cores = 1, write = TRUE, ...)
#           )
# ----------------above works.

### Example per-element data helper #####

#' Example per-element data.frame for user functions
#'
#' @title Example per-element data.frame for user functions
#' @name exampleElementData
#' @rdname exampleElementData
#' @description
#' Generic for constructing a per-element data.frame from a `ModelArray`.
#' See the `ModelArray` method for details.
#'
#' @param x A `ModelArray` object (or compatible type)
#' @param ... Additional arguments (ignored)
#' @export
setGeneric("exampleElementData", function(x, ...) standardGeneric("exampleElementData"))

#' Example per-element data.frame for user functions
#' @rdname exampleElementData
#'
#' @description
#' Returns a copy of `phenotypes` with an extra column named by `scalar` populated
#' with the selected element's values from the `ModelArray`. This mirrors the
#' per-element data that `ModelArray.wrap` passes to user functions (`data = dat`).
#'
#' @param x An ModelArray object
#' @param scalar A character. The name of the element-wise scalar to append
#' @param i_element An integer, the i_th element (1-based)
#' @param phenotypes A data.frame of the cohort with independent variables/covariates
#' @return A data.frame with the additional response column named by `scalar`
#' @examples
#' \dontrun{
#' h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
#' csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
#' ma <- ModelArray(h5_path, scalar_types = c("FD"))
#' phen <- read.csv(csv_path)
#' df <- exampleElementData(ma, scalar = "FD", i_element = 1, phenotypes = phen)
#' }
#' @export
setMethod(
  "exampleElementData",
  "ModelArray",
  function(x, scalar = "FD", i_element = 1L, phenotypes) {
    if (!is.data.frame(phenotypes)) {
      stop("phenotypes must be a data.frame")
    }
    if (!(scalar %in% names(scalars(x)))) {
      stop("scalar not found in modelarray; use one of names(scalars(x))")
    }
    num_elements <- nrow(scalars(x)[[scalar]])
    if (length(i_element) != 1L || is.na(i_element) || i_element < 1L || i_element > num_elements) {
      stop("i_element is out of range")
    }

    dat <- phenotypes
    dat[[scalar]] <- scalars(x)[[scalar]][i_element, ]
    dat
  }
)

### Convenience accessors #####

#' Number of elements in a ModelArray
#'
#' @param x A ModelArray object
#' @param scalar Optional scalar name. Defaults to the first scalar.
#' @return Integer, number of elements (rows)
#' @export
setGeneric("nElements", function(x, scalar = NULL) standardGeneric("nElements"))

#' @rdname nElements
#' @export
setMethod("nElements", "ModelArray", function(x, scalar = NULL) {
  if (is.null(scalar)) scalar <- names(x@scalars)[1]
  nrow(x@scalars[[scalar]])
})

#' Number of input files in a ModelArray
#'
#' @param x A ModelArray object
#' @param scalar Optional scalar name. Defaults to the first scalar.
#' @return Integer, number of input files (columns)
#' @export
setGeneric("nInputFiles", function(x, scalar = NULL) standardGeneric("nInputFiles"))

#' @rdname nInputFiles
#' @export
setMethod("nInputFiles", "ModelArray", function(x, scalar = NULL) {
  if (is.null(scalar)) scalar <- names(x@scalars)[1]
  ncol(x@scalars[[scalar]])
})

#' Names of scalars in a ModelArray
#'
#' @param x A ModelArray object
#' @return Character vector of scalar names
#' @export
setGeneric("scalarNames", function(x) standardGeneric("scalarNames"))

#' @rdname scalarNames
#' @export
setMethod("scalarNames", "ModelArray", function(x) {
  names(x@scalars)
})

#' Names of analyses in a ModelArray
#'
#' @param x A ModelArray object
#' @return Character vector of analysis names
#' @export
setGeneric("analysisNames", function(x) standardGeneric("analysisNames"))

#' @rdname analysisNames
#' @export
setMethod("analysisNames", "ModelArray", function(x) {
  names(x@results)
})

#' Element metadata from a ModelArray
#'
#' @description
#' Reads element metadata (e.g., greyordinates for cifti data) from the h5 file
#' if present. Returns NULL if no element metadata is found.
#'
#' @param x A ModelArray object
#' @return A matrix or data.frame of element metadata, or NULL
#' @export
setGeneric("elementMetadata", function(x) standardGeneric("elementMetadata"))

#' @rdname elementMetadata
#' @export
setMethod("elementMetadata", "ModelArray", function(x) {
  filepath <- x@path
  if (length(filepath) > 1) filepath <- filepath[1]

  # Try known metadata dataset names
  metadata_paths <- c("greyordinates", "fixels", "voxels")
  for (p in metadata_paths) {
    result <- tryCatch(
      rhdf5::h5read(filepath, p),
      error = function(e) NULL
    )
    if (!is.null(result)) return(result)
  }
  NULL
})
# setMethod("lm",
#           signature = c("formula", "ModelArray", "data.frame", "character", "integer"))
