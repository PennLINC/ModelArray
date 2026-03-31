### Methods of "ModelArray" ###

### Show ModelArray #####
#' @rdname ModelArray-class
#'
#' @description
#' Prints a summary of the ModelArray including file path, scalar dimensions,
#' and any saved analysis names.
#'
#' @param object A \linkS4class{ModelArray} object.
#'
#' @return Invisible \code{NULL}. Called for its side effect of printing
#'   to the console.
#'
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
#' @description
#' Retrieve the named list of source filename vectors from a
#' \linkS4class{ModelArray}. Each element of the list corresponds to one
#' scalar and contains a character vector of filenames, one per input
#' file/subject.
#'
#' @param x A \linkS4class{ModelArray} object.
#'
#' @return A named list of character vectors. Names correspond to scalar names
#'   (e.g. \code{"FD"}, \code{"FC"}).
#'
#' @seealso \code{\link{scalars}}, \code{\link{results}},
#'   \code{\link{scalarNames}}, \code{\link{nInputFiles}}
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5", scalar_types = c("FD"))
#' sources(ma) # named list
#' sources(ma)[["FD"]] # character vector of filenames
#' }
#'
#' @name sources
#' @export
setGeneric("sources", function(x) standardGeneric("sources"))

#' @rdname sources
#' @export
setMethod("sources", "ModelArray", function(x) x@sources)


#' Element-wise scalar data of a ModelArray object
#'
#' @description
#' Retrieve scalar matrices from a \linkS4class{ModelArray}. When called
#' with no additional arguments, returns the full named list of all scalar
#' matrices. When called with a single scalar name, returns that one matrix.
#'
#' @param x A \linkS4class{ModelArray} object.
#' @param ... Optional: a single character string giving the scalar name to
#'   extract. If omitted, the entire named list is returned.
#'
#' @return If called with no extra arguments, a named list of
#'   \linkS4class{DelayedArray} matrices (elements as rows, source files as
#'   columns). If called with a scalar name, the corresponding single
#'   \linkS4class{DelayedArray} matrix.
#'
#' @seealso \code{\link{sources}}, \code{\link{results}},
#'   \code{\link{scalarNames}}, \code{\link{nElements}}
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5", scalar_types = c("FD", "FC"))
#' scalars(ma) # named list of all scalars
#' scalars(ma, "FD") # single DelayedArray matrix
#' }
#'
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
#' @description
#' Retrieve previously saved analysis results from a \linkS4class{ModelArray}.
#' When called with no additional arguments, returns the full named list of all
#' result sets. When called with a single analysis name, returns that one
#' result set.
#'
#' @details
#' Each result set is itself a list containing at minimum
#' \code{results_matrix} (a \linkS4class{DelayedArray} with elements as rows
#' and statistics as columns). Column names are stored alongside the matrix
#' in the HDF5 file.
#'
#' Results are only available if \code{analysis_names} was supplied when the
#' \linkS4class{ModelArray} was constructed, or if results have been written
#' back with \code{\link{writeResults}}.
#'
#' @param x A \linkS4class{ModelArray} object.
#' @param ... Optional: a single character string giving the analysis name to
#'   extract. If omitted, the entire named list is returned.
#'
#' @return If called with no extra arguments, a named list of result lists.
#'   If called with an analysis name, the corresponding result list (containing
#'   at minimum \code{results_matrix}).
#'
#' @seealso \code{\link{sources}}, \code{\link{scalars}},
#'   \code{\link{analysisNames}}, \code{\link{writeResults}},
#'   \code{\link{ModelArray.lm}}, \code{\link{ModelArray.gam}}
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5",
#'   scalar_types = c("FD"),
#'   analysis_names = c("lm_age")
#' )
#' results(ma) # named list of all results
#' results(ma, "lm_age") # single result set
#' results(ma, "lm_age")$results_matrix
#' }
#'
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
#' @description
#' Constructs a per-element data.frame from a \linkS4class{ModelArray} that
#' mirrors the \code{data} argument passed to user functions by
#' \code{\link{ModelArray.wrap}}. This is useful for testing and debugging
#' user-supplied functions outside of the full element-wise analysis loop.
#'
#' @details
#' Returns a copy of \code{phenotypes} with an extra column named by
#' \code{scalar} populated with the selected element's values from the
#' \linkS4class{ModelArray}. This mirrors the per-element data that
#' \code{\link{ModelArray.wrap}} passes to user functions (as \code{data = dat}).
#'
#' Use this to verify that your custom function works correctly on a single
#' element before committing to a full \code{\link{ModelArray.wrap}} run.
#'
#' @param x A \linkS4class{ModelArray} object.
#' @param ... Additional arguments passed to the method (currently unused
#'   at the generic level).
#'
#' @seealso \code{\link{ModelArray.wrap}}, \code{\link{scalars}}
#'
#' @rdname exampleElementData
#' @export
setGeneric("exampleElementData", function(x, ...) standardGeneric("exampleElementData"))

#' @rdname exampleElementData
#'
#' @param x A \linkS4class{ModelArray} object.
#' @param scalar Character. The name of the element-wise scalar to append
#'   as a column. Must be one of \code{names(scalars(x))}. Default is
#'   \code{"FD"}.
#' @param i_element Integer. The 1-based index of the element whose values
#'   should be extracted. Must be between 1 and the number of elements for
#'   the given scalar.
#' @param phenotypes A data.frame of the cohort with columns of independent
#'   variables and covariates. Must have the same number of rows as the
#'   number of source files in the \linkS4class{ModelArray}.
#'
#' @return A data.frame: the input \code{phenotypes} with one additional
#'   column named by \code{scalar} containing that element's values across
#'   all subjects.
#'
#' @examples
#' \dontrun{
#' h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
#' csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
#' ma <- ModelArray(h5_path, scalar_types = c("FD"))
#' phen <- read.csv(csv_path)
#' df <- exampleElementData(ma, scalar = "FD", i_element = 1, phenotypes = phen)
#'
#' # Now test your custom function on this single element:
#' my_fun <- function(data, ...) {
#'   mod <- lm(FD ~ age + sex, data = data)
#'   broom::tidy(mod)
#' }
#' my_fun(data = df)
#' }
#'
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
#' @description
#' Returns the number of elements (e.g., fixels or voxels) for a given scalar
#' in a \linkS4class{ModelArray}. This is the row count of the scalar matrix.
#'
#' @param x A \linkS4class{ModelArray} object.
#' @param scalar Optional character string. Name of the scalar to query.
#'   Defaults to the first scalar in \code{names(scalars(x))}.
#'
#' @return Integer. The number of elements (rows) in the scalar matrix.
#'
#' @seealso \code{\link{nInputFiles}}, \code{\link{numElementsTotal}},
#'   \code{\link{scalars}}, \code{\link{scalarNames}}
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5", scalar_types = c("FD"))
#' nElements(ma)
#' nElements(ma, "FD")
#' }
#'
#' @name nElements
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
#' @description
#' Returns the number of input files (i.e., subjects or source files) for a
#' given scalar in a \linkS4class{ModelArray}. This is the column count of the
#' scalar matrix.
#'
#' @param x A \linkS4class{ModelArray} object.
#' @param scalar Optional character string. Name of the scalar to query.
#'   Defaults to the first scalar in \code{names(scalars(x))}.
#'
#' @return Integer. The number of input files (columns) in the scalar matrix.
#'
#' @seealso \code{\link{nElements}}, \code{\link{sources}},
#'   \code{\link{scalarNames}}
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5", scalar_types = c("FD"))
#' nInputFiles(ma)
#' }
#'
#' @name nInputFiles
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
#' @description
#' Returns the names of all scalar datasets loaded in a
#' \linkS4class{ModelArray} (e.g. \code{"FD"}, \code{"FC"}, \code{"log_FC"}).
#'
#' @param x A \linkS4class{ModelArray} object.
#'
#' @return Character vector of scalar names.
#'
#' @seealso \code{\link{scalars}}, \code{\link{analysisNames}},
#'   \code{\link{nElements}}
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5", scalar_types = c("FD", "FC"))
#' scalarNames(ma) # c("FD", "FC")
#' }
#'
#' @name scalarNames
#' @export
setGeneric("scalarNames", function(x) standardGeneric("scalarNames"))

#' @rdname scalarNames
#' @export
setMethod("scalarNames", "ModelArray", function(x) {
  names(x@scalars)
})


#' Names of analyses in a ModelArray
#'
#' @description
#' Returns the names of all analysis result sets currently loaded in a
#' \linkS4class{ModelArray}. These correspond to subfolder names under
#' \code{/results/} in the HDF5 file.
#'
#' @param x A \linkS4class{ModelArray} object.
#'
#' @return Character vector of analysis names. Returns \code{character(0)}
#'   if no analyses have been loaded or saved.
#'
#' @seealso \code{\link{results}}, \code{\link{scalarNames}},
#'   \code{\link{writeResults}}
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5",
#'   scalar_types = c("FD"),
#'   analysis_names = c("lm_age")
#' )
#' analysisNames(ma) # "lm_age"
#' }
#'
#' @name analysisNames
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
#' Reads element metadata (e.g., greyordinates for CIFTI data, or fixel/voxel
#' coordinate information) from the HDF5 file if present. The function searches
#' for known metadata dataset names (\code{"greyordinates"}, \code{"fixels"},
#' \code{"voxels"}) at the top level of the HDF5 file.
#'
#' @param x A \linkS4class{ModelArray} object.
#'
#' @return A matrix or data.frame of element metadata if found, or \code{NULL}
#'   if no known metadata dataset exists in the HDF5 file.
#'
#' @seealso \code{\link{nElements}}, \code{\link{scalars}}
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5", scalar_types = c("FD"))
#' meta <- elementMetadata(ma)
#' if (!is.null(meta)) head(meta)
#' }
#'
#' @name elementMetadata
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
    if (!is.null(result)) {
      return(result)
    }
  }
  NULL
})

# setMethod("lm",
#           signature = c("formula", "ModelArray", "data.frame", "character", "integer"))
