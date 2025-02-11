### Methods of "ModelArray" ###

### Show ModelArray #####
#' Show ModelArray object
#'
#' @description
#' Print the basic information for an ModelArray object, including number of source files, scalar names, and any analysis names.
#'
#' @param object An ModelArray object
#' @export
setMethod("show", "ModelArray", function(object) {
  # , group_name_results="results"
  
  # # check if there is a group of results:
  # flag_results_exist <- flagResultsExist(object, group_name_results)
  # if (flag_results_exist==TRUE) {
  #   str_results <- paste0("There is ", group_name_results, " in this ModelArray")
  # } else {
  #   str_results <- paste0("There is no ", group_name_results, " in this ModelArray")
  # }
  
  cat(
    is(object)[[1]],
    " located at ",
    object@path,
    "\n\n",
    format("  Source files:", justify = "left", width = 20),
    length(sources(object)[[1]]),
    "\n",
    # TODO: print for every scalar_name (instead of [[1]]); add "counts = "
    format("  Scalars:", justify = "left", width = 20),
    paste0(names(scalars(object)), collapse = ", "),
    "\n",
    # format("  Results:", justify = "left", width = 20), str_results, "\n",
    format("  Analyses:", justify = "left", width = 20),
    paste0(names(results(object)), collapse = ", "),
    "\n",
    sep = ""
    
  )
})

### Accessors for ModelArray #####


#' @aliases sources
setGeneric("sources", function(x)
  standardGeneric("sources"))

#' Source filenames of an ModelArray object
#'
#' @param x An ModelArray object
#' @return A list of source filenames
#' @export
setMethod("sources", "ModelArray", function(x)
  x@sources)


#' @aliases scalars
setGeneric("scalars", function(x, ...)
  standardGeneric("scalars"))

#' Element-wise scalar data of an ModelArray object
#'
#' @param x An ModelArray object
#' @param ... Additional arguments. Currently accept scalar name (a character)
#' @return A matrix of element-wise scalar data: elements (row) by source files (column).
#' @export
setMethod("scalars", "ModelArray", function(x, ...) {
  dots <- list(...)
  
  if (length(dots) == 1) {
    scalar <- dots[[1]]
    x@scalars[[scalar]]
    
  } else {
    x@scalars
    
  }
})

#' @aliases results
setGeneric("results", function(x, ...)
  standardGeneric("results"))

#' Statistical results of an ModelArray object
#'
#' @param x An ModelArray object
#' @param ... Additional arguments. Currently accept analysis name (a character)
#' @return Statistical results in this ModelArray object
#' @export
setMethod("results", "ModelArray", function(x, ...) {
  dots <- list(...)
  
  if (length(dots) == 1) {
    analysis_name <- dots[[1]]
    x@results[[analysis_name]]
    
  } else {
    x@results
    # message: if the type of $results_matrix is character, this may not be the case for all columns, but for columns that involve look-up table (lut)
  }
})

### Trying function setGeneric.... #####
# ----------------below works:
# setGeneric("lm", function(x,...) standardGeneric("lm"))
# setMethod("lm",
#           "ModelArray",    # this can be multiple classes e.g. signature(e1 = "foo", e2 = "numeric"), or signature("A1", "A2") - from: http://adv-r.had.co.nz/S4.html
#
#           function(x, ...) {
#             message("run ModelArray.lm!")
#           }
#           # function(formula, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, write = TRUE, ...){
#           #   print("xyz")
#           # }
#           # ModelArray.lm(object, ...)
#           # ModelArray.lm(formula, object, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, write = TRUE, ...)
#           )
# ----------------above works.

# # NOTE: ref: https://stackoverflow.com/questions/56560280/can-i-define-s4-methods-that-dispatch-on-more-than-one-argument-from-an-s3-gener
# setGeneric("lm", function(formula, fixelarray, phenotypes, scalar, idx, ...) standardGeneric("lm"),
#            signature = c(formula, fixelarray, phenotypes, scalar, idx)
#            )
# setMethod("lm",
#           signature = c("formula", "ModelArray", "data.frame", "character", "integer"))
