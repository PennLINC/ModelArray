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

  cat(is(object)[[1]], " located at ", object@path, "\n\n",
    # TODO: print for every scalar_name (instead of [[1]]); add "counts = "
    format("  Source files:", justify = "left", width = 20), length(sources(object)[[1]]), "\n",
    format("  Scalars:", justify = "left", width = 20), paste0(names(scalars(object)), collapse = ", "), "\n",
    # format("  Results:", justify = "left", width = 20), str_results, "\n",
    format("  Analyses:", justify = "left", width = 20), paste0(names(results(object)), collapse = ", "), "\n",
    sep = ""
  )
})

### Accessors for ModelArray #####


#' @aliases sources
setGeneric("sources", function(x) standardGeneric("sources"))

#' Source filenames of an ModelArray object
#'
#' @param x An ModelArray object
#' @return A list of source filenames
#' @export
setMethod("sources", "ModelArray", function(x) x@sources)


#' @aliases scalars
setGeneric("scalars", function(x, ...) standardGeneric("scalars"))

#' Element-wise scalar data of an ModelArray object
#'
#' @param x An ModelArray object
#' @param ... Additional arguments. Currently accept scalar name (a character)
#' @return A matrix of element-wise scalar data: elements (row) by source files (column).
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

#' @aliases results
setGeneric("results", function(x, ...) standardGeneric("results"))

#' Statistical results of an ModelArray object
#'
#' @param x An ModelArray object
#' @param ... Additional arguments. Currently accept analysis name (a character)
#' @return Statistical results in this ModelArray object
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
#' @param scalar A character vector. One or more element-wise scalar names to append
#' @param i_element An integer, the i_th element (1-based)
#' @param phenotypes A data.frame of the cohort with independent variables/covariates
#' @return A data.frame with additional response column(s) named by `scalar`
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
    if (length(scalar) < 1L) {
      stop("scalar must contain at least one scalar name")
    }
    missing_scalars <- setdiff(scalar, names(scalars(x)))
    if (length(missing_scalars) > 0) {
      stop(paste0(
        "scalar not found in modelarray: ",
        paste(missing_scalars, collapse = ", "),
        "; use one of names(scalars(x))"
      ))
    }

    num_elements <- nrow(scalars(x)[[scalar[[1]]]])
    if (length(i_element) != 1L || is.na(i_element) || i_element < 1L || i_element > num_elements) {
      stop("i_element is out of range")
    }

    # If phenotypes is in long format (one row per scalar), subset to a single scalar
    if ("scalar_name" %in% colnames(phenotypes) && length(unique(phenotypes$scalar_name)) > 1) {
      dat <- subset(phenotypes, scalar_name == scalar[[1]])
    } else {
      dat <- phenotypes
    }

    # Align / validate source ordering against the ModelArray sources
    src_modelarray <- sources(x)[[scalar[[1]]]]
    src_pheno <- dat[["source_file"]]
    if (is.null(src_pheno)) {
      stop("phenotypes must contain a 'source_file' column for alignment")
    }
    if (length(src_modelarray) != length(src_pheno)) {
      stop("Length of phenotypes$source_file does not match ModelArray sources")
    }
    if (length(src_modelarray) != length(unique(src_modelarray))) {
      stop("ModelArray sources are not unique; cannot align phenotypes")
    }
    if (length(src_pheno) != length(unique(src_pheno))) {
      stop("phenotypes$source_file entries are not unique; cannot align")
    }
    if (!identical(src_modelarray, src_pheno)) {
      if ((all(src_modelarray %in% src_pheno)) && (all(src_pheno %in% src_modelarray))) {
        reorder_idx <- match(src_modelarray, src_pheno)
        dat <- dat[reorder_idx, , drop = FALSE]
        row.names(dat) <- NULL
      } else {
        stop("phenotypes$source_file entries differ from ModelArray sources; cannot align")
      }
    }

    for (s in scalar) {
      dat[[s]] <- scalars(x)[[s]][i_element, ]
    }
    dat
  }
)

# # NOTE: ref: https://stackoverflow.com/questions/56560280/
# can-i-define-s4-methods-that-dispatch-on-more-than-one-argument-from-an-s3-gener
# setGeneric("lm", function(formula, fixelarray, phenotypes, scalar, idx, ...) standardGeneric("lm"),
#            signature = c(formula, fixelarray, phenotypes, scalar, idx)
#            )
# setMethod("lm",
#           signature = c("formula", "ModelArray", "data.frame", "character", "integer"))
