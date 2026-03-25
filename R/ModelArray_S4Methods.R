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

    src_ref <- sources(x)[[scalar[[1]]]]
    for (s in scalar) {
      src_s <- sources(x)[[s]]
      reorder_idx <- match(src_ref, src_s)
      if (any(is.na(reorder_idx))) {
        stop(paste0("sources for scalar ", s, " are not a permutation of reference scalar ", scalar[[1]]))
      }
      dat[[s]] <- scalars(x)[[s]][i_element, reorder_idx]
    }
    dat
  }
)

#' Multi-scalar example per-element data.frame for user functions
#'
#' @description
#' `exampleElementDataMultiScalar` is a helper for constructing per-element data.frames
#' when working with multiple scalars and long-format phenotypes (one row per scalar).
#' It mirrors the ID-based alignment used internally by `ModelArray.lm`, but returns the
#' per-element data.frame directly for inspection or use in custom functions.
#'
#' @param x A `ModelArray` object
#' @param scalar A character vector of one or more scalar names
#' @param i_element An integer, the i_th element (1-based)
#' @param phenotypes A data.frame of the cohort with independent variables/covariates
#' @param id_cols Optional character vector of identifier columns used to align subjects
#'   across scalars when `phenotypes` is in long format. If `NULL`, defaults to
#'   `participant_id`, `subject_id`, or `source_file` (first present).
#' @return A data.frame with additional response column(s) named by `scalar`
#' @rdname exampleElementData
#' @export
exampleElementDataMultiScalar <- function(x, scalar, i_element = 1L, phenotypes, id_cols = NULL) {
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

  scalar_ref <- scalar[[1]]

  # Resolve id_cols
  if (!is.null(id_cols)) {
    if (!is.character(id_cols) || length(id_cols) < 1L) {
      stop("id_cols must be a non-empty character vector when provided")
    }
    missing_ids <- setdiff(id_cols, colnames(phenotypes))
    if (length(missing_ids) > 0) {
      stop(paste0("id_cols not found in phenotypes: ", paste(missing_ids, collapse = ", ")))
    }
    id_cols_resolved <- id_cols
  } else {
    id_cols_resolved <- if ("participant_id" %in% colnames(phenotypes)) {
      "participant_id"
    } else if ("subject_id" %in% colnames(phenotypes)) {
      "subject_id"
    } else {
      "source_file"
    }
  }

  make_id_vector <- function(df) {
    vals <- df[, id_cols_resolved, drop = FALSE]
    vals[] <- lapply(vals, as.character)
    if (length(id_cols_resolved) == 1) {
      vals[[1]]
    } else {
      do.call(paste, c(vals, sep = "__"))
    }
  }

  # Long-format multi-scalar phenotypes
  if ("scalar_name" %in% colnames(phenotypes) && length(unique(phenotypes$scalar_name)) > 1L && length(scalar) > 1L) {
    ph_list <- lapply(scalar, function(s) subset(phenotypes, scalar_name == s))
    names(ph_list) <- scalar

    # uniqueness checks per scalar on IDs
    for (s in scalar) {
      ids_s <- make_id_vector(ph_list[[s]])
      if (anyDuplicated(ids_s)) {
        stop(paste0(
          "phenotypes has duplicate identifier entries for scalar ",
          s,
          " using id_cols: ",
          paste(id_cols_resolved, collapse = ", ")
        ))
      }
    }

    ids_by_scalar <- lapply(ph_list, make_id_vector)
    common_ids <- Reduce(intersect, ids_by_scalar)
    if (length(common_ids) == 0) {
      stop("No common ids across requested scalars; cannot align.")
    }

    # Base data from reference scalar
    ph_ref <- ph_list[[scalar_ref]]
    ids_ref <- make_id_vector(ph_ref)
    ph_idx_ref <- match(common_ids, ids_ref)
    dat <- ph_ref[ph_idx_ref, , drop = FALSE]
    row.names(dat) <- NULL

    # Attach scalar values for each requested scalar
    for (s in scalar) {
      ph_s <- ph_list[[s]]
      ids_s <- make_id_vector(ph_s)
      ph_idx_s <- match(common_ids, ids_s)
      sf_vec <- ph_s[["source_file"]][ph_idx_s]

      src_s <- sources(x)[[s]]
      col_idx <- match(sf_vec, src_s)
      if (any(is.na(col_idx))) {
        stop(paste0("Some ids for scalar ", s, " are missing in ModelArray sources"))
      }

      dat[[s]] <- scalars(x)[[s]][i_element, col_idx]
    }
    return(dat)
  }

  # Otherwise (wide phenotypes or single scalar): fall back to exampleElementData
  exampleElementData(x, scalar = scalar, i_element = i_element, phenotypes = phenotypes)
}

# # NOTE: ref: https://stackoverflow.com/questions/56560280/
# can-i-define-s4-methods-that-dispatch-on-more-than-one-argument-from-an-s3-gener
# setGeneric("lm", function(formula, fixelarray, phenotypes, scalar, idx, ...) standardGeneric("lm"),
#            signature = c(formula, fixelarray, phenotypes, scalar, idx)
#            )
# setMethod("lm",
#           signature = c("formula", "ModelArray", "data.frame", "character", "integer"))

