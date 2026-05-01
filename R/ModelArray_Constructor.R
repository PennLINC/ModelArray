# ModelArray class ----
#' ModelArray class
#'
#' ModelArray is an S4 class that represents element-wise scalar data and
#' associated statistical results backed by an HDF5 file on disk.
#'
#' @description
#' A ModelArray wraps one or more element-wise scalar matrices (e.g., FD, FC,
#' log_FC for fixel data) read lazily via \pkg{DelayedArray}, along with any
#' previously saved analysis results. The object holds references to the
#' underlying HDF5 file and reads data on demand, making it suitable for
#' large-scale neuroimaging datasets.
#'
#' @details
#' Each scalar in the HDF5 file is stored at \code{/scalars/<name>/values}
#' as a matrix of elements (rows) by source files (columns). Source filenames
#' are read from HDF5 attributes or companion datasets. Analysis results, if
#' present, live under \code{/results/<analysis_name>/results_matrix}.
#'
#' ModelArray objects are typically created with the \code{\link{ModelArray}}
#' constructor function. Element-wise models are fit with
#' \code{\link{ModelArray.lm}}, \code{\link{ModelArray.gam}}, or
#' \code{\link{ModelArray.wrap}}.
#'
#' @slot sources A named list of character vectors. Each element corresponds
#'   to a scalar and contains the source filenames (one per input file/subject).
#' @slot scalars A named list of [DelayedArray::DelayedArray][DelayedArray-class] matrices.
#'   Each matrix has elements as rows and source files as columns.
#' @slot results A named list of analysis results. Each element is itself a
#'   list containing at minimum \code{results_matrix} (a
#'   [DelayedArray::DelayedArray][DelayedArray-class]).
#' @slot path Character. Path(s) to the HDF5 file(s) on disk.
#'
#' @seealso \code{\link{ModelArray}} for the constructor,
#'   \code{\link{ModelArray.lm}}, \code{\link{ModelArray.gam}},
#'   \code{\link{ModelArray.wrap}} for analysis,
#'   \code{\link{scalars}}, \code{\link{sources}}, \code{\link{results}} for
#'   accessors.
#'
#' @name ModelArray-class
#' @aliases ModelArray-class
#' @rdname ModelArray-class
#' @exportClass ModelArray
ModelArray <- setClass(
  "ModelArray",
  slots = c(
    results = "list",
    sources = "list",
    scalars = "list",
    path = "character"
  )
)

#' ModelArraySeed
#'
#' Generates a "seed" for the h5 file format. A wrapper around HDF5ArraySeed
#' used to instantiate a delayed array
#'
#' @param filepath Path to an existing h5 file.
#' @param name Name of the group/field in the h5 file.
#' @param type Type of DelayedArray object, used as an argument for `HDF5Array::HDF5ArraySeed`.
#' @noRd
#'
ModelArraySeed <- function(filepath, name, type = NA) {
  # NOTE: the checker for if h5 groups fixels/voxels/scalars exist
  # (a.k.a valid fixel-wise data) is deleted, as ModelArray is generalized to any modality.

  seed <- HDF5Array::HDF5ArraySeed(filepath, name = name, type = type) # HDF5Array is also from BioConductor...

  seed
}


#' Load element-wise data from an HDF5 file
#'
#' @description
#' Reads scalar matrices and (optionally) saved analysis results from
#' an HDF5 file and returns a \linkS4class{ModelArray} object.
#'
#' @details
#' The constructor reads each scalar listed in \code{scalar_types} from
#' \code{/scalars/<scalar_type>/values}, wrapping them as
#' [DelayedArray::DelayedArray][DelayedArray-class] objects. Source filenames are extracted
#' from HDF5 attributes or companion datasets.
#'
#' If \code{analysis_names} is non-empty, saved results are loaded from
#' \code{/results/<name>/results_matrix}.
#'
#' \strong{Debugging tip:} If you encounter
#' \code{"error in evaluating the argument 'seed'..."}, check that
#' \code{scalar_types} matches groups in the file. Inspect with
#' \code{rhdf5::h5ls(filepath)}.
#'
#' @param filepath Character. Path to an existing HDF5 (\code{.h5})
#'   file containing element-wise scalar data.
#' @param scalar_types Character vector. Names of scalar groups to read
#'   from \code{/scalars/} in the HDF5 file. Default is \code{c("FD")}.
#'   Must match group names in the file.
#' @param analysis_names Character vector. Subfolder names under
#'   \code{/results/} to load. Default is \code{character(0)} (none).
#'
#' @return A \linkS4class{ModelArray} object.
#'
#' @seealso \linkS4class{ModelArray} for the class definition,
#'   \code{\link{h5summary}} for inspecting an HDF5 file.
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("path/to/data.h5", scalar_types = c("FD"))
#' ma
#' }
#'
#' @rdname ModelArray
#' @aliases ModelArray
#' @export
ModelArray <- function(filepath,
                       scalar_types = c("FD"),
                       analysis_names = character(0)) {
  # TODO: try and use hdf5r instead of rhdf5 and delayedarray here
  # fn.h5 <- H5File$new(filepath, mode="a")
  # open; "a": creates a new file or opens an existing one for read/write
  # fixel_data <- fn.h5[["fixels"]]
  # NOTE: without DelayedArray (Bioconductor),
  # the fixel_data won't look like a regular matrix in R or get transposed;
  # NOTE: I also need to test if only using hdf5r can still extract scalars(modelarray)[["FD"]]

  # TODO: IN THE FUTURE, THE SCALAR_TYPES AND ANALYSIS_NAMES ARE AUTOMATICALLY DETECTED
  # (at least detect + provide some options)

  ## scalar_data:
  sources <- vector("list", length(scalar_types))
  scalar_data <- vector("list", length(scalar_types))

  for (x in seq_along(scalar_types)) {
    # TODO: IT'S BETTER TO CHECK IF THIS SCALAR_TYPE EXISTS OR NOT..... - Chenying

    # /scalars/<scalar_type>/values:
    scalar_data[[x]] <- ModelArraySeed(
      filepath,
      name = sprintf("scalars/%s/values", scalar_types[x]),
      type = NA
    ) %>% DelayedArray::DelayedArray()

    # load source filenames (column_names): prefer attribute; fallback to dataset
    attrs <- rhdf5::h5readAttributes(filepath, name = sprintf("scalars/%s/values", scalar_types[x]))
    colnames_attr <- attrs$column_names
    if (is.null(colnames_attr)) {
      # Fallback: attempt to read from dataset-based column names
      # Try multiple plausible locations for compatibility across writers
      paths_to_try <- c(
        sprintf("scalars/%s/column_names", scalar_types[x]),
        sprintf("scalars/%s/values/column_names", scalar_types[x]),
        sprintf("scalars/scalars/%s/values/column_names", scalar_types[x]),
        sprintf("scalars/scalars/%s/column_names", scalar_types[x])
      )

      colnames_ds <- NULL
      last_error <- NULL
      for (p in paths_to_try) {
        tmp <- tryCatch(
          {
            rhdf5::h5read(filepath, p)
          },
          error = function(e) {
            last_error <<- e
            NULL
          }
        )
        if (!is.null(tmp)) {
          colnames_ds <- tmp
          if (grepl("^scalars/scalars/", p)) {
            warning(
              "Column names found at nested path '", p, "'. ",
              "This is a known quirk from some converters (e.g., concifti).",
              call. = FALSE
            )
          }
          break
        }
      }
      if (is.null(colnames_ds)) {
        stop(paste0(
          "Neither attribute 'column_names' nor a dataset with column names found. Tried: ",
          paste(paths_to_try, collapse = ", "),
          if (!is.null(last_error)) paste0(". Last error: ", conditionMessage(last_error)) else ""
        ))
      }
      # Ensure character vector, not list/matrix; trim potential null terminators and whitespace
      if (is.list(colnames_ds)) {
        colnames_ds <- unlist(colnames_ds, use.names = FALSE)
      }
      colnames_ds <- as.vector(colnames_ds)
      colnames_ds <- as.character(colnames_ds)
      # Trim any trailing NULs (hex 00) and surrounding whitespace for cross-language string compatibility
      # Use escaped hex in pattern to avoid embedding a NUL in the source code
      colnames_ds <- gsub("[\\x00]+$", "", colnames_ds, perl = TRUE, useBytes = TRUE)
      colnames_ds <- trimws(colnames_ds)
      sources[[x]] <- colnames_ds
    } else {
      sources[[x]] <- as.character(colnames_attr)
    }

    # transpose scalar_data[[x]] if needed:
    if (dim(scalar_data[[x]])[2] == length(sources[[x]])) {
      # do nothing
    } else if (dim(scalar_data[[x]])[1] == length(sources[[x]])) {
      scalar_data[[x]] <- t(scalar_data[[x]])
    } else {
      stop(
        paste0(
          "the dimension of scalar_data[[",
          toString(x),
          "]] does not match to length of sources[[",
          toString(x),
          "]]"
        )
      )
    }

    # add sources as colnames:
    colnames(scalar_data[[x]]) <- sources[[x]]
  }

  names(scalar_data) <- scalar_types
  names(sources) <- scalar_types


  ## results:
  if (length(analysis_names) == 0) {
    # user did not request any analyses; do not touch /results
    results_data <- list()
  } else {
    # user requested analyses; check if results group exists in this .h5 file
    flag_results_exist <- flagResultsGroupExistInh5(filepath)
    # message(flag_results_exist)
    if (flag_results_exist == FALSE) {
      results_data <- list()
    } else {
      # results group exist --> to load subfolders
      results_data <- vector("list", length(analysis_names))

      for (x in seq_along(analysis_names)) {
        analysis_name <- analysis_names[x]

        # we need to check if this subfolder exists in this .h5 file:
        flag_analysis_exist <- flagAnalysisExistInh5(filepath, analysis_name = analysis_name)
        if (flag_analysis_exist == FALSE) {
          stop(paste0("This analysis: ", analysis_name, " does not exist..."))
        } else {
          # exists
          # Load column names for results: prefer attribute; fallback to dataset
          attrs <- rhdf5::h5readAttributes(filepath,
            name = sprintf("results/%s/results_matrix", analysis_name)
          )
          names_results_matrix <- attrs$colnames
          if (is.null(names_results_matrix)) {
            # Fallback to dataset-based column names (similar to scalar handling)
            paths_to_try <- c(
              sprintf("results/%s/column_names", analysis_name),
              sprintf("results/%s/results_matrix/column_names", analysis_name)
            )
            colnames_ds <- NULL
            last_error <- NULL
            for (p in paths_to_try) {
              tmp <- tryCatch(
                {
                  rhdf5::h5read(filepath, p)
                },
                error = function(e) {
                  last_error <<- e
                  NULL
                }
              )
              if (!is.null(tmp)) {
                colnames_ds <- tmp
                break
              }
            }
            if (is.null(colnames_ds)) {
              stop(paste0(
                "Neither attribute 'colnames' nor a dataset with column names found for results. Tried: ",
                paste(paths_to_try, collapse = ", "),
                if (!is.null(last_error)) paste0(". Last error: ", conditionMessage(last_error)) else ""
              ))
            }
            if (is.list(colnames_ds)) {
              colnames_ds <- unlist(colnames_ds, use.names = FALSE)
            }
            colnames_ds <- as.vector(colnames_ds)
            colnames_ds <- as.character(colnames_ds)
            # Trim trailing NULs and whitespace
            colnames_ds <- gsub("[\\x00]+$", "", colnames_ds, perl = TRUE, useBytes = TRUE)
            colnames_ds <- trimws(colnames_ds)
            names_results_matrix <- colnames_ds
          }

          # names_results_matrix <- ModelArraySeed(filepath, name = sprintf(
          #   "results/%s/has_names", analysis_name), type = NA) %>%
          #   DelayedArray::DelayedArray()
          # if (dim(names_results_matrix)[1]<dim(names_results_matrix[2]){
          #   names_results_matrix <- t(names_results_matrix)
          # }

          # /results/<analysis_name>/results_matrix:
          results_data[[x]]$results_matrix <- ModelArraySeed(
            filepath,
            name = sprintf("results/%s/results_matrix", analysis_name),
            type = NA
          ) %>% DelayedArray::DelayedArray()

          if (dim(results_data[[x]]$results_matrix)[2] != length(names_results_matrix)) {
            # transpose if needed
            results_data[[x]]$results_matrix <- t(results_data[[x]]$results_matrix)
          }

          # designate the column names
          colnames(results_data[[x]]$results_matrix) <- as.character(
            names_results_matrix
          )


          # /results/<analysis_name>/lut_col?:   # LOOP OVER # OF COL OF $RESULTS_MATRIX, AND SEE IF THERE IS LUT_COL
          for (i_col in seq_along(names_results_matrix)) {
            object_name <- paste0("lut_forcol", as.character(i_col))
            flag_lut_exist <- flagObjectExistInh5(
              filepath,
              group_name = paste0("/results/", analysis_name),
              object_name = object_name
            )
            if (flag_lut_exist == TRUE) {
              lut <- ModelArraySeed(
                filepath,
                name = paste0("results/", analysis_name, "/", object_name),
                type = NA
              ) %>% DelayedArray::DelayedArray()

              # results_data[[x]]$lut[[i_col]] <- lut

              # turn values in results_matrix into factors |
              # HOWEVER, this also makes the entire $results_matrix into type "character"....
              lut <- lut %>% as.character()
              for (j_lut in seq_along(lut)) {
                str_lut <- lut[j_lut]
                idx_list <- results_data[[x]]$results_matrix[, i_col] %in% c(j_lut)
                results_data[[x]]$results_matrix[idx_list, i_col] <- lut[j_lut]
              }

              # } else {  # the lut for this column does not exist
              #   results_data[[x]]$lut[[i_col]] <- NULL
            }
          }

          # name the analysis:
          names(results_data)[[x]] <- analysis_name


          # NOTES:
          # if there is no "$lut", we can remove "$results_matrix", so that results(ModelArray)
          # would look like: $<myAnalysis>, instead of $<myAnalysis>$results_matrix
        }
      }
    }
  }


  new(
    "ModelArray",
    sources = sources,
    scalars = scalar_data,
    results = results_data,
    # TODO: issue: LHS SHOULD BE THE SAME AS THE NAME IN THE H5 FILE, NOT NECESSARY CALLED "results"
    path = filepath
  )
}


#' Number of elements in ModelArray
#'
#' @description
#' Returns the number of elements in ModelArray, for a specific scalar
#'
#' @param modelarray ModelArray class
#' @param scalar_name A character, the scalar name (one of the existing scalar in \code{modelarray})
#' @return number of elements in ModelArray, for this specific scalar
#' @export
numElementsTotal <- function(modelarray, scalar_name = "FD") {
  if (!(scalar_name %in% names(scalars(modelarray)))) {
    # not an existing scalar
    stop("scalar_name requested in not in modelarray! Please check out: scalars(modelarray)")
  }

  numElementsTotal <- nrow(scalars(modelarray)[[scalar_name]])

  numElementsTotal
}

#' Fit a linear model for a single element
#'
#' @description
#' Fits \code{\link[stats]{lm}} on one element's data. When a precomputed
#' context (\code{ctx}) is provided, all loop-invariant work (formula parsing,
#' collision checks, source alignment) is skipped. When \code{ctx} is
#' \code{NULL}, the function falls back to computing everything inline
#' (legacy behaviour for direct calls / debugging).
#'
#' If the number of subjects with finite scalar values (not \code{NaN},
#' \code{NA}, or \code{Inf}) does not exceed \code{num.subj.lthr}, the
#' element is skipped and all statistics are set to \code{NaN}.
#'
#' @param i_element Integer. The 1-based index of the element to analyse.
#' @param formula A \code{\link[stats]{formula}} passed to
#'   \code{\link[stats]{lm}}. Ignored when \code{ctx} is provided (the
#'   formula is taken from the context).
#' @param modelarray A \linkS4class{ModelArray} object. Ignored when
#'   \code{ctx} is provided.
#' @param phenotypes A data.frame of the cohort with columns of independent
#'   variables and covariates. Must contain a \code{"source_file"} column
#'   matching \code{sources(modelarray)[[scalar]]}. Ignored when \code{ctx} is provided.
#' @param scalar Character. The name of the element-wise scalar to analyse.
#'   Must be one of \code{names(scalars(modelarray))}. Ignored when \code{ctx} is provided.
#' @param var.terms Character vector. Statistics to extract per term from
#'   \code{\link[broom]{tidy.lm}} (e.g. \code{"estimate"}, \code{"statistic"},
#'   \code{"p.value"}).
#' @param var.model Character vector. Statistics to extract for the overall
#'   model from \code{\link[broom]{glance.lm}} (e.g. \code{"adj.r.squared"},
#'   \code{"p.value"}).
#' @param num.subj.lthr Numeric. The pre-computed minimum number of subjects
#'   with finite values required for this element to be analysed. Elements
#'   below this threshold are skipped. This value is typically computed by
#'   the parent function from \code{num.subj.lthr.abs} and
#'   \code{num.subj.lthr.rel}.
#' @param num.stat.output  Integer or \code{NULL}. The total number of output
#'   columns (including \code{element_id}). Used when
#'   \code{flag_initiate = FALSE} to generate an all-\code{NaN} row for
#'   skipped elements. Must be \code{NULL} when \code{flag_initiate = TRUE}.
#' @param flag_initiate Logical. If \code{TRUE}, fit the model once and return
#'   metadata for initialising the output data.frame (column names and term
#'   names). If \code{FALSE}, return a numeric vector of results for this
#'   element.
#' @param on_error Character. One of \code{"stop"}, \code{"skip"}, or
#'   \code{"debug"}. When an error occurs fitting the model: \code{"stop"}
#'   halts execution; \code{"skip"} returns all-\code{NaN} for this element;
#'   \code{"debug"} drops into \code{\link{browser}} (if interactive) then
#'   skips. Default: \code{"stop"}.
#' @param ctx A precomputed context list from \code{.build_lm_context()},
#'   or \code{NULL} for legacy inline computation.
#' @param ... Additional arguments passed to \code{\link[stats]{lm}}.
#'
#' @return If \code{flag_initiate = TRUE}, a list with components:
#'   \describe{
#'     \item{column_names}{Character vector. The column names for the output
#'       data.frame, with \code{"element_id"} first.}
#'     \item{list.terms}{Character vector. The names of the model terms
#'       (from \code{\link[broom]{tidy.lm}}).}
#'   }
#'   If \code{flag_initiate = FALSE}, a numeric vector of length
#'   \code{num.stat.output} with \code{element_id} (0-based) as the first
#'   value and the requested statistics in subsequent positions. All-\code{NaN}
#'   (except \code{element_id}) if the element had insufficient valid subjects
#'   or if an error occurred with \code{on_error = "skip"}.
#'
#' @seealso \code{\link{ModelArray.lm}}, \code{.build_lm_context},
#'   \code{\link{analyseOneElement.gam}} for the GAM equivalent,
#'   \code{\link{analyseOneElement.wrap}} for user-supplied functions
#'
#' @keywords internal
#' @rdname analyseOneElement.lm
#' @export
analyseOneElement.lm <- function(i_element,
                                 formula = NULL,
                                 modelarray = NULL,
                                 phenotypes = NULL,
                                 scalar = NULL,
                                 var.terms = c("estimate", "statistic", "p.value"),
                                 var.model = c("adj.r.squared", "p.value"),
                                 num.subj.lthr,
                                 num.stat.output = NULL,
                                 flag_initiate = FALSE,
                                 on_error = "stop",
                                 ctx = NULL,
                                 ...) {


  # Resolve context ----
  # Use precomputed if available, else build inline

  if (!is.null(ctx)) {
    # Fast path: all invariant work already done
    effective_formula   <- ctx$formula
    effective_scalar    <- ctx$scalar
  } else {
    # Legacy path: compute everything inline (for direct calls / debugging)
    effective_formula   <- formula
    effective_scalar    <- scalar
    ctx <- .build_lm_context(formula, modelarray, phenotypes, scalar)
  }


  # Per-element data assembly ----
  elem <- .assemble_element_data(i_element, ctx, num.subj.lthr)

  if (!elem$sufficient) {
    if (flag_initiate) {
      return(list(column_names = NaN, list.terms = NaN))
    } else {
      onerow <- c(i_element - 1, cheapr::rep_len_(NaN, num.stat.output - 1L))
      return(onerow)
    }
  }

  dat <- elem$dat


  # Fit model ----
  arguments_lm <- list(...)
  arguments_lm$formula <- effective_formula
  arguments_lm$data <- dat

  onemodel <- tryCatch(
    {
      do.call(stats::lm, arguments_lm)
    },
    error = function(e) {
      msg <- paste0(
        "analyseOneElement.lm error at element ",
        i_element, ": ", conditionMessage(e)
      )
      if (on_error == "debug" && interactive()) {
        message(msg)
        browser()
      }
      if (on_error == "skip" || on_error == "debug") {
        warning(msg)
        if (flag_initiate) {
          return(structure(list(.lm_error_initiate = TRUE), class = "lm_error"))
        } else {
          return(structure(list(.lm_error_runtime = TRUE), class = "lm_error"))
        }
      }
      stop(e)
    }
  )

  if (inherits(onemodel, "lm_error")) {
    if (flag_initiate) {
      return(list(column_names = NaN, list.terms = NaN))
    } else {
      onerow <- c(i_element - 1, cheapr::rep_len_(NaN, num.stat.output - 1L))
      return(onerow)
    }
  }


  # Extract statistics ----
  onemodel.tidy <- onemodel %>% broom::tidy()
  onemodel.glance <- onemodel %>% broom::glance()

  # Normalize intercept term name to match expected output convention
  onemodel.tidy$term[onemodel.tidy$term == "(Intercept)"] <- "Intercept"

  list.terms <- onemodel.tidy$term

  ## Term-level statistics ----
  onerow.terms <- tibble::tibble()
  for (i_term in seq_along(list.terms)) {
    for (i_var in seq_along(var.terms)) {
      temp <- tibble::tibble(
        placeholder = onemodel.tidy[[var.terms[i_var]]][i_term]
      )
      colnames(temp) <- paste0(list.terms[i_term], ".", var.terms[i_var])
      onerow.terms <- bind_cols_check_emptyTibble(onerow.terms, temp)
    }
  }

  ## Model-level statistics ----
  onerow.model <- tibble::tibble()
  for (i_var in seq_along(var.model)) {
    temp <- tibble::tibble(placeholder = onemodel.glance[[var.model[i_var]]])
    colnames(temp) <- paste0("model.", var.model[i_var])
    onerow.model <- bind_cols_check_emptyTibble(onerow.model, temp)
  }

  onerow.all <- bind_cols_check_emptyTibble(onerow.terms, onerow.model)
  column_names <- c("element_id", colnames(onerow.all))

  if (flag_initiate) {
    toreturn <- list(
      column_names = column_names,
      list.terms = list.terms
    )
    return(toreturn)
  } else {
    onerow <- c(i_element - 1, as.numeric(onerow.all))
    return(onerow)
  }
}

#' Fit a GAM for a single element
#'
#' #' @description
#' Returns metadata (column names, smooth term names, parametric term names,
#' and the smoothing parameter criterion attribute name) used by
#' \code{\link{ModelArray.gam}} to initialise the output data.frame. When
#' \code{flag_initiate = FALSE}, it returns a numeric vector representing one
#' row of the final results matrix.
#'
#' If the number of subjects with finite scalar values does not exceed
#' \code{num.subj.lthr}, the element is skipped and all statistics are set
#' to \code{NaN}.
#'
#' @inheritParams analyseOneElement.lm
#'
#' @param formula A \code{\link[stats]{formula}} passed to
#'   \code{\link[mgcv]{gam}}.
#' @param var.smoothTerms Character vector. Statistics to extract for smooth
#'   terms from \code{\link[broom]{tidy.gam}} with \code{parametric = FALSE}
#'   (e.g. \code{"edf"}, \code{"ref.df"}, \code{"statistic"},
#'   \code{"p.value"}).
#' @param var.parametricTerms Character vector. Statistics to extract for
#'   parametric terms from \code{\link[broom]{tidy.gam}} with
#'   \code{parametric = TRUE} (e.g. \code{"estimate"}, \code{"std.error"},
#'   \code{"statistic"}, \code{"p.value"}).
#' @param var.model Character vector. Statistics to extract for the overall
#'   model from \code{\link[broom]{glance.gam}} and
#'   \code{\link[mgcv]{summary.gam}} (e.g. \code{"adj.r.squared"},
#'   \code{"dev.expl"}, \code{"sp.criterion"}).
#' @param flag_sse Logical. If \code{TRUE}, also compute the error sum of
#'   squares (\code{model.sse}) for the model, which is needed for
#'   partial R-squared calculations in \code{\link{ModelArray.gam}}.
#'   Default: \code{FALSE}.
#' @param ctx A precomputed context list from \code{.build_gam_context()},
#'   or \code{NULL}.
#' @param ... Additional arguments passed to \code{\link[mgcv]{gam}}.
#'
#' @return If \code{flag_initiate = TRUE}, a list with components:
#'   \describe{
#'     \item{column_names}{Character vector of output column names.}
#'     \item{list.smoothTerms}{Character vector of smooth term names.}
#'     \item{list.parametricTerms}{Character vector of parametric term names.}
#'     \item{sp.criterion.attr.name}{Character. The name attribute of the
#'       smoothing parameter selection criterion (e.g. \code{"REML"} or
#'       \code{"GCV.Cp"}).}
#'   }
#'   If \code{flag_initiate = FALSE}, a numeric vector of length
#'   \code{num.stat.output} with \code{element_id} (0-based) first and
#'   requested statistics in subsequent positions. All-\code{NaN} (except
#'   \code{element_id}) if the element was skipped.
#'
#' @seealso \code{\link{ModelArray.gam}} which calls this function iteratively,
#'   \code{\link{analyseOneElement.lm}} for the linear model equivalent,
#'   \code{\link{analyseOneElement.wrap}} for user-supplied functions.
#'
#' @keywords internal
#' @rdname analyseOneElement.gam
#' @export

analyseOneElement.gam <- function(i_element,
                                  formula = NULL,
                                  modelarray = NULL,
                                  phenotypes = NULL,
                                  scalar = NULL,
                                  var.smoothTerms = c("statistic", "p.value"),
                                  var.parametricTerms = c("estimate", "statistic", "p.value"),
                                  var.model = c("dev.expl"),
                                  num.subj.lthr,
                                  num.stat.output = NULL,
                                  flag_initiate = FALSE,
                                  flag_sse = FALSE,
                                  on_error = "stop",
                                  ctx = NULL,
                                  ...) {


  # Resolve context ----
  if (!is.null(ctx)) {
    effective_formula <- ctx$formula
    effective_scalar  <- ctx$scalar
  } else {
    effective_formula <- formula
    effective_scalar  <- scalar
    ctx <- .build_gam_context(formula, modelarray, phenotypes, scalar)
  }

  # Per-element data assembly ----
  elem <- .assemble_element_data(i_element, ctx, num.subj.lthr)

  if (!elem$sufficient) {
    if (flag_initiate) {
      return(list(
        column_names = NaN,
        list.smoothTerms = NaN,
        list.parametricTerms = NaN,
        sp.criterion.attr.name = NaN
      ))
    } else {
      onerow <- c(i_element - 1, cheapr::rep_len_(NaN, num.stat.output - 1L))
      return(onerow)
    }
  }

  dat <- elem$dat

  # Fit GAM ----
  use_G_path <- !is.null(ctx$G_template) &&
    nrow(dat) == nrow(ctx$phenotypes)
  
  arguments <- list(...)
  
  onemodel <- tryCatch(
    {
      if (use_G_path) {
        # Fast path: reuse pre-built smooth bases and penalties
        G_elem <- ctx$G_template
        G_elem$y <- dat[[as.character(effective_formula[[2]])]]
        G_elem$sp <- rep(-1, ctx$n_sp)
        
        # Pass any additional ... arguments as gam() control args
        # (e.g., method = "REML") via the G_elem object
        if (!is.null(arguments$method)) {
          G_elem$method <- arguments$method
        }
        
        mgcv::gam(G = G_elem)
      } else {
        # Standard path: full gam() call (missing subjects or no G_template)
        arguments$formula <- effective_formula
        arguments$data <- dat
        do.call(mgcv::gam, arguments)
      }
    },
    error = function(e) {
      # If G= path failed, try standard path as fallback
      if (use_G_path) {
        fallback <- tryCatch({
          args_fb <- arguments
          args_fb$formula <- effective_formula
          args_fb$data <- dat
          do.call(mgcv::gam, args_fb)
        }, error = function(e2) NULL)
        if (!is.null(fallback)) return(fallback)
      }
      
      msg <- paste0(
        "analyseOneElement.gam error at element ",
        i_element, ": ", conditionMessage(e)
      )
      if (on_error == "debug" && interactive()) {
        message(msg)
        browser()
      }
      if (on_error == "skip" || on_error == "debug") {
        warning(msg)
        if (flag_initiate) {
          return(structure(list(.gam_error_initiate = TRUE), class = "gam_error"))
        } else {
          return(structure(list(.gam_error_runtime = TRUE), class = "gam_error"))
        }
      }
      stop(e)
    }
  )

  if (inherits(onemodel, "gam_error")) {
    if (flag_initiate) {
      return(list(
        column_names = NaN,
        list.smoothTerms = NaN,
        list.parametricTerms = NaN,
        sp.criterion.attr.name = NaN
      ))
    } else {
      onerow <- c(i_element - 1, cheapr::rep_len_(NaN, num.stat.output - 1L))
      return(onerow)
    }
  }


  # Extract statistics ----

  onemodel.tidy.smoothTerms <- onemodel %>% broom::tidy(parametric = FALSE)
  onemodel.tidy.parametricTerms <- onemodel %>% broom::tidy(parametric = TRUE)
  onemodel.glance <- onemodel %>% broom::glance()
  onemodel.summary <- onemodel %>% summary()

  # ---- Normalize term names (must happen before column name assembly) ----

  # Parametric: (Intercept) → Intercept
  if (nrow(onemodel.tidy.parametricTerms) > 0) {
    onemodel.tidy.parametricTerms$term[
      onemodel.tidy.parametricTerms$term == "(Intercept)"
    ] <- "Intercept"
  }

  # Smooth: s(age) → s_age, ti(x,z) → ti_x_z, s(x):oFactor → s_x_BYoFactor
  num.smoothTerms <- onemodel.summary$m
  if (num.smoothTerms > 0) {
    if (nrow(onemodel.tidy.smoothTerms) > 0) {
      onemodel.tidy.smoothTerms$term[
        onemodel.tidy.smoothTerms$term == "(Intercept)"
      ] <- "Intercept"
    }

    for (i_row in seq_len(nrow(onemodel.tidy.smoothTerms))) {
      term_name <- onemodel.tidy.smoothTerms$term[i_row]
      str_list <- strsplit(term_name, split = "[()]")[[1]]

      str <- str_list[2]
      smooth_name <- str_list[1]
      str_valid <- paste0(smooth_name, "_", str)

      if (length(str_list) > 2) {
        str_valid <- paste0(
          str_valid, "_",
          paste(str_list[3:length(str_list)], collapse = "")
        )
      }

      # : → BY
      str_valid <- gsub(":", "BY", str_valid, fixed = TRUE)
      # , → _
      str_valid <- gsub(",", "_", str_valid, fixed = TRUE)

      onemodel.tidy.smoothTerms$term[i_row] <- str_valid
    }
  }

  list.smoothTerms <- onemodel.tidy.smoothTerms$term
  list.parametricTerms <- onemodel.tidy.parametricTerms$term

  ## Smooth term statistics ----
  onerow.smoothTerms <- tibble::tibble()
  for (i_term in seq_along(list.smoothTerms)) {
    for (i_var in seq_along(var.smoothTerms)) {
      temp <- tibble::tibble(
        placeholder = onemodel.tidy.smoothTerms[[var.smoothTerms[i_var]]][i_term]
      )
      colnames(temp) <- paste0(list.smoothTerms[i_term], ".", var.smoothTerms[i_var])
      onerow.smoothTerms <- bind_cols_check_emptyTibble(onerow.smoothTerms, temp)
    }
  }

  ## Parametric term statistics ----
  onerow.parametricTerms <- tibble::tibble()
  for (i_term in seq_along(list.parametricTerms)) {
    for (i_var in seq_along(var.parametricTerms)) {
      temp <- tibble::tibble(
        placeholder = onemodel.tidy.parametricTerms[[var.parametricTerms[i_var]]][i_term]
      )
      colnames(temp) <- paste0(list.parametricTerms[i_term], ".", var.parametricTerms[i_var])
      onerow.parametricTerms <- bind_cols_check_emptyTibble(onerow.parametricTerms, temp)
    }
  }

  ## Model-level statistics ----
  onerow.model <- tibble::tibble()
  for (i_var in seq_along(var.model)) {
    val <- NULL
    if (var.model[i_var] %in% colnames(onemodel.glance)) {
      val <- onemodel.glance[[var.model[i_var]]]
    } else if (var.model[i_var] == "adj.r.squared") {
      val <- onemodel.summary$r.sq
    } else if (var.model[i_var] == "dev.expl") {
      val <- onemodel.summary$dev.expl
    } else if (var.model[i_var] == "sp.criterion") {
      val <- onemodel.summary$sp.criterion
      if (flag_initiate) {
        sp.criterion.attr.name <- names(onemodel.summary$sp.criterion)
      }
    } else if (var.model[i_var] == "scale") {
      val <- onemodel.summary$scale
    } else {
      val <- onemodel.glance[[var.model[i_var]]]
    }
    temp <- tibble::tibble(placeholder = val)
    colnames(temp) <- paste0("model.", var.model[i_var])
    onerow.model <- bind_cols_check_emptyTibble(onerow.model, temp)
  }

  ## SSE if requested ----
  onerow.sse <- tibble::tibble()
  if (flag_sse) {
    sse_val <- sum(onemodel$residuals^2)
    temp <- tibble::tibble(placeholder = sse_val)
    colnames(temp) <- "model.sse"
    onerow.sse <- temp
  }

  onerow.all <- bind_cols_check_emptyTibble(onerow.smoothTerms, onerow.parametricTerms)
  onerow.all <- bind_cols_check_emptyTibble(onerow.all, onerow.model)
  onerow.all <- bind_cols_check_emptyTibble(onerow.all, onerow.sse)

  column_names <- c("element_id", colnames(onerow.all))

  if (flag_initiate) {
    if (!exists("sp.criterion.attr.name")) {
      sp.criterion.attr.name <- NA_character_
    }
    return(list(
      column_names = column_names,
      list.smoothTerms = list.smoothTerms,
      list.parametricTerms = list.parametricTerms,
      sp.criterion.attr.name = sp.criterion.attr.name
    ))
  } else {
    onerow <- c(i_element - 1, as.numeric(onerow.all))
    return(onerow)
  }
}


#' Run a user-supplied function for a single element
#'
#' @description
#' Runs a user-supplied function on one element's data, preparing the
#' per-element data.frame by attaching all scalar values as new columns to
#' the provided \code{phenotypes}. This is the per-element workhorse called
#' iteratively by \code{\link{ModelArray.wrap}}.
#'
#' @details
#' Most users should call \code{\link{ModelArray.wrap}} directly, which handles
#' looping, parallelisation, and result assembly.
#' \code{analyseOneElement.wrap} is exported for advanced use cases such as
#' debugging a single element or building custom analysis loops.
#'
#' The user-supplied \code{user_fun} is called as
#' \code{user_fun(data = dat, ...)} where \code{dat} is \code{phenotypes} with
#' columns appended for \strong{all} scalars in the \linkS4class{ModelArray}
#' (not just the one named by \code{scalar}). The function should return a
#' one-row data.frame/tibble, a named list, or a named atomic vector. The
#' result is coerced to a numeric vector for assembly into the final results
#' matrix.
#'
#' Unlike \code{\link{analyseOneElement.lm}} and
#' \code{\link{analyseOneElement.gam}}, this function appends \strong{all}
#' scalar columns (not just the response) and checks for column name
#' collisions between scalar names and existing \code{phenotypes} columns.
#'
#' If the number of subjects with finite values across all scalars does not
#' exceed \code{num.subj.lthr}, the element is skipped and all statistics
#' are set to \code{NaN}.
#'
#' @inheritParams analyseOneElement.lm
#'
#' @param user_fun A function that accepts at least an argument named
#'   \code{data} (a data.frame: \code{phenotypes} with scalar columns
#'   appended for the current element) and returns one of:
#'   \itemize{
#'     \item A one-row \code{data.frame} or \code{tibble}. Multi-row
#'       returns will error.
#'     \item A named list. Unnamed lists are accepted and auto-named
#'       as \code{v1}, \code{v2}, etc.
#'     \item A named atomic vector. Unnamed vectors are accepted and
#'       auto-named as \code{v1}, \code{v2}, etc.
#'   }
#'   All return values are coerced to a numeric vector internally.
#'   Unsupported types (e.g., environments, S4 objects) will error.
#' @param ctx A precomputed context list from \code{.build_wrap_context()},
#'   or \code{NULL}.
#' @param ... Additional arguments forwarded to \code{user_fun}.
#'
#' @return If \code{flag_initiate = TRUE}, a list with one component:
#'   \describe{
#'     \item{column_names}{Character vector. The column names derived from
#'       the return value of \code{user_fun}, with \code{"element_id"}
#'       prepended.}
#'   }
#'   If \code{flag_initiate = FALSE}, a numeric vector of length
#'   \code{num.stat.output} with \code{element_id} (0-based) first and
#'   the coerced output of \code{user_fun} in subsequent positions.
#'   All-\code{NaN} (except \code{element_id}) if the element was skipped
#'   or if an error occurred with \code{on_error = "skip"}.
#'
#' @seealso \code{\link{ModelArray.wrap}} which calls this function
#'   iteratively, \code{\link{exampleElementData}} for building a test
#'   data.frame matching the format passed to \code{user_fun},
#'   \code{\link{analyseOneElement.lm}} for the linear model equivalent,
#'   \code{\link{analyseOneElement.gam}} for the GAM equivalent.
#'
#' @keywords internal
#' @rdname analyseOneElement.wrap
#' @export
analyseOneElement.wrap <- function(i_element,
                                   user_fun,
                                   modelarray = NULL,
                                   phenotypes = NULL,
                                   scalar = NULL,
                                   num.subj.lthr,
                                   num.stat.output = NULL,
                                   flag_initiate = FALSE,
                                   on_error = "stop",
                                   ctx = NULL,
                                   ...) {

  # Resolve context ----

  if (is.null(ctx)) {
    ctx <- .build_wrap_context(modelarray, phenotypes, scalar)
  }


  # Per-element data assembly ----

  elem <- .assemble_element_data(i_element, ctx, num.subj.lthr)

  if (!elem$sufficient) {
    if (flag_initiate) {
      return(list(column_names = NaN))
    } else {
      onerow <- c(i_element - 1, cheapr::rep_len_(NaN, num.stat.output - 1L))
      return(onerow)
    }
  }

  dat <- elem$dat


  # Execute user function ----

  arguments <- list(...)
  arguments$data <- dat

  result <- tryCatch(
    {
      do.call(user_fun, arguments)
    },
    error = function(e) {
      msg <- paste0(
        "analyseOneElement.wrap error at element ",
        i_element, ": ", conditionMessage(e)
      )
      if (on_error == "debug" && interactive()) {
        message(msg)
        browser()
      }
      if (on_error == "skip" || on_error == "debug") {
        warning(msg)
        if (flag_initiate) {
          return(structure(list(.wrap_error_initiate = TRUE), class = "wrap_error"))
        } else {
          return(structure(list(.wrap_error_runtime = TRUE), class = "wrap_error"))
        }
      }
      stop(e)
    }
  )

  # Coerce result to named numeric vector ----
  if (inherits(result, "wrap_error")) {
    if (flag_initiate) {
      return(list(column_names = NaN))
    } else {
      onerow <- c(i_element - 1, cheapr::rep_len_(NaN, num.stat.output - 1L))
      return(onerow)
    }
  }

  if (is.data.frame(result) || tibble::is_tibble(result)) {
    if (nrow(result) != 1) {
      msg <- paste0(
        "The user function must return a one-row data.frame/tibble, ",
        "a named list, or a named vector."
      )
      if (on_error == "skip" || on_error == "debug") {
        warning(paste0("analyseOneElement.wrap at element ", i_element, ": ", msg))
        if (flag_initiate) return(list(column_names = NaN))
        return(c(i_element - 1, cheapr::rep_len_(NaN, num.stat.output - 1L)))
      }
      stop(msg)
    }
    result_names <- colnames(result)
    result_values <- as.numeric(result[1, ])
  } else if (is.list(result)) {
    result_names <- names(result)
    if (is.null(result_names) || any(result_names == "")) {
      result_names <- paste0("v", seq_along(result))
    }
    result_values <- as.numeric(unlist(result))
  } else if (is.atomic(result)) {
    if (is.null(names(result))) {
      names(result) <- paste0("v", seq_along(result))
    }
    result_names <- names(result)
    result_values <- as.numeric(result)
  } else {
    msg <- paste0("user_fun must return a one-row data.frame/tibble, a named list, ",
                  "or a named atomic vector. Got: ", class(result)[1])
    if (on_error == "skip" || on_error == "debug") {
      warning(paste0("analyseOneElement.wrap at element ", i_element, ": ", msg))
      if (flag_initiate) return(list(column_names = NaN))
      return(c(i_element - 1, cheapr::rep_len_(NaN, num.stat.output - 1L)))
    }
    stop(msg)
  }

  column_names <- c("element_id", result_names)

  if (flag_initiate) {
    return(list(column_names = column_names))
  } else {
    onerow <- c(i_element - 1, result_values)
    return(onerow)
  }
}

#' Write outputs from element-wise statistical analysis to an HDF5 file
#'
#' @description
#' Creates a group named \code{analysis_name} under \code{/results/} in the
#' HDF5 file, then writes the statistical results data.frame (i.e. for one
#' analysis) into it as \code{results_matrix} along with column names.
#'
#' @details
#' The results are stored at
#' \code{/results/<analysis_name>/results_matrix} with column names saved
#' as a separate dataset at
#' \code{/results/<analysis_name>/column_names}.
#'
#' If any column of \code{df.output} is not numeric or integer, it is
#' coerced to numeric via \code{factor()} and the factor levels are saved
#' as a look-up table at
#' \code{/results/<analysis_name>/lut_forcol<i>}.
#'
#' \strong{Debugging tip:} If you encounter
#' \code{"Error in H5File.open(filename, mode, file_create_pl, file_access_pl)"},
#' check if the message mentions "No such file or directory". Try using an
#' absolute path for the \code{fn.output} argument.
#'
#' @param fn.output Character. The HDF5 (\code{.h5}) filename for the output.
#'   The file must already exist; use an absolute path if you encounter
#'   file-not-found errors.
#' @param df.output A data.frame of element-wise statistical results, as
#'   returned by \code{\link{ModelArray.lm}},
#'   \code{\link{ModelArray.gam}}, or \code{\link{ModelArray.wrap}}.
#'   Must inherit from \code{data.frame}.
#' @param analysis_name Character. The name for this set of results. Used
#'   as the group name under \code{/results/} in the HDF5 file.
#'   Default is \code{"myAnalysis"}.
#' @param overwrite Logical. If a group with the same \code{analysis_name}
#'   already exists in the HDF5 file, whether to overwrite it (\code{TRUE})
#'   or skip with a warning (\code{FALSE}). Default is \code{TRUE}.
#'
#' @return Invisible \code{NULL}. Called for its side effect of writing
#'   results to the HDF5 file.
#'
#' @seealso \code{\link{ModelArray.lm}}, \code{\link{ModelArray.gam}},
#'   \code{\link{ModelArray.wrap}} which produce the \code{df.output},
#'   \code{\link{results}} for reading results back from a
#'   \linkS4class{ModelArray}, \code{\link{h5summary}} for inspecting what
#'   has been written.
#'
#' @examples
#' \dontrun{
#' ma <- ModelArray("data.h5", scalar_types = c("FD"))
#' phenotypes <- read.csv("cohort.csv")
#'
#' results <- ModelArray.lm(
#'   FD ~ age + sex,
#'   data = ma,
#'   phenotypes = phenotypes,
#'   scalar = "FD"
#' )
#'
#' writeResults(
#'   fn.output = "data.h5",
#'   df.output = results,
#'   analysis_name = "lm_age_sex",
#'   overwrite = TRUE
#' )
#'
#' # Verify
#' h5summary("data.h5")
#' }
#'
#' @rdname writeResults
#' @export
writeResults <- function(fn.output,
                         df.output,
                         analysis_name = "myAnalysis",
                         overwrite = TRUE) {
  # This is enhanced version with: 1) change to hdf5r; 2) write results with only one row for one element

  # check "df.output"
  if (!("data.frame" %in% class(df.output))) {
    stop("Results dataset is not correct; must be data of type `data.frame`")
  }

  fn.output.h5 <- hdf5r::H5File$new(fn.output, mode = "a")
  # open; "a": creates a new file or opens an existing one for read/write

  # check if group "results" already exists!
  if (fn.output.h5$exists("results") == TRUE) {
    # group "results" exist
    results.grp <- fn.output.h5$open("results")
  } else {
    results.grp <- fn.output.h5$create_group("results")
  }

  # check if group "results\<analysis_name>" exists:
  exists_no_overwrite <- results.grp$exists(analysis_name) == TRUE &&
    overwrite == FALSE
  if (exists_no_overwrite) {
    warning(paste0(analysis_name, " exists but not to overwrite!"))
    # TODO: add checker for exisiting analysis_name, esp the matrix size
    results.analysis.grp <- results.grp$open(analysis_name)
    results_matrix_ds <- results.analysis.grp[["results_matrix"]]
  } else {
    # not exist; or exist && overwrite: to create
    exists_and_overwrite <- results.grp$exists(analysis_name) == TRUE &&
      overwrite == TRUE
    if (exists_and_overwrite) {
      # delete existing one first
      results.grp$link_delete(analysis_name)
      # NOTE: the file size will not shrink after your deletion..
      # this is because of HDF5, regardless of package of hdf5r or rhdf5
      # TODO: add a garbage collector after saving the results
    }

    # create:
    results.analysis.grp <- results.grp$create_group(analysis_name)
    # create a subgroup called analysis_name under results.grp

    # check "df.output": make sure all columns are floats (i.e. numeric)
    for (i_col in seq(1, ncol(df.output), by = 1)) {
      # for each column of df.output
      col_class <- as.character(sapply(df.output, class)[i_col]) # class of this column

      not_numeric_or_int <- (col_class != "numeric") &&
        (col_class != "integer")
      if (not_numeric_or_int) {
        # the column class is not numeric or integer
        message(
          paste0(
            "the column #",
            as.character(i_col),
            " of df.output to save: ",
            "data class is not numeric or integer...fixing it"
          )
        )

        # turn into numeric && write the notes in .h5 file...:
        factors <- df.output %>%
          dplyr::pull(., var = i_col) %>%
          factor()
        df.output[, i_col] <- df.output %>%
          dplyr::pull(., var = i_col) %>%
          factor() %>%
          as.numeric(.) # change into numeric of 1,2,3....

        # write a LUT for this column:
        results.analysis.grp[[paste0("lut_forcol", as.character(i_col))]] <- levels(factors)
        # save lut to .h5/results/<myAnalysis>/lut_col<?>
      }
    }

    # save:
    results.analysis.grp[["results_matrix"]] <- as.matrix(df.output)
    # results_matrix_ds <- results.analysis.grp[["results_matrix"]]   # name it

    # write column names as a dataset to avoid attribute size limits
    results.analysis.grp[["column_names"]] <- as.character(colnames(df.output))
    # NOTES: keep dataset-based names similar to scalar input handling
  }

  # # return:   # will not work if fn.output.h5$close_all()
  # output_list <- list(results.grp = results.grp,
  #                     results.analysis.grp = results.analysis.grp,
  #                     results_matrix_ds = results_matrix_ds)
  # return(output_list)

  fn.output.h5$close_all()


  # message("Results file written!")
}
