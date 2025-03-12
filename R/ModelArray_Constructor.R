# Exported Functions

### setClass of "ModelArray" #####
#' An S4 class to represent element-wise scalar data and statistics.
#'
#' @slot sources A list of source filenames
#' @slot scalars A list of element-wise scalar matrix
#' @slot results A list of statistical result matrix
#' @slot path Path to the h5 file on disk
#' @importClassesFrom DelayedArray DelayedArray
ModelArray <- setClass(
  "ModelArray",
  # contains="DelayedArray",
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
#' @importFrom HDF5Array HDF5ArraySeed
#' @noRd
#'
ModelArraySeed <- function(filepath, name, type = NA) {
  # NOTE: the checker for if h5 groups fixels/voxels/scalars exist (a.k.a valid fixel-wise data)
  # is deleted, as ModelArray is generalized to any modality.

  seed <- HDF5Array::HDF5ArraySeed(
    filepath,
    name = name, type = type
  ) # HDF5Array is also from BioConductor...

  seed
}


#' Load element-wise data from .h5 file as an ModelArray object
#'
#' @details:
#' Tips for debugging:
#' if you run into this error: "Error in h(simpleError(msg, call)) :
#' error in evaluating the argument 'seed' in selecting a method for function
#' 'DelayedArray': HDF5. Symbol table. Can't open object."
#' Then please check if you give correct "scalar_types" - check via rhdf5::h5ls(filename_for_h5)
#'
#' @param filepath file
#' @param scalar_types expected scalars
#' @param analysis_names the subfolder names for results in .h5 file
#' @return ModelArray object
#' @export
#' @import methods
#' @importFrom dplyr %>%
#' @importFrom DelayedArray DelayedArray realize
#' @importFrom rhdf5 h5readAttributes
ModelArray <- function(filepath, scalar_types = c("FD"), analysis_names = c("myAnalysis")) {
  # TODO: try and use hdf5r instead of rhdf5 and delayedarray here
  # fn.h5 <- H5File$new(filepath, mode="a")    # open; "a": creates a new file or opens an existing one for read/write
  # fixel_data <- fn.h5[["fixels"]]
  # NOTE: without DelayedArray (Bioconductor), the fixel_data won't look like a regular matrix in R or get transposed;
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
      filepath, name = sprintf("scalars/%s/values", scalar_types[x]),
      type = NA
    ) %>%
      DelayedArray::DelayedArray()

    # load attribute "column_names", i.e. source filenames:
    sources[[x]] <- rhdf5::h5readAttributes(
      filepath,
      name = sprintf("scalars/%s/values", scalar_types[x])
    )$column_names %>% as.character()

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
  # first, we need to check if results group exists in this .h5 file
  flag_results_exist <- flagResultsGroupExistInh5(filepath)
  # message(flag_results_exist)
  if (flag_results_exist == FALSE) {
    results_data <- list()
  } else { # results group exist --> to load subfolders
    results_data <- vector("list", length(analysis_names))

    for (x in seq_along(analysis_names)) {
      analysis_name <- analysis_names[x]

      # we need to check if this subfolder exists in this .h5 file:
      flag_analysis_exist <- flagAnalysisExistInh5(filepath, analysis_name = analysis_name)
      if (flag_analysis_exist == FALSE) {
        stop(paste0("This analysis: ", analysis_name, " does not exist..."))
      } else { # exists
        # /results/<analysis_name>/has_names:
        names_results_matrix <- rhdf5::h5readAttributes(
          filepath,
          name = sprintf("results/%s/results_matrix", analysis_name)
        )$colnames # after updating writeResults()

        # names_results_matrix <- ModelArraySeed(
        # filepath, name = sprintf("results/%s/has_names", analysis_name), type = NA) %>%
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

        if (dim(results_data[[x]]$results_matrix)[2] != length(names_results_matrix)) { # transpose if needed
          results_data[[x]]$results_matrix <- t(results_data[[x]]$results_matrix)
        }

        colnames(results_data[[x]]$results_matrix) <- as.character(
          DelayedArray::realize(names_results_matrix)
        ) # designate the column names


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
        # if there is no "$lut", we can remove "$results_matrix",
        # so that results(ModelArray) would look like: $<myAnalysis>,
        # instead of $<myAnalysis>$results_matrix
      }
    }
  }


  new(
    "ModelArray",
    sources = sources,
    scalars = scalar_data,
    # TODO: issue: LHS SHOULD BE THE SAME AS THE NAME IN THE H5 FILE, NOT NECESSARY CALLED "results"
    results = results_data,
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
  if (!(scalar_name %in% names(scalars(modelarray)))) { # not an existing scalar
    stop("scalar_name requested in not in modelarray! Please check out: scalars(modelarray)")
  }

  numElementsTotal <- nrow(scalars(modelarray)[[scalar_name]])

  numElementsTotal
}

#' Fit linear model for one element.
#'
#' @description
#' `analyseOneElement.lm` fits a linear model for one element data, and returns requested model statistics.
#'
#' @details
#' `ModelArray.lm` iteratively calls this function to get statistics for all requested elements.
#'
#' @param i_element An integer, the i_th element, starting from 1.
#'     For initiating (flag_initiate = TRUE), use i_element=1
#' @param formula Formula (passed to `stats::lm()`)
#' @param modelarray ModelArray class
#' @param phenotypes A data.frame of the cohort with columns of independent variables and covariates
#'     to be added to the model.
#' @param scalar A character. The name of the element-wise scalar to be analysed
#' @param var.terms A list of characters. The list of variables to save for terms (got from `broom::tidy()`).
#' @param var.model A list of characters. The list of variables to save for the model (got from `broom::glance()`).
#' @param num.subj.lthr The minimal number of subjects with valid val e in input h5 file, i.e.
#'     number of subjects with finite values (defined by `is.finite()`, i.e. not NaN or NA or Inf)
#'     in h5 file > \code{num.subj.lthr}, then this element will be run normally; otherwise, this element will
#'     be skipped and statistical outputs will be set as NaN.
#' @param num.stat.output The number of output stat metrics (for generating all NaN stat when # subjects does
#'     not meet criteria). This includes column `element_id`. This is required when flag_initiate = TRUE.
#' @param flag_initiate TRUE or FALSE, Whether this is to initiate the new analysis.
#'     If TRUE, it will return column names etc to be used for initiating data.frame;
#'     if FALSE, it will return the list of requested statistic values.
#' @param ... Additional arguments for `stats::lm()`
#'
#' @return If flag_initiate==TRUE, returns column names, and list of term names of final results;
#'     if flag_initiate==FALSE, it will return the list of requested statistic values for a element.
#' @export
#' @importFrom stats lm
#' @import broom
#' @importFrom dplyr %>% select bind_cols
#' @import tibble

analyseOneElement.lm <- function(i_element,
                                 formula, modelarray, phenotypes, scalar,
                                 var.terms, var.model,
                                 num.subj.lthr, num.stat.output = NULL,
                                 flag_initiate = FALSE,
                                 ...) {
  values <- scalars(modelarray)[[scalar]][i_element, ]

  ## check number of subjects with (in)valid values:
  flag_sufficient <- NULL # whether number of subjects with valid values are sufficient
  num.subj.valid <- length(values[is.finite(values)])
  if (num.subj.valid > num.subj.lthr) {
    flag_sufficient <- TRUE
  } else {
    flag_sufficient <- FALSE
  }

  if (flag_sufficient == TRUE) {
    dat <- phenotypes
    dat[[scalar]] <- values

    # dots <- list(...)
    # dots_names <- names(dots)
    # if ("weights" %in% dots_names) {
    #   message(dots$weights)
    #   myWeights <- dots$weights
    #   dots$weights <- NULL  # remove weights from
    #
    #   arguments_lm <- dots
    #
    # }
    arguments_lm <- list(...)
    arguments_lm$formula <- formula
    arguments_lm$data <- dat

    # onemodel <- stats::lm(formula, data = dat, ...)
    # onemodel <- stats::lm(formula, data = dat, weights = myWeights,...)
    # explicitly passing arguments into lm, to avoid error of argument "weights"
    onemodel <- do.call(stats::lm, arguments_lm)

    onemodel.tidy <- onemodel %>% broom::tidy()
    onemodel.glance <- onemodel %>% broom::glance()

    # delete columns you don't want:
    var.terms.full <- names(onemodel.tidy)

    var.model.full <- names(onemodel.glance)

    # list to remove:
    var.terms.orig <- var.terms
    var.terms <- list("term", var.terms) %>% unlist() # we will always keep "term" column
    var.terms.remove <- list()
    for (l in var.terms.full) {
      if (!(l %in% var.terms)) {
        var.terms.remove <- var.terms.remove %>%
          append(., l) %>%
          unlist() # the order will still be kept
      }
    }

    var.model.remove <- list()
    for (l in var.model.full) {
      if (!(l %in% var.model)) {
        var.model.remove <- var.model.remove %>%
          append(., l) %>%
          unlist() # the order will still be kept
      }
    }

    # remove those columns:
    if (length(var.terms.remove) != 0) { # if length=0, it's list(), nothing to remove
      onemodel.tidy <- dplyr::select(onemodel.tidy, -all_of(var.terms.remove))
    }
    if (length(var.model.remove) != 0) {
      onemodel.glance <- dplyr::select(onemodel.glance, -all_of(var.model.remove))
    }

    # adjust:
    # change the term name from "(Intercept)" to "Intercept"
    onemodel.tidy$term[onemodel.tidy$term == "(Intercept)"] <- "Intercept"
    onemodel.glance <- onemodel.glance %>% mutate(term = "model") # add a column

    # get the list of terms:
    list.terms <- onemodel.tidy$term

    # check if the onemodel.* does not have real statistics (but only a column of 'term')
    temp_colnames <- onemodel.tidy %>% colnames()
    # union of colnames and "term"; if colnames only has "term" or lengt of 0 (tibble()),
    # union = "term", all(union)=TRUE;
    # otherwise, if there is colnames other than "term", all(union) = c(TRUE, FALSE, ...)
    temp <- union(temp_colnames, "term")
    # just an empty tibble (so below, all(dim(onemodel.tidy)) = FALSE)
    if (all(temp == "term")) onemodel.tidy <- tibble()

    temp_colnames <- onemodel.glance %>% colnames()
    temp <- union(temp_colnames, "term") # union of colnames and "term";
    if (all(temp == "term")) onemodel.glance <- tibble()



    # flatten .tidy results into one row:
    if (all(dim(onemodel.tidy))) { # not empty | if any dim is 0, all=FALSE
      onemodel.tidy.onerow <- onemodel.tidy %>% tidyr::pivot_wider(
        names_from = term,
        values_from = all_of(var.terms.orig),
        names_glue = "{term}.{.value}"
      )
    } else {
      onemodel.tidy.onerow <- onemodel.tidy
    }

    if (all(dim(onemodel.glance))) { # not empty
      onemodel.glance.onerow <- onemodel.glance %>% tidyr::pivot_wider(
        names_from = term,
        values_from = all_of(var.model),
        names_glue = "{term}.{.value}"
      )
    } else {
      onemodel.glance.onerow <- onemodel.glance
    }


    # combine the tables: check if any of them is empty (tibble())
    onemodel.onerow <- bind_cols_check_emptyTibble(onemodel.tidy.onerow, onemodel.glance.onerow)

    # add a column of element ids:
    colnames.temp <- colnames(onemodel.onerow)
    onemodel.onerow <- onemodel.onerow %>% tibble::add_column(
      element_id = i_element - 1,
      .before = colnames.temp[1]
    ) # add as the first column

    # now you can get the headers, # of columns, etc of the output results


    if (flag_initiate == TRUE) { # return the column names:

      # return:
      column_names <- colnames(onemodel.onerow)
      toreturn <- list(
        column_names = column_names,
        list.terms = list.terms
      )
      toreturn
    } else if (flag_initiate == FALSE) { # return the one row results:

      # return:
      onerow <- as.numeric(onemodel.onerow) # change from tibble to numeric to save some space
      onerow
    }
  } else { # if flag_sufficient==FALSE
    if (flag_initiate == TRUE) {
      toreturn <- list(
        column_names = NaN,
        list.terms = NaN
      )
      toreturn
    } else if (flag_initiate == FALSE) {
      onerow <- c(
        i_element - 1, # first is element_id, which will still be a finite number
        rep(NaN, (num.stat.output - 1))
      ) # other columns will be NaN

      onerow
    }
  }
}
#' Fit GAM for one element
#'
#' @description
#' `analyseOneElement.gam` fits a GAM model for one element data, and returns requested model statistics.
#'
#' @details
#' `ModelArray.gam` iteratively calls this function to get statistics for all requested elements.
#'
#' @param i_element An integer, the i_th element, starting from
#'     1. For initiating (flag_initiate = TRUE), use i_element=1
#' @param formula A formula (passed to `mgcv::gam()`)
#' @param modelarray ModelArray class
#' @param phenotypes A data.frame of the cohort with columns of independent variables and covariates
#'     to be added to the model
#' @param scalar A character. The name of the element-wise scalar to be analysed
#' @param var.smoothTerms The list of variables to save for smooth terms
#'     (got from broom::tidy(parametric = FALSE)). Example smooth term: age in formula "outcome ~ s(age)".
#' @param var.parametricTerms The list of variables to save for parametric terms (got from
#'     broom::tidy(parametric = TRUE)). Example parametric term: sex in formula "outcome ~ s(age) + sex".
#' @param var.model The list of variables to save for the model (got from broom::glance() and summary()).
#' @param num.subj.lthr The minimal number of subjects with valid value in input h5 file, i.e. number of
#'     subjects with finite values (defined by `is.finite()`, i.e. not NaN or NA or Inf) in h5 file >
#'     \code{num.subj.lthr}, then this element will be run normally; otherwise, this element will be skipped
#'     and statistical outputs will be set as NaN.
#' @param num.stat.output The number of output stat metrics (for generating all NaN stat when
#'     # subjects does not meet criteria). This includes column `element_id`.
#'     This is required when flag_initiate = TRUE.
#' @param flag_initiate TRUE or FALSE, Whether this is to initiate the new analysis.
#'     If TRUE, it will return column names etc to be used for initiating data.frame;
#'     if FALSE, it will return the list of requested statistic values.
#' @param flag_sse TRUE or FALSE, Whether to calculate SSE (sum of squared error) for the model
#'     (`model.sse`). SSE is needed for calculating partial R-squared.
#' @param ... Additional arguments for `mgcv::gam()`
#' @return If flag_initiate==TRUE, returns column names, list of term names of final results, and
#'     attr.name of sp.criterion; if flag_initiate==FALSE, it will return the list of requested
#'     statistic values for a element.
#' @export
#' @import mgcv
#' @import broom
#' @importFrom dplyr select %>% bind_cols
#' @import tibble

analyseOneElement.gam <- function(i_element, formula, modelarray, phenotypes, scalar,
                                  var.smoothTerms, var.parametricTerms, var.model,
                                  num.subj.lthr, num.stat.output = NULL,
                                  flag_initiate = FALSE, flag_sse = FALSE,
                                  ...) {
  values <- scalars(modelarray)[[scalar]][i_element, ]

  ## check number of subjects with (in)valid values:
  flag_sufficient <- NULL # whether number of subjects with valid values are sufficient
  num.subj.valid <- length(values[is.finite(values)])
  if (num.subj.valid > num.subj.lthr) {
    flag_sufficient <- TRUE
  } else {
    flag_sufficient <- FALSE
  }

  if (flag_sufficient == TRUE) {
    dat <- phenotypes
    dat[[scalar]] <- values

    arguments <- list(...)
    arguments$formula <- formula
    arguments$data <- dat

    # explicitly passing arguments into command, to avoid error of argument "weights"
    onemodel <- do.call(mgcv::gam, arguments)

    onemodel.tidy.smoothTerms <- onemodel %>% broom::tidy(parametric = FALSE)
    onemodel.tidy.parametricTerms <- onemodel %>% broom::tidy(parametric = TRUE)
    onemodel.glance <- onemodel %>% broom::glance()
    onemodel.summary <- onemodel %>% summary()
    # add additional model's stat to onemodel.glance():
    onemodel.glance[["adj.r.squared"]] <- onemodel.summary$r.sq
    onemodel.glance[["dev.expl"]] <- onemodel.summary$dev.expl

    sp.criterion.attr.name <- onemodel.summary$sp.criterion %>% attr(which = "name")
    onemodel.glance[["sp.criterion"]] <- onemodel.summary$sp.criterion[[sp.criterion.attr.name]]
    onemodel.glance[["scale"]] <- onemodel.summary$scale # scale estimate

    num.smoothTerms <- onemodel.summary$m # The number of smooth terms in the model.



    # delete columns you don't want:
    var.smoothTerms.full <- names(onemodel.tidy.smoothTerms)
    var.parametricTerms.full <- names(onemodel.tidy.parametricTerms)
    var.model.full <- names(onemodel.glance)

    # list to remove:
    var.smoothTerms.orig <- var.smoothTerms
    var.smoothTerms <- list("term", var.smoothTerms) %>% unlist() # we will always keep "term" column
    var.smoothTerms.remove <- list()
    for (l in var.smoothTerms.full) {
      if (!(l %in% var.smoothTerms)) {
        var.smoothTerms.remove <- var.smoothTerms.remove %>%
          append(., l) %>%
          unlist() # the order will still be kept
      }
    }

    var.parametricTerms.orig <- var.parametricTerms
    var.parametricTerms <- list("term", var.parametricTerms) %>% unlist() # we will always keep "term" column
    var.parametricTerms.remove <- list()
    for (l in var.parametricTerms.full) {
      if (!(l %in% var.parametricTerms)) {
        var.parametricTerms.remove <- var.parametricTerms.remove %>%
          append(., l) %>%
          unlist() # the order will still be kept
      }
    }

    var.model.remove <- list()
    for (l in var.model.full) {
      if (!(l %in% var.model)) {
        var.model.remove <- var.model.remove %>%
          append(., l) %>%
          unlist() # the order will still be kept
      }
    }

    # remove those columns:
    if (length(var.smoothTerms.remove) != 0) { # if length=0, it's list(), nothing to remove
      onemodel.tidy.smoothTerms <- dplyr::select(onemodel.tidy.smoothTerms, -all_of(var.smoothTerms.remove))
    }
    if (length(var.parametricTerms.remove) != 0) { # if length=0, it's list(), nothing to remove
      onemodel.tidy.parametricTerms <- dplyr::select(onemodel.tidy.parametricTerms, -all_of(var.parametricTerms.remove))
    }
    if (length(var.model.remove) != 0) {
      onemodel.glance <- dplyr::select(onemodel.glance, -all_of(var.model.remove))
    }

    # adjust:
    if (num.smoothTerms > 0) { # if there is any smooth term
      # change the term name from "(Intercept)" to "Intercept"
      onemodel.tidy.smoothTerms$term[onemodel.tidy.smoothTerms$term == "(Intercept)"] <- "Intercept"
    }
    if (nrow(onemodel.tidy.parametricTerms) > 0) { # if there is any parametric term
      # change the term name from "(Intercept)" to "Intercept"
      onemodel.tidy.parametricTerms$term[onemodel.tidy.parametricTerms$term == "(Intercept)"] <- "Intercept"
    }

    # change from s(x) to s_x: (could be s, te, etc); from s(x):oFactor to s_x_BYoFactor; from ti(x,z) to ti_x_z
    if (num.smoothTerms > 0) { # if there is any smooth term
      for (i_row in seq_along(onemodel.tidy.smoothTerms)) {
        # step 1: change from s(x) to s_x
        term_name <- onemodel.tidy.smoothTerms$term[i_row]
        str_list <- strsplit(term_name, split = "[()]")[[1]]

        str <- str_list[2] # extract string between ()
        smooth_name <- str_list[1] # "s" or some other smooth method type such as "te"
        str_valid <- paste0(smooth_name, "_", str)

        if (length(str_list) > 2) { # there is string after variable name
          str_valid <- paste0(
            str_valid, "_",
            paste(str_list[3:length(str_list)], collapse = "")
          ) # combine rest of strings
        }

        # detect ":", and change to "BY"   # there is "_" replacing for ")" in "s()" already
        str_valid <- gsub(":", "BY", str_valid, fixed = TRUE)

        # detect ",", and change to "_"
        str_valid <- gsub(",", "_", str_valid, fixed = TRUE)

        onemodel.tidy.smoothTerms$term[i_row] <- str_valid
      }
    }


    onemodel.glance <- onemodel.glance %>% mutate(term = "model") # add a column

    # get the list of terms:
    if (num.smoothTerms > 0) {
      list.smoothTerms <- onemodel.tidy.smoothTerms$term # if empty, gives warning
    } else {
      list.smoothTerms <- NULL
    }

    if (nrow(onemodel.tidy.parametricTerms) > 0) {
      list.parametricTerms <- onemodel.tidy.parametricTerms$term
    } else {
      list.parametricTerms <- NULL
    }


    # check if the onemodel.* does not have real statistics (but only a column of 'term')
    temp_colnames <- onemodel.tidy.smoothTerms %>% colnames()
    # union of colnames and "term"; if colnames only has "term" or lengt of 0 (tibble()),
    # union = "term", all(union)=TRUE;
    # otherwise, if there is colnames other than "term", all(union) = c(TRUE, FALSE, ...)
    temp <- union(temp_colnames, "term")
    # just an empty tibble (so below, all(dim(onemodel.tidy.smoothTerms)) = FALSE)
    if (all(temp == "term")) onemodel.tidy.smoothTerms <- tibble()

    temp_colnames <- onemodel.tidy.parametricTerms %>% colnames()
    temp <- union(temp_colnames, "term")
    if (all(temp == "term")) onemodel.tidy.parametricTerms <- tibble() # just an empty tibble

    temp_colnames <- onemodel.glance %>% colnames()
    temp <- union(temp_colnames, "term")
    if (all(temp == "term")) onemodel.glance <- tibble() # just an empty tibble


    # flatten .tidy results into one row:
    if (all(dim(onemodel.tidy.smoothTerms))) { # not empty | if any dim is 0, all=FALSE
      onemod_smoothTerms_1row <- onemodel.tidy.smoothTerms %>% tidyr::pivot_wider(
        names_from = term,
        values_from = all_of(var.smoothTerms.orig),
        names_glue = "{term}.{.value}"
      )
    } else {
      onemod_smoothTerms_1row <- onemodel.tidy.smoothTerms
    }

    if (all(dim(onemodel.tidy.parametricTerms))) { # not empty
      onemod_paramTerms_1row <- onemodel.tidy.parametricTerms %>% tidyr::pivot_wider(
        names_from = term,
        values_from = all_of(var.parametricTerms.orig),
        names_glue = "{term}.{.value}"
      )
    } else {
      onemod_paramTerms_1row <- onemodel.tidy.parametricTerms
    }

    if (all(dim(onemodel.glance))) { # not empty
      onemodel.glance.onerow <- onemodel.glance %>% tidyr::pivot_wider(
        names_from = term,
        values_from = all_of(var.model),
        names_glue = "{term}.{.value}"
      )
    } else {
      onemodel.glance.onerow <- onemodel.glance
    }


    # combine the tables: check if any of them is empty (tibble())
    onemodel.onerow <- bind_cols_check_emptyTibble(
      onemod_smoothTerms_1row,
      onemod_paramTerms_1row
    )
    onemodel.onerow <- bind_cols_check_emptyTibble(
      onemodel.onerow,
      onemodel.glance.onerow
    )


    # add a column of element ids:
    colnames.temp <- colnames(onemodel.onerow)
    onemodel.onerow <- onemodel.onerow %>% tibble::add_column(
      element_id = i_element - 1,
      .before = colnames.temp[1]
    ) # add as the first column

    # add sse if requested:
    if (flag_sse == TRUE) {
      # using values from model itself, where NAs in y have been excluded
      # --> sse won't be NA --> partial R-squared won't be NA
      onemodel.onerow[["model.sse"]] <- sum((onemodel$y - onemodel$fitted.values)^2)
    }


    # now you can get the headers, # of columns, etc of the output results


    if (flag_initiate == TRUE) { # return the column names:

      # return:
      column_names <- colnames(onemodel.onerow)
      toreturn <- list(
        column_names = column_names,
        list.smoothTerms = list.smoothTerms,
        list.parametricTerms = list.parametricTerms,
        sp.criterion.attr.name = sp.criterion.attr.name
      )
      toreturn
    } else if (flag_initiate == FALSE) { # return the one row results:

      # return:
      onerow <- as.numeric(onemodel.onerow) # change from tibble to numeric to save some space
      onerow
    }
  } else { # if flag_sufficient==FALSE
    if (flag_initiate == TRUE) {
      toreturn <- list(
        column_names = NaN,
        list.smoothTerms = NaN,
        list.parametricTerms = NaN,
        sp.criterion.attr.name = NaN
      )
      toreturn
    } else if (flag_initiate == FALSE) {
      onerow <- c(
        i_element - 1, # first is element_id, which will still be a finite number
        rep(NaN, (num.stat.output - 1))
      ) # other columns will be NaN

      onerow
    }
  }
}




#' Write outputs from element-wise statistical analysis to the HDF5 file.
#'
#' @description
#' Create a group named `analysis_name` in HDF5 file, then write the statistical results
#' data.frame (i.e. for one analysis) in it.
#'
#' @details
#' debug tip: For "Error in H5File.open(filename, mode, file_create_pl, file_access_pl)",
#' check if there is message 'No such file or directory'. Try absolute .h5 filename.
#'
#' @param fn.output A character, The HDF5 (.h5) filename for the output
#' @param df.output A data.frame object with element-wise statistical results, returned from `ModelArray.lm()` etc
#' @param analysis_name A character, the name of the results
#' @param overwrite If a group with the same analysis_name exists in HDF5 file, whether overwrite it (TRUE)
#'     or not (FALSE)
#' @import hdf5r
#' @export
writeResults <- function(fn.output, df.output, analysis_name = "myAnalysis", overwrite = TRUE) {
  # This is enhanced version with: 1) change to hdf5r; 2) write results with only one row for one element

  # check "df.output"
  if (!("data.frame" %in% class(df.output))) {
    stop("Results dataset is not correct; must be data of type `data.frame`")
  }

  fn.output.h5 <- hdf5r::H5File$new(fn.output, mode = "a")
  # open; "a": creates a new file or opens an existing one for read/write

  # check if group "results" already exists!
  if (fn.output.h5$exists("results") == TRUE) { # group "results" exist
    results.grp <- fn.output.h5$open("results")
  } else {
    results.grp <- fn.output.h5$create_group("results")
  }

  # check if group "results\<analysis_name>" exists:
  if (results.grp$exists(analysis_name) == TRUE & overwrite == FALSE) {
    warning(paste0(analysis_name, " exists but not to overwrite!"))
    # TODO: add checker for exisiting analysis_name, esp the matrix size
    results.analysis.grp <- results.grp$open(analysis_name)
    results_matrix_ds <- results.analysis.grp[["results_matrix"]]
  } else { # not exist; or exist & overwrite: to create
    if (results.grp$exists(analysis_name) == TRUE & overwrite == TRUE) { # delete existing one first
      results.grp$link_delete(analysis_name)
      # NOTE: the file size will not shrink after your deletion.. this is because of HDF5,
      # regardless of package of hdf5r or rhdf5
      # TODO: add a garbage collector after saving the results
    }

    # create:
    results.analysis.grp <- results.grp$create_group(analysis_name)
    # create a subgroup called analysis_name under results.grp

    # check "df.output": make sure all columns are floats (i.e. numeric)
    for (i_col in seq(1, ncol(df.output), by = 1)) { # for each column of df.output
      col_class <- as.character(sapply(df.output, class)[i_col]) # class of this column

      if ((col_class != "numeric") & (col_class != "integer")) { # the column class is not numeric or integer
        message(
          paste0(
            "the column #",
            as.character(i_col),
            " of df.output to save: data class is not numeric or integer...fixing it"
          )
        )

        # turn into numeric & write the notes in .h5 file...:
        factors <- df.output %>%
          pull(., var = i_col) %>%
          factor()
        df.output[, i_col] <- df.output %>%
          pull(., var = i_col) %>%
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

    # attach column names:
    hdf5r::h5attr(results.analysis.grp[["results_matrix"]], "colnames") <- colnames(df.output)
    # NOTES: update ConFixel correspondingly
  }

  # # return:   # will not work if fn.output.h5$close_all()
  # output_list <- list(results.grp = results.grp,
  #                     results.analysis.grp = results.analysis.grp,
  #                     results_matrix_ds = results_matrix_ds)
  # return(output_list)

  fn.output.h5$close_all()


  # message("Results file written!")
}

#' @export
#'
generator_gamFormula_factorXsmooth <- function(response.var, factor.var, smooth.var, phenotypes,
                                         reference.group = NULL, prefix.ordered.factor = "o",
                                         fx = TRUE, k = NULL) {
  # ... existing code ...
}

#' @export
#'
generator_gamFormula_continuousInteraction <- function(response.var, cont1.var, cont2.var,
                                  fx = TRUE, k = NULL) {
  # ... existing code ...
}
