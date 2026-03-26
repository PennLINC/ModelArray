# Internal shared helpers for ModelArray analysis functions
# These are not exported — used by ModelArray.lm, ModelArray.gam, ModelArray.wrap


#' Validate that data is a ModelArray
#' @noRd
.validate_modelarray_input <- function(data) {
  if (!inherits(data, "ModelArray")) {
    stop("data's class is not ModelArray!")
  }
}


#' Validate and default element.subset
#' @return Integer vector of element indices
#' @noRd
.validate_element_subset <- function(element.subset, data, scalar) {
  if (is.null(element.subset)) {
    num.element.total <- numElementsTotal(modelarray = data, scalar_name = scalar)
    element.subset <- 1:num.element.total
  }
  if (min(element.subset) < 1) {
    stop("Minimal value in element.subset should >= 1")
  }
  if (max(element.subset) > nrow(scalars(data)[[scalar]])) {
    stop(
      paste0(
        "Maximal value in element.subset should <= number of elements = ",
        as.character(nrow(scalars(data)[[scalar]]))
      )
    )
  }
  if (!is.integer(element.subset)) {
    stop("Please enter integers for element.subset!")
  }
  element.subset
}


#' Align phenotypes to ModelArray sources
#'
#' Checks that source_file column exists, lengths match, entries are unique,
#' and reorders phenotypes to match ModelArray source order if needed.
#' @return Reordered phenotypes data.frame
#' @noRd
.align_phenotypes <- function(data, phenotypes, scalar) {
  sources.modelarray <- sources(data)[[scalar]]
  sources.phenotypes <- phenotypes[["source_file"]]
  if (is.null(sources.phenotypes)) {
    stop("Did not find column 'source_file' in argument 'phenotypes'. Please check!")
  }

  if (length(sources.modelarray) != length(sources.phenotypes)) {
    stop(
      paste0(
        "The length of the source file list in phenotypes's column 'source_file' ",
        "is not the same as the length of the source file list in ModelArray 'data'! Please check out! ",
        "The latter can be accessed by: sources(data)[[scalar]]"
      )
    )
  }

  if (length(sources.modelarray) != length(unique(sources.modelarray))) {
    stop(
      paste0(
        "The source files in ModelArray 'data' are not unique! Please investigate! ",
        "It can be accessed by: sources(data)[[scalar]]"
      )
    )
  }
  if (length(sources.phenotypes) != length(unique(sources.phenotypes))) {
    stop(
      paste0(
        "The source files from phenotypes's column 'source_file' ",
        "are not unique! Please investigate and remove the duplicated one!"
      )
    )
  }

  if (!identical(sources.modelarray, sources.phenotypes)) {
    both_match <- (all(sources.modelarray %in% sources.phenotypes)) &&
      (all(sources.phenotypes %in% sources.modelarray))
    if (both_match) {
      reorder_idx <- match(sources.modelarray, sources.phenotypes)
      phenotypes <- phenotypes[reorder_idx, ]
      row.names(phenotypes) <- NULL
      if (!identical(phenotypes[["source_file"]], sources.modelarray)) {
        stop("matching source file names were not successful...")
      }
    } else {
      stop(
        paste0(
          "phenotypes's column 'source_file' has different element(s) than the source file list",
          " in ModelArray 'data'! Please investigate! ",
          "The latter can be accessed by: sources(data)[[scalar]]"
        )
      )
    }
  }

  phenotypes
}


#' Compute subject threshold from absolute and relative parameters
#' @return Numeric threshold value
#' @noRd
.compute_subject_threshold <- function(phenotypes, num.subj.lthr.abs, num.subj.lthr.rel) {
  num.subj.total <- nrow(phenotypes)
  max(num.subj.total * num.subj.lthr.rel, num.subj.lthr.abs)
}


#' Find a valid initiator element by searching middle → forward → backward
#'
#' @param analyse_one_fn The per-element analysis function
#' @param num.elements.total Total number of elements
#' @param extra_args Named list of extra arguments for analyse_one_fn
#'   (must NOT include i_element, num.stat.output, or flag_initiate)
#' @param verbose Whether to print progress messages
#' @return The outputs_initiator list from the first successful element
#' @noRd
.find_initiator_element <- function(analyse_one_fn, num.elements.total,
                                    extra_args, verbose = TRUE) {
  i_element_try <- floor(num.elements.total / 2)

  call_init <- function(i_element) {
    do.call(analyse_one_fn, c(
      list(i_element = i_element),
      extra_args,
      list(num.stat.output = NULL, flag_initiate = TRUE)
    ))
  }

  outputs_initiator <- call_init(i_element_try)
  if (!is.nan(outputs_initiator$column_names[1])) {
    return(list(outputs = outputs_initiator, i_element = i_element_try))
  }

  # Try forward from middle to end
  if (verbose) {
    message(
      "There are insufficient valid subjects for initiating with the middle element; ",
      "trying other elements; this may take a while...."
    )
  }
  for (i_element_temp in (i_element_try + 1):num.elements.total) {
    if (i_element_temp %% 100 == 0) {
      message(
        "Trying element #", toString(i_element_temp),
        " and the following elements for initiating...."
      )
    }
    outputs_initiator <- call_init(i_element_temp)
    if (!is.nan(outputs_initiator$column_names[1])) {
      return(list(outputs = outputs_initiator, i_element = i_element_temp))
    }
  }

  # Try backward from start to middle
  message(
    "there no elements with sufficient valid ",
    "subjects for initiating the process... ",
    "trying element #1 and the following elements for initiating; ",
    "this may take a while...."
  )
  for (i_element_temp in 1:(i_element_try - 1)) {
    if (i_element_temp %% 100 == 0) {
      message(
        "trying element #", toString(i_element_temp),
        " and the following elements for initiating...."
      )
    }
    outputs_initiator <- call_init(i_element_temp)
    if (!is.nan(outputs_initiator$column_names[1])) {
      return(list(outputs = outputs_initiator, i_element = i_element_temp))
    }
  }

  stop(
    "Have tried all elements, but there are no elements with sufficient valid, ",
    "finite h5 scalar values (i.e. not NaN or NA, not infinite). ",
    "Please check if thresholds 'num.subj.lthr.abs' and 'num.subj.lthr.rel' were set too high, ",
    "or there are problems in the group mask or individual masks!"
  )
}


#' Dispatch parallel/serial apply with optional progress bar
#'
#' @param element.subset Integer vector of element indices
#' @param FUN The per-element function (first arg must be i_element)
#' @param n_cores Number of cores
#' @param pbar Whether to show progress bar
#' @param ... Additional arguments passed to FUN
#' @return List of results from FUN
#' @noRd
.parallel_dispatch <- function(element.subset, FUN, n_cores, pbar, ...) {
  if (pbar) {
    old_pb_opts <- pbapply::pboptions()
    pbapply::pboptions(type = "txt")
    on.exit(pbapply::pboptions(old_pb_opts), add = TRUE)
  }

  if (n_cores > 1) {
    if (pbar) {
      fits <- pbmcapply::pbmclapply(element.subset,
        FUN,
        mc.cores = n_cores,
        ignore.interactive = TRUE,
        ...
      )
    } else {
      fits <- parallel::mclapply(element.subset,
        FUN,
        mc.cores = n_cores,
        ...
      )
    }
  } else {
    if (pbar) {
      fits <- pbapply::pblapply(element.subset,
        FUN,
        ...
      )
    } else {
      fits <- lapply(element.subset,
        FUN,
        ...
      )
    }
  }

  fits
}


#' Correct p-values for a set of terms and append corrected columns
#'
#' @param df_out Data.frame of results
#' @param term_list Character vector of term names (e.g., c("Age", "Sex") or c("model"))
#' @param correct_methods Character vector of correction methods (e.g., c("fdr"))
#' @param var_list Character vector of requested variables (checked for "p.value")
#' @return Modified df_out with corrected p-value columns inserted
#' @importFrom dplyr %>%
#' @noRd
.correct_pvalues <- function(df_out, term_list, correct_methods, var_list) {
  if (all(correct_methods == "none")) {
    return(df_out)
  }
  if (!("p.value" %in% var_list)) {
    return(df_out)
  }

  for (methodstr in correct_methods) {
    for (tempstr in term_list) {
      tempstr.raw <- paste0(tempstr, ".p.value")
      tempstr.corrected <- paste0(tempstr.raw, ".", methodstr)
      temp.corrected <- stats::p.adjust(df_out[[tempstr.raw]], method = methodstr)
      df_out <- df_out %>%
        tibble::add_column("{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
    }
  }

  df_out
}
