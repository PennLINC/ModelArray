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
      fits <- pbapply::pblapply(
        element.subset,
        FUN,
        ...
      )
    } else {
      fits <- lapply(
        element.subset,
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


#' Initialize an incremental writer for /results datasets
#' @noRd
.init_results_stream_writer <- function(write_results_name,
                                        write_results_file,
                                        n_rows,
                                        column_names,
                                        flush_every = 1000L,
                                        storage_mode = "double",
                                        compression_level = 4L) {
  if (is.null(write_results_name)) {
    return(NULL)
  }
  if (!is.character(write_results_name) || length(write_results_name) != 1L || write_results_name == "") {
    stop("write_results_name must be a non-empty character string when provided")
  }
  if (is.null(write_results_file) || !is.character(write_results_file) || length(write_results_file) != 1L) {
    stop("write_results_file must be a single character path when write_results_name is provided")
  }
  if (!is.numeric(flush_every) || length(flush_every) != 1L || flush_every <= 0) {
    stop("write_results_flush_every must be a positive integer")
  }
  flush_every <- as.integer(flush_every)

  if (!file.exists(write_results_file)) {
    rhdf5::h5createFile(write_results_file)
  }

  h5_write <- hdf5r::H5File$new(write_results_file, mode = "a")
  if (!h5_write$exists("results")) {
    h5_write$create_group("results")
  }
  results_grp <- h5_write$open("results")
  if (results_grp$exists(write_results_name)) {
    results_grp$link_delete(write_results_name)
  }
  results_grp$create_group(write_results_name)
  h5_write$close_all()

  dataset_path <- paste0("results/", write_results_name, "/results_matrix")
  chunk_rows <- min(flush_every, n_rows)

  rhdf5::h5createDataset(
    file = write_results_file,
    dataset = dataset_path,
    dims = c(n_rows, length(column_names)),
    storage.mode = storage_mode,
    chunk = c(chunk_rows, length(column_names)),
    level = as.integer(compression_level)
  )

  list(
    file = write_results_file,
    dataset_path = dataset_path,
    names_path = paste0("results/", write_results_name, "/column_names"),
    column_names = as.character(column_names),
    n_cols = length(column_names),
    write_row_cursor = 1L
  )
}


#' Append one block to an incremental /results writer
#' @noRd
.results_stream_write_block <- function(writer, block_df) {
  if (is.null(writer)) {
    return(writer)
  }
  block <- as.matrix(block_df)
  row_idx <- writer$write_row_cursor:(writer$write_row_cursor + nrow(block) - 1L)
  rhdf5::h5write(
    obj = block,
    file = writer$file,
    name = writer$dataset_path,
    index = list(row_idx, seq_len(writer$n_cols))
  )
  writer$write_row_cursor <- writer$write_row_cursor + nrow(block)
  writer
}


#' Finalize an incremental /results writer
#' @noRd
.finalize_results_stream_writer <- function(writer) {
  if (is.null(writer)) {
    return(invisible(NULL))
  }
  rhdf5::h5write(
    obj = writer$column_names,
    file = writer$file,
    name = writer$names_path
  )
  rhdf5::h5closeAll()
  invisible(NULL)
}
