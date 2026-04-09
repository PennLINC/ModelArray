# Internal shared helpers for ModelArray analysis functions
# These are not exported — used by ModelArray.lm, ModelArray.gam, ModelArray.wrap

# Validation and dataset preparation ----
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

# Create analysis context ----
# Precomputed context builders for element-wise analysis functions.
# These hoist loop-invariant computations out of the per-element functions
# so they execute once before the loop rather than millions of times inside it.

# Level 1: Base context shared by all analysis functions (lm, gam, wrap).
# Handles scalar discovery, collision checks, and cross-scalar source alignment.

#' Build base context shared by all element-wise analysis functions
#'
#' Precomputes scalar name lookups, column-name collision checks, and
#' cross-scalar source-alignment indices that are invariant across elements.
#' This is the foundation on which model-specific contexts are built.
#'
#' @param modelarray A \linkS4class{ModelArray} object.
#' @param phenotypes A data.frame, already aligned to \code{modelarray} via
#'   \code{.align_phenotypes()}.
#' @param scalar Character. The name of the response scalar.
#' @param scalar_subset Character vector or \code{NULL}. Which scalars to
#'   attach and align. If \code{NULL}, all scalars in the ModelArray are used
#'   (the \code{wrap} behaviour). If a character vector, only those scalars
#'   are included (the \code{lm}/\code{gam} behaviour where we detect
#'   formula-referenced scalars).
#' @return A named list with components:
#'   \describe{
#'     \item{modelarray}{The ModelArray object (read-only reference).}
#'     \item{phenotypes}{The aligned phenotypes data.frame.}
#'     \item{scalar}{Character. The response scalar name.}
#'     \item{all_scalar_names}{Character vector. Names of all scalars in the
#'       ModelArray.}
#'     \item{attached_scalars}{Character vector. Names of scalars that will
#'       be attached to the per-element data.frame (may be a subset of
#'       \code{all_scalar_names} or all of them).}
#'     \item{predictor_reorder}{A named list of integer vectors. For each
#'       scalar in \code{attached_scalars}, the precomputed index vector
#'       from \code{match(phen_sources, scalar_sources)} that reorders
#'       scalar columns to match phenotype rows. For the response scalar
#'       (which is already aligned by \code{.align_phenotypes()}), this
#'       entry is \code{NULL} (no reordering needed).}
#'   }
#' @noRd
.build_base_context <- function(modelarray, phenotypes, scalar,
                                scalar_subset = NULL) {
  all_scalar_names <- names(scalars(modelarray))


  # Determine which scalars to attach

  if (is.null(scalar_subset)) {
    # wrap mode: attach all scalars
    attached_scalars <- all_scalar_names
  } else {
    # lm/gam mode: attach only the explicitly requested scalars
    # (response + any predictor scalars detected from the formula)
    attached_scalars <- unique(scalar_subset)
  }


  # Collision check: scalar names vs phenotype column names.
  # This is invariant — the same collision would occur at every element,

  # so we fail fast once rather than millions of times.

  collisions <- intersect(attached_scalars, colnames(phenotypes))
  if (length(collisions) > 0) {
    stop(
      "Column name collision between phenotypes and scalar names: ",
      paste(collisions, collapse = ", "),
      ". Please rename or remove these columns from phenotypes before ",
      "modeling."
    )
  }

  # Precompute source-alignment reorder indices for each attached scalar.
  # The response scalar is already aligned by .align_phenotypes(), so it

  # needs no reordering. Predictor scalars from mergeModelArrays() may
  # have different column orderings and need match()-based reordering.
  phen_sources <- phenotypes[["source_file"]]
  predictor_reorder <- stats::setNames(
    vector("list", length(attached_scalars)),
    attached_scalars
  )

  for (sname in attached_scalars) {
    s_sources <- sources(modelarray)[[sname]]

    if (identical(s_sources, phen_sources)) {
      # Already aligned (typical for the response scalar, or when all
      # scalars share the same source order after mergeModelArrays)
      predictor_reorder[[sname]] <- NULL
    } else {
      # Validate that the sets match (once, not per element)
      if (!(all(s_sources %in% phen_sources) &&
            all(phen_sources %in% s_sources))) {
        stop(
          "The source files for scalar '", sname,
          "' do not match phenotypes$source_file."
        )
      }
      predictor_reorder[[sname]] <- match(phen_sources, s_sources)
    }
  }

  list(
    modelarray        = modelarray,
    phenotypes        = phenotypes,
    scalar            = scalar,
    all_scalar_names  = all_scalar_names,
    attached_scalars  = attached_scalars,
    predictor_reorder = predictor_reorder
  )
}


# Level 2: Model-specific context builders that add formula-derived
# information on top of the base context.

#' Build context for element-wise linear model fitting
#'
#' Extends the base context with formula parsing results: the LHS variable
#' name, RHS variable names, and which RHS variables are scalars (as opposed
#' to phenotype columns). All of this is invariant across elements.
#'
#' @param formula A \code{\link[stats]{formula}} to be passed to
#'   \code{\link[stats]{lm}}.
#' @param modelarray A \linkS4class{ModelArray} object.
#' @param phenotypes A data.frame, already aligned via
#'   \code{.align_phenotypes()}.
#' @param scalar Character. The response scalar name.
#' @return A named list inheriting all components from
#'   \code{.build_base_context()} plus:
#'   \describe{
#'     \item{formula}{The model formula.}
#'     \item{lhs_name}{Character. The response variable name from the
#'       formula LHS.}
#'     \item{rhs_vars}{Character vector. All variable names on the RHS.}
#'     \item{scalar_predictors}{Character vector. RHS variables that are
#'       scalar names in the ModelArray (i.e., cross-scalar predictors).}
#'   }
#' @noRd
.build_lm_context <- function(formula, modelarray, phenotypes, scalar) {
  # Formula parsing — currently repeated per element inside
  # analyseOneElement.lm [4]
  all_vars <- all.vars(formula)
  lhs_name <- tryCatch(
    as.character(formula[[2]]),
    error = function(e) NULL
  )
  rhs_vars <- setdiff(all_vars, lhs_name)

  # Detect which RHS variables are scalars in the ModelArray
  all_scalar_names <- names(scalars(modelarray))
  scalar_predictors <- intersect(rhs_vars, all_scalar_names)

  # The scalars to attach: response + any scalar predictors
  scalar_subset <- unique(c(scalar, scalar_predictors))

  # Build the base context with only the needed scalars
  base_ctx <- .build_base_context(
    modelarray    = modelarray,
    phenotypes    = phenotypes,
    scalar        = scalar,
    scalar_subset = scalar_subset
  )

  # Extend with formula-specific fields
  base_ctx$formula            <- formula
  base_ctx$lhs_name           <- lhs_name
  base_ctx$rhs_vars           <- rhs_vars
  base_ctx$scalar_predictors  <- scalar_predictors

  base_ctx
}


#' Build context for element-wise GAM fitting
#'
#' Extends the base context with formula parsing results and GAM-specific
#' formula validation that currently runs inside the ModelArray.gam()
#' preamble but whose results are never passed to the per-element function.
#'
#' @param formula A \code{\link[stats]{formula}} to be passed to
#'   \code{\link[mgcv]{gam}}.
#' @param modelarray A \linkS4class{ModelArray} object.
#' @param phenotypes A data.frame, already aligned via
#'   \code{.align_phenotypes()}.
#' @param scalar Character. The response scalar name.
#' @return A named list inheriting all components from
#'   \code{.build_lm_context()} (which itself inherits from
#'   \code{.build_base_context()}) plus:
#'   \describe{
#'     \item{gam_formula_breakdown}{The result of
#'       \code{mgcv::interpret.gam(formula)}, cached for reuse.}
#'   }
#' @noRd
.build_gam_context <- function(formula, modelarray, phenotypes, scalar) {
  # GAM formula validation — currently runs in ModelArray.gam() preamble [1]
  # but the breakdown result is discarded. We cache it here.
  gam_breakdown <- tryCatch(
    mgcv::interpret.gam(formula),
    error = function(cond) {
      stop("The formula is not valid for mgcv::gam()! Please check and revise.")
    }
  )

  # The formula structure is the same as lm for variable detection purposes
  ctx <- .build_lm_context(formula, modelarray, phenotypes, scalar)

  # Add GAM-specific cached data
  ctx$gam_formula_breakdown <- gam_breakdown

  ctx
}


#' Build context for element-wise user-supplied function execution
#'
#' Uses the base context in "attach all scalars" mode since
#' \code{analyseOneElement.wrap} attaches every scalar in the ModelArray
#' to the per-element data.frame, not just formula-referenced ones.
#'
#' @param modelarray A \linkS4class{ModelArray} object.
#' @param phenotypes A data.frame, already aligned via
#'   \code{.align_phenotypes()}.
#' @param scalar Character. The primary scalar name (used for the initial
#'   validity check).
#' @return A named list: the base context with \code{scalar_subset = NULL}
#'   (all scalars attached).
#' @noRd
.build_wrap_context <- function(modelarray, phenotypes, scalar) {
  # scalar_subset = NULL triggers "attach all" mode in .build_base_context()
  ctx <- .build_base_context(
    modelarray    = modelarray,
    phenotypes    = phenotypes,
    scalar        = scalar,
    scalar_subset = NULL
  )

  ctx
}


# Level 3: Shared per-element data assembly helper.
# Extracts scalar values for one element using the precomputed context,
# builds the validity mask, and returns the filtered data.frame.
# This replaces duplicated logic across analyseOneElement.lm,
# analyseOneElement.gam, and analyseOneElement.wrap [4].

#' Assemble per-element data.frame from precomputed context
#'
#' Reads scalar rows from the ModelArray, applies precomputed reorder
#' indices, builds the intersection validity mask, and returns the
#' filtered data.frame ready for model fitting or user function execution.
#'
#' @param i_element Integer. 1-based element index.
#' @param ctx A context list from one of the \code{.build_*_context()}
#'   functions.
#' @param num.subj.lthr Numeric. Minimum number of subjects with finite
#'   values required.
#' @return A list with components:
#'   \describe{
#'     \item{dat}{A data.frame: the filtered phenotypes with scalar columns
#'       attached. \code{NULL} if the element has insufficient valid
#'       subjects.}
#'     \item{sufficient}{Logical. Whether the element passed the subject
#'       threshold.}
#'     \item{num_valid}{Integer. Number of subjects with finite values
#'       across all attached scalars.}
#'   }
#' @noRd
.assemble_element_data <- function(i_element, ctx, num.subj.lthr) {
  # Read the response scalar row — the only mandatory per-element I/O
  response_vals <- scalars(ctx$modelarray)[[ctx$scalar]][i_element, ]

  # Start the validity mask with the response scalar
  masks <- list(is.finite(response_vals))

  # Read and reorder additional attached scalars
  scalar_values <- list()
  scalar_values[[ctx$scalar]] <- response_vals

  other_scalars <- setdiff(ctx$attached_scalars, ctx$scalar)
  for (sname in other_scalars) {
    s_vals <- scalars(ctx$modelarray)[[sname]][i_element, ]

    # Apply precomputed reorder index (or use as-is if NULL)
    reorder_idx <- ctx$predictor_reorder[[sname]]
    if (!is.null(reorder_idx)) {
      s_vals <- s_vals[reorder_idx]
    }

    scalar_values[[sname]] <- s_vals
    masks[[length(masks) + 1L]] <- is.finite(s_vals)
  }

  # Intersection mask across all scalars
  valid_mask <- Reduce("&", masks)
  num_valid <- sum(valid_mask)

  if (!(num_valid > num.subj.lthr)) {
    return(list(dat = NULL, sufficient = FALSE, num_valid = num_valid))
  }

  # Build filtered data.frame
  dat <- ctx$phenotypes[valid_mask, , drop = FALSE]
  for (sname in ctx$attached_scalars) {
    dat[[sname]] <- scalar_values[[sname]][valid_mask]
  }

  list(dat = dat, sufficient = TRUE, num_valid = num_valid)
}

# Find initiator element ----
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


# Parallelization ----
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

# P-value correction ----
#' Correct p-values for a set of terms and append corrected columns
#'
#' @param df_out Data.frame of results
#' @param term_list Character vector of term names (e.g., c("Age", "Sex") or c("model"))
#' @param correct_methods Character vector of correction methods (e.g., c("fdr"))
#' @param var_list Character vector of requested variables (checked for "p.value")
#' @return Modified df_out with corrected p-value columns inserted
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
      df_out[[tempstr.corrected]] <- stats::p.adjust(df_out[[tempstr.raw]], method = methodstr)
    }
  }

  df_out
}

# Streaming writes ----
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
