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


#' Resolve the response scalar for formula-based analyses
#'
#' Uses the left-hand side of the formula when `scalar` is omitted.
#' If `scalar` is supplied, validates that it exists and matches the formula
#' response when that response can be mapped to a scalar name.
#' @noRd
.resolve_formula_scalar <- function(formula, data, scalar = NULL) {
  scalar_names <- names(scalars(data))
  lhs_expr <- tryCatch(formula[[2]], error = function(e) NULL)
  lhs_vars <- if (is.null(lhs_expr)) character(0) else all.vars(lhs_expr)
  lhs_scalar <- intersect(lhs_vars, scalar_names)

  if (is.null(scalar)) {
    if (length(lhs_scalar) == 1L) {
      return(lhs_scalar)
    }
    stop(
      "Could not infer the response scalar from the formula. ",
      "Please supply `scalar=` explicitly."
    )
  }

  if (!(scalar %in% scalar_names)) {
    stop(
      "scalar '", scalar, "' not found in data. Available scalars: ",
      paste(scalar_names, collapse = ", ")
    )
  }

  if (length(lhs_scalar) == 1L && !identical(lhs_scalar, scalar)) {
    stop(
      "The formula response refers to scalar '", lhs_scalar,
      "' but `scalar=` was set to '", scalar, "'."
    )
  }

  scalar
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
  # analyseOneElement.lm
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
  # GAM formula validation — currently runs in ModelArray.gam() preamble
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
# analyseOneElement.gam, and analyseOneElement.wrap.

#' Return the storage path associated with one scalar
#' @noRd
.scalar_storage_path <- function(modelarray, scalar) {
  paths <- modelarray@path
  if (length(paths) == 0L) {
    return(NA_character_)
  }
  if (!is.null(names(paths)) && scalar %in% names(paths)) {
    return(paths[[scalar]])
  }
  paths[[1]]
}


#' Return TRUE when any attached scalar is TileDB-backed
#' @noRd
.has_tiledb_attached_scalar <- function(ctx) {
  any(vapply(ctx$attached_scalars, function(sname) {
    path <- .scalar_storage_path(ctx$modelarray, sname)
    if (is.na(path) || is.null(path)) {
      return(FALSE)
    }
    identical(.modelarray_backend(ctx$modelarray, path), "tiledb")
  }, logical(1)))
}


#' Resolve the scalar row cache block size for a context
#'
#' TileDB row-at-a-time reads are expensive in element-wise analyses. When a
#' TileDB-backed scalar participates in the model, materialize rows in chunks
#' before dispatching to workers. Set option `ModelArray.tiledb_read_block_size`
#' to override the memory-bounded default, or set it to 0 to disable. Set
#' `ModelArray.tiledb_read_block_mb` to tune the default memory target.
#' @noRd
.scalar_read_block_size <- function(ctx) {
  if (!.has_tiledb_attached_scalar(ctx)) {
    return(NULL)
  }

  block_size <- getOption("ModelArray.tiledb_read_block_size", NULL)
  if (is.null(block_size)) {
    block_mb <- getOption("ModelArray.tiledb_read_block_mb", 512)
    if (!is.numeric(block_mb) || length(block_mb) != 1L || is.na(block_mb) || block_mb <= 0) {
      stop("Option ModelArray.tiledb_read_block_mb must be a positive single number")
    }
    n_cols <- max(vapply(ctx$attached_scalars, function(sname) {
      ncol(scalars(ctx$modelarray)[[sname]])
    }, numeric(1)), 1L)
    bytes_per_row <- max(1L, n_cols) * max(1L, length(ctx$attached_scalars)) * 8L
    target_bytes <- block_mb * 1024^2
    block_size <- max(1L, min(.Machine$integer.max, floor(target_bytes / bytes_per_row)))
  }

  if (!is.numeric(block_size) || length(block_size) != 1L || is.na(block_size)) {
    stop("Option ModelArray.tiledb_read_block_size must be a single number")
  }
  block_size <- as.integer(min(block_size, .Machine$integer.max))
  if (block_size <= 0L) {
    return(NULL)
  }
  max(block_size, 1L)
}


#' Add a materialized scalar row cache to an analysis context
#'
#' The cache is intentionally chunk-local. In forked parallel execution this
#' lets the parent process perform the TileDB I/O once and workers reuse the
#' in-memory block through copy-on-write.
#' @noRd
.with_scalar_row_cache <- function(ctx, element.subset) {
  if (is.null(.scalar_read_block_size(ctx))) {
    return(ctx)
  }

  element.subset <- as.integer(element.subset)
  if (length(element.subset) == 0L) {
    return(ctx)
  }

  cache <- stats::setNames(
    vector("list", length(ctx$attached_scalars)),
    ctx$attached_scalars
  )
  for (sname in ctx$attached_scalars) {
    cache[[sname]] <- as.matrix(
      scalars(ctx$modelarray)[[sname]][element.subset, , drop = FALSE]
    )
  }

  ctx$scalar_row_cache <- cache
  ctx$scalar_row_cache_elements <- element.subset
  ctx$scalar_row_cache_contiguous <- length(element.subset) == 1L ||
    all(diff(element.subset) == 1L)
  ctx$scalar_row_cache_first <- element.subset[[1]]
  ctx$scalar_row_cache_position <- stats::setNames(
    seq_along(element.subset),
    as.character(element.subset)
  )
  ctx
}


#' Return a cached or storage-backed scalar row
#' @noRd
.read_scalar_row <- function(ctx, sname, i_element) {
  cache <- ctx$scalar_row_cache
  if (!is.null(cache) && !is.null(cache[[sname]])) {
    if (isTRUE(ctx$scalar_row_cache_contiguous)) {
      pos <- as.integer(i_element - ctx$scalar_row_cache_first + 1L)
    } else {
      pos <- ctx$scalar_row_cache_position[[as.character(i_element)]]
      if (is.null(pos)) {
        pos <- NA_integer_
      }
      pos <- unname(pos)
    }
    if (!is.na(pos) && pos >= 1L && pos <= nrow(cache[[sname]])) {
      return(as.vector(cache[[sname]][pos, ]))
    }
  }

  as.vector(scalars(ctx$modelarray)[[sname]][i_element, ])
}


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
  response_vals <- .read_scalar_row(ctx, ctx$scalar, i_element)

  # Start the validity mask with the response scalar
  valid_mask <- is.finite(response_vals)

  # Read and reorder additional attached scalars
  scalar_values <- list()
  scalar_values[[ctx$scalar]] <- response_vals

  other_scalars <- setdiff(ctx$attached_scalars, ctx$scalar)
  for (sname in other_scalars) {
    s_vals <- .read_scalar_row(ctx, sname, i_element)

    # Apply precomputed reorder index (or use as-is if NULL)
    reorder_idx <- ctx$predictor_reorder[[sname]]
    if (!is.null(reorder_idx)) {
      s_vals <- s_vals[reorder_idx]
    }

    scalar_values[[sname]] <- s_vals
    valid_mask <- valid_mask & is.finite(s_vals)
  }

  # Intersection mask across all scalars
  which_valid <- cheapr::which_(valid_mask)
  num_valid <- length(which_valid)

  if (!(num_valid > num.subj.lthr)) {
    return(list(dat = NULL, sufficient = FALSE, num_valid = num_valid))
  }

  # Build filtered data.frame
  dat <- cheapr::sset(ctx$phenotypes, which_valid)
  scalar_col_list <- lapply(
    stats::setNames(ctx$attached_scalars, ctx$attached_scalars),
    function(sname) scalar_values[[sname]][which_valid]
  )
  scalar_cols <- cheapr::fast_df(.args = scalar_col_list)
  dat <- cheapr::col_c(dat, scalar_cols)

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


#' Print progress context for chunked element dispatch
#' @noRd
.report_iteration_progress <- function(iteration, total_iterations,
                                       chunk_start, chunk_end,
                                       total_elements, pbar) {
  if (!isTRUE(pbar) || total_iterations <= 1L) {
    return(invisible(NULL))
  }
  message(
    "iteration ", iteration, "/", total_iterations,
    " (elements ", chunk_start, "-", chunk_end,
    " of ", total_elements, ")"
  )
  invisible(NULL)
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
#' Validate a gzip compression level used by a streaming writer
#' @noRd
.validate_stream_compression_level <- function(compression_level, argument) {
  if (!is.numeric(compression_level) ||
        length(compression_level) != 1L ||
        is.na(compression_level) ||
        !is.finite(compression_level) ||
        compression_level != floor(compression_level) ||
        compression_level < 0L ||
        compression_level > 9L) {
    stop(argument, " must be a single integer between 0 and 9", call. = FALSE)
  }
  as.integer(compression_level)
}


#' Create a dense TileDB values array for streaming
#' @noRd
.create_tiledb_values_array <- function(uri,
                                        n_rows,
                                        n_cols,
                                        chunk_rows,
                                        storage_mode,
                                        compression_level) {
  .require_tiledb_support()

  supported_modes <- c("logical", "integer", "double")
  if (!is.character(storage_mode) ||
        length(storage_mode) != 1L ||
        !storage_mode %in% supported_modes) {
    stop(
      "TileDB streaming supports storage modes: ",
      paste(supported_modes, collapse = ", "),
      call. = FALSE
    )
  }

  parent_uri <- dirname(uri)
  if (dir.exists(parent_uri) || .tiledb_object_exists(parent_uri)) {
    unlink(parent_uri, recursive = TRUE, force = TRUE)
  }
  dir.create(parent_uri, recursive = TRUE, showWarnings = FALSE)

  dimensions <- list(
    tiledb::tiledb_dim(
      name = "d1",
      domain = c(1L, as.integer(n_rows)),
      tile = as.integer(chunk_rows),
      type = "INT32"
    ),
    tiledb::tiledb_dim(
      name = "d2",
      domain = c(1L, as.integer(n_cols)),
      tile = as.integer(n_cols),
      type = "INT32"
    )
  )
  domain <- tiledb::tiledb_domain(dims = dimensions)

  gzip_filter <- tiledb::tiledb_filter("GZIP")
  gzip_filter <- tiledb::tiledb_filter_set_option(
    gzip_filter,
    "COMPRESSION_LEVEL",
    compression_level
  )
  value_filters <- tiledb::tiledb_filter_list(list(gzip_filter))
  value_attribute <- tiledb::tiledb_attr(
    name = "values",
    type = tiledb::r_to_tiledb_type(vector(storage_mode)),
    filter_list = value_filters
  )
  schema <- tiledb::tiledb_array_schema(
    domain = domain,
    attrs = list(value_attribute),
    cell_order = "COL_MAJOR",
    tile_order = "COL_MAJOR",
    sparse = FALSE
  )
  tiledb::tiledb_array_create(uri, schema)

  .write_tiledb_array_metadata(uri, list(type = storage_mode))
  invisible(NULL)
}


#' Write one matrix block to a dense TileDB values array
#' @noRd
.write_tiledb_values_block <- function(uri, row_idx, n_cols, block) {
  array <- tiledb::tiledb_array(uri)
  on.exit(try(tiledb::tiledb_array_close(array), silent = TRUE), add = TRUE)
  query <- tiledb::tiledb_query(array, "WRITE")
  query <- tiledb::tiledb_query_set_subarray(
    query,
    as.integer(c(range(row_idx), 1L, n_cols)),
    "INT32"
  )
  query <- tiledb::tiledb_query_set_layout(query, "COL_MAJOR")
  query <- tiledb::tiledb_query_set_buffer(
    query,
    "values",
    as.vector(block)
  )
  query <- tiledb::tiledb_query_submit(query)
  query <- tiledb::tiledb_query_finalize(query)
  if (!identical(tiledb::tiledb_query_status(query), "COMPLETE")) {
    stop("TileDB block write did not complete", call. = FALSE)
  }
  invisible(NULL)
}


#' Initialize an incremental writer for /scalars datasets
#' @noRd
.init_scalar_stream_writer <- function(write_scalar_name,
                                       write_scalar_file,
                                       n_rows,
                                       column_names,
                                       flush_every = 1000L,
                                       storage_mode = "double",
                                       compression_level = 4L) {
  compression_level <- .validate_stream_compression_level(
    compression_level,
    "write_scalar_compression_level"
  )
  if (!is.character(write_scalar_file) ||
        length(write_scalar_file) != 1L ||
        is.na(write_scalar_file) ||
        !nzchar(write_scalar_file)) {
    stop("write_scalar_file must be a non-empty character path")
  }
  backend <- .resolve_storage_backend(write_scalar_file, "auto")
  n_cols <- length(column_names)
  chunk_rows <- min(flush_every, n_rows)

  if (identical(backend, "tiledb")) {
    scalar_uri <- .storage_child_path(
      .modelarray_path(write_scalar_file, "scalars"),
      write_scalar_name,
      "write_scalar_name"
    )
    values_uri <- .modelarray_path(scalar_uri, "values")
    .create_tiledb_values_array(
      uri = values_uri,
      n_rows = n_rows,
      n_cols = n_cols,
      chunk_rows = chunk_rows,
      storage_mode = storage_mode,
      compression_level = compression_level
    )
    return(list(
      backend = "tiledb",
      file = write_scalar_file,
      scalar_name = write_scalar_name,
      values_uri = values_uri,
      column_names = as.character(column_names),
      storage_mode = storage_mode,
      n_cols = n_cols,
      write_row_cursor = 1L
    ))
  }

  if (!file.exists(write_scalar_file)) {
    rhdf5::h5createFile(write_scalar_file)
  }
  scalar_group <- paste0("scalars/", write_scalar_name)
  h5_write <- hdf5r::H5File$new(write_scalar_file, mode = "a")
  if (!h5_write$exists("scalars")) {
    h5_write$create_group("scalars")
  }
  scalars_group <- h5_write$open("scalars")
  if (scalars_group$exists(write_scalar_name)) {
    scalars_group$link_delete(write_scalar_name)
  }
  scalars_group$create_group(write_scalar_name)
  h5_write$close_all()

  dataset_path <- paste0(scalar_group, "/values")
  rhdf5::h5createDataset(
    file = write_scalar_file,
    dataset = dataset_path,
    dims = c(n_rows, n_cols),
    storage.mode = storage_mode,
    chunk = c(chunk_rows, n_cols),
    level = compression_level
  )

  list(
    backend = "hdf5",
    file = write_scalar_file,
    scalar_name = write_scalar_name,
    dataset_path = dataset_path,
    column_names = as.character(column_names),
    storage_mode = storage_mode,
    n_cols = n_cols,
    write_row_cursor = 1L
  )
}


#' Append one block to an incremental /scalars writer
#' @noRd
.scalar_stream_write_block <- function(writer, block) {
  if (is.null(writer)) {
    return(writer)
  }
  block <- as.matrix(block)
  if (ncol(block) != writer$n_cols) {
    stop("Scalar block does not match the initialized writer column count")
  }
  storage.mode(block) <- writer$storage_mode
  row_idx <- writer$write_row_cursor:(writer$write_row_cursor + nrow(block) - 1L)

  if (identical(writer$backend, "tiledb")) {
    .write_tiledb_values_block(
      uri = writer$values_uri,
      row_idx = row_idx,
      n_cols = writer$n_cols,
      block = block
    )
  } else {
    rhdf5::h5write(
      obj = block,
      file = writer$file,
      name = writer$dataset_path,
      index = list(row_idx, seq_len(writer$n_cols))
    )
  }

  writer$write_row_cursor <- writer$write_row_cursor + nrow(block)
  writer
}


#' Finalize an incremental /scalars writer
#' @noRd
.finalize_scalar_stream_writer <- function(writer) {
  if (is.null(writer)) {
    return(invisible(NULL))
  }
  if (identical(writer$backend, "tiledb")) {
    .write_tiledb_array_metadata(
      writer$values_uri,
      list(
        column_names = jsonlite::toJSON(
          writer$column_names,
          auto_unbox = TRUE
        )
      )
    )
    return(invisible(NULL))
  }

  .write_hdf5_attribute(
    writer$file,
    writer$dataset_path,
    "column_names",
    writer$column_names
  )
  rhdf5::h5write(
    obj = writer$column_names,
    file = writer$file,
    name = paste0("scalars/", writer$scalar_name, "/column_names")
  )
  rhdf5::h5closeAll()
  invisible(NULL)
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
  if (!is.character(write_results_name) ||
        length(write_results_name) != 1L ||
        is.na(write_results_name) ||
        !nzchar(write_results_name)) {
    stop("write_results_name must be a non-empty character string when provided")
  }
  if (is.null(write_results_file) ||
        !is.character(write_results_file) ||
        length(write_results_file) != 1L ||
        is.na(write_results_file) ||
        !nzchar(write_results_file)) {
    stop("write_results_file must be a single character path when write_results_name is provided")
  }
  if (!is.numeric(flush_every) || length(flush_every) != 1L || flush_every <= 0) {
    stop("write_results_flush_every must be a positive integer")
  }
  flush_every <- as.integer(flush_every)

  backend <- .resolve_storage_backend(write_results_file, "auto")
  if (identical(backend, "tiledb")) {
    .require_tiledb_support()
    analysis_uri <- .storage_child_path(
      .modelarray_path(write_results_file, "results"),
      write_results_name,
      "write_results_name"
    )
    compression_level <- .validate_stream_compression_level(
      compression_level,
      "write_results_compression_level"
    )
    values_uri <- .modelarray_path(analysis_uri, "results_matrix")
    .create_tiledb_values_array(
      uri = values_uri,
      n_rows = n_rows,
      n_cols = length(column_names),
      chunk_rows = min(flush_every, n_rows),
      storage_mode = storage_mode,
      compression_level = compression_level
    )
    return(list(
      backend = "tiledb",
      file = write_results_file,
      analysis_name = write_results_name,
      values_uri = values_uri,
      column_names = as.character(column_names),
      storage_mode = storage_mode,
      n_cols = length(column_names),
      write_row_cursor = 1L
    ))
  }

  if (!identical(backend, "hdf5")) {
    stop(
      "Unsupported write_results_file backend: ",
      backend,
      call. = FALSE
    )
  }

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
    backend = "hdf5",
    file = write_results_file,
    analysis_name = write_results_name,
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
  if (ncol(block_df) != writer$n_cols) {
    stop("Result block does not match the initialized writer column count")
  }
  block <- as.matrix(block_df)
  row_idx <- writer$write_row_cursor:(writer$write_row_cursor + nrow(block) - 1L)
  if (identical(writer$backend, "tiledb")) {
    storage.mode(block) <- writer$storage_mode
    .write_tiledb_values_block(
      uri = writer$values_uri,
      row_idx = row_idx,
      n_cols = writer$n_cols,
      block = block
    )
    writer$write_row_cursor <- writer$write_row_cursor + nrow(block)
    return(writer)
  }

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
  if (identical(writer$backend, "tiledb")) {
    .write_tiledb_array_metadata(
      writer$values_uri,
      list(
        column_names = jsonlite::toJSON(
          writer$column_names,
          auto_unbox = TRUE
        )
      )
    )
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
