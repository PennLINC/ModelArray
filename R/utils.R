#' check if an object in .h5 exists
#' @param fn_h5 filename of the .h5 file
#' @param group_name full directory of this object in .h5 name
#' @param object_name name of the object, should be a string without "/"
#' @noRd
flagObjectExistInh5 <- function(fn_h5, group_name = "/results", object_name = "myAnalysis") {
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  h5.nrow <- nrow(h5[h5$group == group_name & h5$name == object_name, ])
  h5.nrow > 0
}


#' check if h5 group "results" exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @noRd
flagResultsGroupExistInh5 <- function(fn_h5) {
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  h5.nrow <- nrow(h5[h5$group == "/" & h5$name == "results", ])
  h5.nrow > 0
}

#' check if a subfolder of results exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @param analysis_name The subfolder name in "results" in .h5 file
#' @noRd
flagAnalysisExistInh5 <- function(fn_h5, analysis_name) {
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  h5.nrow <- nrow(h5[h5$group == "/results" & h5$name == analysis_name, ])
  h5.nrow > 0
}


#' Write an HDF5 attribute with hdf5r's current API
#' @noRd
.write_hdf5_attribute <- function(filepath, object_path, attr_name, value) {
  h5 <- hdf5r::H5File$new(filepath, mode = "a")
  on.exit(h5$close_all(), add = TRUE)
  object <- h5[[object_path]]
  # create_attr() errors if the attribute is already there, unlike the
  # rhdf5::h5writeAttribute() it replaces, which overwrote silently.
  if (object$attr_exists(attr_name)) {
    object$attr_delete(attr_name)
  }
  object$create_attr(attr_name, robj = value)
  invisible(NULL)
}


#' Does a directory look like a ModelArray storage tree?
#'
#' Filesystem-only on purpose: this runs for HDF5 paths too, so it must not
#' reach for the optional TileDB packages.
#' @noRd
.looks_like_modelarray_store <- function(filepath) {
  dir.exists(file.path(filepath, "scalars")) ||
    dir.exists(file.path(filepath, "results"))
}


#' Resolve a storage backend from a path and user preference
#' @noRd
.resolve_storage_backend <- function(filepath, backend = c("auto", "hdf5", "tiledb")) {
  backend <- match.arg(backend)
  if (backend != "auto") {
    return(backend)
  }

  # A `.tdb` suffix declares intent, whether or not the store exists yet, so it
  # is also how a *new* TileDB store gets created without passing `backend`.
  if (grepl("\\.tdb/?$", filepath, ignore.case = TRUE)) {
    return("tiledb")
  }
  # Otherwise only claim a directory when it actually holds a store. Writers
  # delete the group they are about to replace, so a mistyped output directory
  # must not be mistaken for one.
  if (dir.exists(filepath) && .looks_like_modelarray_store(filepath)) {
    return("tiledb")
  }
  "hdf5"
}


#' Return the backend retained by a ModelArray object
#' @noRd
.modelarray_backend <- function(x, filepath = x@path[1]) {
  stored_backend <- tryCatch(
    methods::slot(x, "backend"),
    error = function(e) character(0)
  )

  if (length(stored_backend) > 1L && length(x@path) == length(stored_backend)) {
    path_index <- match(filepath, x@path)
    if (!is.na(path_index)) {
      stored_backend <- stored_backend[path_index]
    }
  }

  # unname(): mergeModelArrays() stores one backend per scalar name, and every
  # caller compares the result with identical(., "tiledb") — a surviving name
  # makes that FALSE, which would silently disable TileDB handling on any
  # merged object.
  stored_backend <- unname(stored_backend[1])
  if (
    length(stored_backend) == 0L ||
      is.na(stored_backend) ||
      !stored_backend %in% c("hdf5", "tiledb")
  ) {
    return(.resolve_storage_backend(filepath, "auto"))
  }
  stored_backend
}


#' Require TileDB packages for TileDB-backed operations
#' @noRd
.require_tiledb_support <- function() {
  missing <- c(
    if (!requireNamespace("tiledb", quietly = TRUE)) "tiledb",
    if (!requireNamespace("TileDBArray", quietly = TRUE)) "TileDBArray",
    if (!requireNamespace("jsonlite", quietly = TRUE)) "jsonlite"
  )
  if (length(missing) > 0) {
    stop(
      "TileDB support requires the optional package(s): ",
      paste(missing, collapse = ", "),
      ". Install them before using a TileDB-backed ModelArray.",
      call. = FALSE
    )
  }
}


#' Join path segments for local paths and TileDB URIs
#' @noRd
.modelarray_path <- function(...) {
  parts <- c(...)
  parts <- gsub("/+$", "", parts)
  parts <- gsub("^/+", "", parts)
  first <- c(...)[[1]]
  first <- gsub("/+$", "", first)
  paste(c(first, parts[-1]), collapse = "/")
}


#' Validate a storage object name as one safe path component
#' @noRd
.validate_storage_component <- function(value, argument) {
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !nzchar(value)
  ) {
    stop(argument, " must be a non-empty character string", call. = FALSE)
  }
  if (value %in% c(".", "..") || grepl("[/\\\\]", value)) {
    stop(
      argument,
      " must be a single storage name, cannot be '.' or '..', and cannot contain '/' or '\\'",
      call. = FALSE
    )
  }
  value
}


#' Build a storage child path and verify local-path containment
#' @noRd
.storage_child_path <- function(parent, component, argument) {
  component <- .validate_storage_component(component, argument)
  child <- .modelarray_path(parent, component)

  is_uri <- grepl("^[[:alpha:]][[:alnum:]+.-]*://", parent)
  if (!is_uri) {
    normalized_parent <- normalizePath(parent, winslash = "/", mustWork = FALSE)
    normalized_child <- if (file.exists(child) || dir.exists(child)) {
      normalizePath(child, winslash = "/", mustWork = FALSE)
    } else {
      file.path(normalized_parent, component)
    }
    parent_prefix <- paste0(sub("/+$", "", normalized_parent), "/")
    if (!startsWith(normalized_child, parent_prefix)) {
      stop(argument, " resolves outside its storage parent", call. = FALSE)
    }
  }

  child
}


#' Normalize vectors read from metadata or companion datasets
#' @noRd
.clean_name_vector <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.list(x) && !is.data.frame(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  x <- as.vector(x)
  x <- as.character(x)
  x <- gsub("[\\x00]+$", "", x, perl = TRUE, useBytes = TRUE)
  trimws(x)
}


#' Read names from HDF5 attributes, falling back to datasets
#' @noRd
.read_hdf5_names <- function(filepath, attr_path, attr_name, dataset_paths, label) {
  attrs <- rhdf5::h5readAttributes(filepath, name = attr_path)
  names_attr <- attrs[[attr_name]]
  if (!is.null(names_attr)) {
    return(.clean_name_vector(names_attr))
  }

  colnames_ds <- NULL
  last_error <- NULL
  for (p in dataset_paths) {
    tmp <- tryCatch(
      rhdf5::h5read(filepath, p),
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
      "Neither attribute '", attr_name, "' nor a dataset with ", label, " found. Tried: ",
      paste(dataset_paths, collapse = ", "),
      if (!is.null(last_error)) paste0(". Last error: ", conditionMessage(last_error)) else ""
    ))
  }
  .clean_name_vector(colnames_ds)
}


#' Return a TileDB object type, or INVALID when lookup fails
#' @noRd
.tiledb_object_type <- function(uri) {
  .require_tiledb_support()
  tryCatch(tiledb::tiledb_object_type(uri), error = function(e) "INVALID")
}


#' Does a TileDB object exist?
#' @noRd
.tiledb_object_exists <- function(uri, type = NULL) {
  object_type <- .tiledb_object_type(uri)
  if (is.null(type)) {
    return(toupper(object_type) %in% c("ARRAY", "GROUP"))
  }
  identical(toupper(object_type), toupper(type))
}


#' Convert TileDB metadata values to a character vector
#' @noRd
.parse_tiledb_names <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.raw(x)) {
    x <- rawToChar(x)
  }
  if (is.list(x) && length(x) == 1L) {
    x <- x[[1]]
  }
  x <- .clean_name_vector(x)
  if (length(x) == 1L && jsonlite::validate(x)) {
    x <- jsonlite::fromJSON(x)
  }
  .clean_name_vector(x)
}


#' Read metadata from a TileDB array
#' @noRd
.read_tiledb_array_metadata <- function(uri, key) {
  .require_tiledb_support()
  if (!.tiledb_object_exists(uri, "ARRAY")) {
    return(NULL)
  }
  arr <- tiledb::tiledb_array(uri, query_type = "READ", keep_open = TRUE)
  on.exit(try(tiledb::tiledb_array_close(arr), silent = TRUE), add = TRUE)
  tryCatch(tiledb::tiledb_get_metadata(arr, key), error = function(e) NULL)
}


#' Read metadata from a TileDB group
#' @noRd
.read_tiledb_group_metadata <- function(uri, key) {
  .require_tiledb_support()
  if (!.tiledb_object_exists(uri, "GROUP")) {
    return(NULL)
  }
  grp <- tiledb::tiledb_group(uri, type = "READ")
  on.exit(try(tiledb::tiledb_group_close(grp), silent = TRUE), add = TRUE)
  tryCatch(tiledb::tiledb_group_get_metadata(grp, key), error = function(e) NULL)
}


#' Read the "values" attribute from a 1-D TileDB array
#' @noRd
.read_tiledb_values_array <- function(uri) {
  .require_tiledb_support()
  if (!.tiledb_object_exists(uri, "ARRAY")) {
    return(NULL)
  }
  arr <- tiledb::tiledb_array(
    uri,
    query_type = "READ",
    attrs = "values",
    return_as = "data.frame"
  )
  on.exit(try(tiledb::tiledb_array_close(arr), silent = TRUE), add = TRUE)
  values <- tryCatch(arr[], error = function(e) NULL)
  if (is.null(values)) {
    return(NULL)
  }
  if (is.data.frame(values) && "values" %in% names(values)) {
    return(.clean_name_vector(values[["values"]]))
  }
  if (is.list(values) && "values" %in% names(values)) {
    return(.clean_name_vector(values[["values"]]))
  }
  .clean_name_vector(values)
}


#' Return dimensions for a dense TileDB array
#' @noRd
.tiledb_array_dims <- function(uri) {
  .require_tiledb_support()
  arr <- tiledb::tiledb_array(uri, query_type = "READ", keep_open = TRUE)
  on.exit(try(tiledb::tiledb_array_close(arr), silent = TRUE), add = TRUE)

  dom <- tiledb::domain(tiledb::schema(arr))
  ndim <- tiledb::tiledb_ndim(dom)
  vapply(seq_len(ndim), function(i) {
    extent <- tiledb::tiledb_array_get_non_empty_domain_from_index(arr, i)
    as.integer(as.numeric(extent[2]) - as.numeric(extent[1]) + 1)
  }, integer(1))
}


#' Read scalar column names for HDF5 or TileDB backends
#' @noRd
.read_scalar_column_names <- function(filepath, scalar_type, backend) {
  if (identical(backend, "hdf5")) {
    return(.read_hdf5_names(
      filepath = filepath,
      attr_path = sprintf("scalars/%s/values", scalar_type),
      attr_name = "column_names",
      dataset_paths = c(
        sprintf("scalars/%s/column_names", scalar_type),
        sprintf("scalars/%s/values/column_names", scalar_type),
        sprintf("scalars/scalars/%s/values/column_names", scalar_type),
        sprintf("scalars/scalars/%s/column_names", scalar_type)
      ),
      label = "column names"
    ))
  }

  .require_tiledb_support()
  value_uri <- .modelarray_path(filepath, "scalars", scalar_type, "values")
  scalar_uri <- .modelarray_path(filepath, "scalars", scalar_type)
  names <- .parse_tiledb_names(.read_tiledb_array_metadata(value_uri, "column_names"))
  if (is.null(names)) {
    names <- .parse_tiledb_names(.read_tiledb_group_metadata(scalar_uri, "column_names"))
  }
  if (is.null(names)) {
    names <- .read_tiledb_values_array(.modelarray_path(scalar_uri, "column_names"))
  }
  if (is.null(names)) {
    stop(
      "No TileDB column names found for scalar '", scalar_type, "'. ",
      "Expected column_names metadata on ", value_uri,
      " or a values array at ", .modelarray_path(scalar_uri, "column_names"),
      call. = FALSE
    )
  }
  names
}


#' Read result column names for HDF5 or TileDB backends
#' @noRd
.read_result_column_names <- function(filepath, analysis_name, backend) {
  if (identical(backend, "hdf5")) {
    return(.read_hdf5_names(
      filepath = filepath,
      attr_path = sprintf("results/%s/results_matrix", analysis_name),
      attr_name = "colnames",
      dataset_paths = c(
        sprintf("results/%s/column_names", analysis_name),
        sprintf("results/%s/results_matrix/column_names", analysis_name)
      ),
      label = "column names for results"
    ))
  }

  .require_tiledb_support()
  matrix_uri <- .modelarray_path(filepath, "results", analysis_name, "results_matrix")
  analysis_uri <- .modelarray_path(filepath, "results", analysis_name)
  names <- .parse_tiledb_names(.read_tiledb_array_metadata(matrix_uri, "colnames"))
  if (is.null(names)) {
    names <- .parse_tiledb_names(.read_tiledb_array_metadata(matrix_uri, "column_names"))
  }
  if (is.null(names)) {
    names <- .parse_tiledb_names(.read_tiledb_group_metadata(analysis_uri, "column_names"))
  }
  if (is.null(names)) {
    names <- .read_tiledb_values_array(.modelarray_path(analysis_uri, "column_names"))
  }
  if (is.null(names)) {
    stop("No TileDB column names found for analysis '", analysis_name, "'.", call. = FALSE)
  }
  names
}


#' Test for result groups/arrays in either supported backend
#' @noRd
.results_group_exists <- function(filepath, backend) {
  if (identical(backend, "hdf5")) {
    return(flagResultsGroupExistInh5(filepath))
  }
  .tiledb_object_exists(.modelarray_path(filepath, "results"), "GROUP") ||
    dir.exists(.modelarray_path(filepath, "results"))
}


#' Test for a named analysis in either supported backend
#' @noRd
.analysis_exists <- function(filepath, analysis_name, backend) {
  if (identical(backend, "hdf5")) {
    return(flagAnalysisExistInh5(filepath, analysis_name = analysis_name))
  }
  analysis_uri <- .modelarray_path(filepath, "results", analysis_name)
  matrix_uri <- .modelarray_path(analysis_uri, "results_matrix")
  .tiledb_object_exists(analysis_uri) ||
    dir.exists(analysis_uri) ||
    .tiledb_object_exists(matrix_uri, "ARRAY")
}


#' Test for a named object in either supported backend
#' @noRd
.stored_object_exists <- function(filepath, group_name, object_name, backend) {
  if (identical(backend, "hdf5")) {
    return(flagObjectExistInh5(filepath, group_name = group_name, object_name = object_name))
  }
  group_uri <- .modelarray_path(filepath, gsub("^/", "", group_name))
  object_uri <- .modelarray_path(group_uri, object_name)
  .tiledb_object_exists(object_uri, "ARRAY") ||
    !is.null(.read_tiledb_array_metadata(.modelarray_path(group_uri, "results_matrix"), object_name)) ||
    !is.null(.read_tiledb_group_metadata(group_uri, object_name))
}


#' Read a result look-up table for either supported backend
#' @noRd
.read_result_lut <- function(filepath, analysis_name, object_name, backend) {
  if (identical(backend, "hdf5")) {
    if (!flagObjectExistInh5(
      filepath,
      group_name = paste0("/results/", analysis_name),
      object_name = object_name
    )) {
      return(NULL)
    }
    return(ModelArraySeed(
      filepath,
      name = paste0("results/", analysis_name, "/", object_name),
      type = NA,
      backend = "hdf5"
    ) %>%
      DelayedArray::DelayedArray() %>%
      as.character())
  }

  analysis_uri <- .modelarray_path(filepath, "results", analysis_name)
  matrix_uri <- .modelarray_path(analysis_uri, "results_matrix")
  lut <- .parse_tiledb_names(.read_tiledb_array_metadata(matrix_uri, object_name))
  if (is.null(lut)) {
    lut <- .parse_tiledb_names(.read_tiledb_group_metadata(analysis_uri, object_name))
  }
  if (is.null(lut)) {
    lut <- .read_tiledb_values_array(.modelarray_path(analysis_uri, object_name))
  }
  lut
}


#' Recode one results matrix column from stored LUT indices to labels
#' @noRd
.recode_result_lut_column <- function(results_matrix, i_col, lut) {
  lut <- as.character(lut)
  stored_values <- as.character(as.vector(results_matrix[, i_col]))
  lut_index <- match(stored_values, as.character(seq_along(lut)))
  recoded_values <- stored_values
  recoded_values[!is.na(lut_index)] <- lut[lut_index[!is.na(lut_index)]]
  results_matrix[, i_col] <- recoded_values
  results_matrix
}


#' Coerce result data.frame columns to storage-friendly numeric columns
#' @noRd
.prepare_results_for_storage <- function(df.output) {
  luts <- list()
  for (i_col in seq_len(ncol(df.output))) {
    col <- df.output[[i_col]]
    not_numeric_or_int <- !is.numeric(col) && !is.integer(col)
    if (not_numeric_or_int) {
      message(
        paste0(
          "the column #",
          as.character(i_col),
          " of df.output to save: ",
          "data class is not numeric or integer...fixing it"
        )
      )

      factors <- factor(col)
      df.output[[i_col]] <- as.numeric(factors)
      luts[[paste0("lut_forcol", as.character(i_col))]] <- levels(factors)
    }
  }
  list(data = df.output, luts = luts)
}


#' Write metadata to a TileDB array, ignoring close failures
#' @noRd
.write_tiledb_array_metadata <- function(uri, values) {
  .require_tiledb_support()
  arr <- tiledb::tiledb_array(uri, query_type = "WRITE", keep_open = TRUE)
  on.exit(try(tiledb::tiledb_array_close(arr), silent = TRUE), add = TRUE)
  for (key in names(values)) {
    value <- values[[key]]
    if (inherits(value, "json")) {
      value <- as.character(value)
    }
    tiledb::tiledb_put_metadata(arr, key, value)
  }
  invisible(NULL)
}


#' Write ModelArray results to a TileDB store
#' @noRd
.write_results_tiledb <- function(fn.output, df.output, analysis_name, overwrite) {
  .require_tiledb_support()

  if (
    !is.character(fn.output) ||
      length(fn.output) != 1L ||
      is.na(fn.output) ||
      !nzchar(fn.output)
  ) {
    stop("fn.output must be a non-empty character path", call. = FALSE)
  }
  results_uri <- .modelarray_path(fn.output, "results")
  analysis_uri <- .storage_child_path(results_uri, analysis_name, "analysis_name")

  prepared <- .prepare_results_for_storage(df.output)
  df.output <- prepared$data

  if (!dir.exists(fn.output)) {
    dir.create(fn.output, recursive = TRUE, showWarnings = FALSE)
  }

  matrix_uri <- .modelarray_path(analysis_uri, "results_matrix")
  analysis_exists <- dir.exists(analysis_uri) || .tiledb_object_exists(analysis_uri)

  if (analysis_exists && !overwrite) {
    warning(paste0(analysis_name, " exists but not to overwrite!"))
    return(invisible(NULL))
  }
  if (analysis_exists && overwrite) {
    unlink(analysis_uri, recursive = TRUE, force = TRUE)
  }

  dir.create(analysis_uri, recursive = TRUE, showWarnings = FALSE)
  TileDBArray::writeTileDBArray(
    as.matrix(df.output),
    path = matrix_uri,
    attr = "values"
  )

  metadata <- c(
    list(column_names = jsonlite::toJSON(as.character(colnames(df.output)), auto_unbox = TRUE)),
    stats::setNames(
      lapply(prepared$luts, jsonlite::toJSON, auto_unbox = TRUE),
      names(prepared$luts)
    )
  )
  .write_tiledb_array_metadata(matrix_uri, metadata)
  invisible(NULL)
}


#' print the additional arguments settings
#' @param FUN The function, e.g. mgcv::gam, without "()"
#' @param argu_name The argument name of the function
#' @param dots list of additional arguments
#' @param message_default The message for default
#' @param message_usr_input The message describing user's input
#' @noRd
printAdditionalArgu <- function(FUN, argu_name, dots, message_default = NULL, message_usr_input = NULL) {
  dots_names <- names(dots)
  if (argu_name %in% dots_names) {
    if (is.null(message_default)) {
      message_default <- invisible(eval(formals(FUN)[[argu_name]]))
    }
    if (is.null(message_default)) { # if the default is NULL:
      message_default <- "NULL"
    }

    if (is.null(message_usr_input)) {
      m1 <- paste0(argu_name, " = ", dots[[argu_name]], " (default: ", message_default, ")") %>%
        crayon::black() %>%
        cat() # or, %>% message()
    } else { # specified the message:
      m1 <- paste0(argu_name, " = ", message_usr_input, " (default: ", message_default, ")") %>%
        crayon::black() %>%
        cat()
    }
  } else {
    m1 <- paste0(argu_name, ": default") %>%
      crayon::black() %>%
      cat()
  }

  cat(m1, "\n")
}


#' Check if the list of p-value correction methods are valid for a specific type of term/model.
#' Can be used for any statistical model. As long as the p.value to be correct is named as "p.value".
#'
#' @param correct.list The list of correction methods for this type of term/model
#' @param name.correct.list The name of the list of correction methods for this type of term/model
#' @param var.list The list of statistics to be saved for this type of term/model
#' @param name.var.list The name of the list of statistics to be saved for this type of term/model
#' @noRd
check_validity_correctPValue <- function(correct.list, name.correct.list,
                                         var.list, name.var.list) {
  p.adjust.methods.full <- stats::p.adjust.methods[stats::p.adjust.methods != "none"]
  if (all(correct.list == "none") == FALSE) { # any element is not "none"
    checker.method.in <- correct.list %in% p.adjust.methods.full
    if (all(checker.method.in) == FALSE) { # not all "TRUE"
      stop(
        paste0(
          "Some of elements in ",
          name.correct.list,
          " are not valid. Valid inputs are: ",
          paste(p.adjust.methods.full, collapse = ", ")
        )
      )
    }

    if ("p.value" %in% var.list == FALSE) { # not in the list | # check whether there is "p.value" in var.list
      warning(
        paste0(
          "p.value was not included in ",
          name.var.list,
          ", so not to perform its p.value corrections"
        ),
        immediate. = TRUE
      )
    }
  }
}


#' Print out important arguments in smooth terms s() in mgcv::gam() formula
#'
#' @details
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/s
#'
#' @param ofInterest got via: `gam.formula.breakdown <- mgcv::interpret.gam(formula)`;
#' `ofInterest <- gam.formula.breakdown$smooth.spec[[i]]`
#' @noRd
#'
checker_gam_s <- function(ofInterest) {
  FUN <- mgcv::s

  # add "by=?" back to term name:
  if (ofInterest$by == "NA") {
    term_name <- ofInterest$label
  } else {
    term_name <- paste0(substr(ofInterest$label, 1, nchar(ofInterest$label) - 1), ", by=", ofInterest$by, ")")
  }

  paste0(term_name, ": ") %>%
    crayon::black() %>%
    cat()

  ### k (or bs.dim):   # could be multiple values
  m1 <- invisible(eval(formals(FUN)[["k"]])) # default
  m2 <- ofInterest$bs.dim # could be a list of multiple values
  if ((length(unique(m2)) == 1) && (m1 %in% unique(m2))) { # default
    msg_k <- " (default)"
  } else {
    msg_k <- ""
  }

  paste0(
    "  k = ",
    paste(as.character(m2), collapse = ", "), msg_k, "; "
  ) %>%
    crayon::black() %>%
    cat()

  ### fx:
  m1 <- invisible(eval(formals(FUN)[["fx"]])) %>% as.character() # default
  m2 <- ofInterest$fixed # actual

  if (as.character(!as.logical(m1)) %in% m2) { # there is an opposite logical value in m2 (actual)
    msg_fx <- ""
  } else {
    msg_fx <- " (default)"
  }

  paste0("  fx = ", toString(ofInterest$fixed), msg_fx, "; ") %>%
    crayon::black() %>%
    cat()

  ### bs:
  m1 <- invisible(eval(formals(FUN)[["bs"]])) # default
  # actual:
  mybs <- gsub(".smooth.spec", "", class(ofInterest)) # if there are multiple elements in bs, mybs will be a list
  if ((length(unique(mybs)) == 1) && (m1 %in% unique(mybs))) { # default
    msg_bs <- " (default)"
  } else {
    msg_bs <- ""
  }

  paste0(
    "  bs = ",
    paste(as.character(mybs), collapse = ", "), msg_bs
  ) %>%
    crayon::black() %>%
    cat()

  cat("\n")
}

#' Print out important arguments in smooth term te() or ti() or t2() in mgcv::gam() formula
#'
#' @details
#' Why a separate function is needed for t(), cannot using s(): in ofInterest,
#' "fx" is "fx" for t(), but "fixed" for s() - so they are different.
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/te or /t2()
#'
#' @param FUN could be mgcv::te(), ti() or t2()
#' @param ofInterest got via: `gam.formula.breakdown <- mgcv::interpret.gam(formula)`;
#' `ofInterest <- gam.formula.breakdown$smooth.spec[[i]]`
#' @noRd
#'
checker_gam_t <- function(FUN, ofInterest) {
  paste0(ofInterest$label, ": ") %>%
    crayon::black() %>%
    cat()

  # t()'s interpret.gam()'s smooth.spec does not have k or bs.dim as s() does;
  # may provided in ofInterest$margin[[i-xx]]$bs.dim but not fully clear

  ### fx:
  m1 <- invisible(eval(formals(FUN)[["fx"]])) %>% as.character() # default
  m2 <- ofInterest$fx # actual

  if (as.character(!as.logical(m1)) %in% m2) { # there is an opposite logical value in m2 (actual)
    msg_fx <- ""
  } else {
    msg_fx <- " (default)"
  }

  paste0("  fx = ", toString(ofInterest$fx), msg_fx, "; ") %>%
    crayon::black() %>%
    cat()

  ### bs:  # also different way of extracting from s()'s
  m1 <- invisible(eval(formals(FUN)[["bs"]])) %>% as.character() # default

  mybs <- list() # actual, as a list
  for (i in seq_along(ofInterest$margin)) {
    temp <- gsub(".smooth.spec", "", ofInterest$margin[[i]] %>% class())
    mybs[i] <- temp
  }

  # then check if all elements are default value of bs
  if ((length(unique(mybs)) == 1) && (m1 %in% unique(mybs))) { # all elements are the same, and = default
    msg_bs <- " (default)"
  } else {
    msg_bs <- ""
  }
  # print out
  paste0(
    "  bs = ",
    paste(as.character(mybs), collapse = ", "), msg_bs
  ) %>%
    crayon::black() %>%
    cat()

  cat("\n")
}

#' A checker for formula in gam for ModelArray.gam()
#' @param formula The formula
#' @param gam.formula.breakdown Got from mgcv::interpret.gam(formula)
#' @param onemodel The model of one element got from mgcv::gam()
#' @noRd
#'
checker_gam_formula <- function(formula, gam.formula.breakdown, onemodel = NULL) {
  # print out the formula:
  temp <- formula %>% as.character()
  str_formula <- paste0(temp[2], " ~ ", temp[3])

  m1 <- paste0("The formula requested: ", str_formula) %>%
    crayon::black() %>%
    cat()
  cat(m1, "\n")


  if (length(gam.formula.breakdown$smooth.spec) != 0) { # if there is smooth term
    list_smooth_terms <- character(length(gam.formula.breakdown$smooth.spec))
    for (i_smoothTerm in seq_along(gam.formula.breakdown$smooth.spec)) {
      ofInterest <- gam.formula.breakdown$smooth.spec[[i_smoothTerm]]
      list_smooth_terms[i_smoothTerm] <- ofInterest$label

      # check what class of smooth term (s or ti or ?); then call the corresponding function
      # - throw out the message for important argument for this class
      # ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/smooth.terms
      smooth.class <- strsplit(ofInterest$label, "[(]")[[1]][1]
      if (smooth.class == "s") {
        checker_gam_s(ofInterest)
      } else if (smooth.class == "te") {
        checker_gam_t(mgcv::te, ofInterest)
      } else if (smooth.class == "ti") {
        checker_gam_t(mgcv::ti, ofInterest)
      } else if (smooth.class == "t2") {
        checker_gam_t(mgcv::t2, ofInterest)
      } else {
        stop(paste0("invalid smooth class for term ", ofInterest$label))
      }
    }
  } else { # no smooth term
    message("Warning: there is no smooth term in the requested formula")
  }


  # # fit for one element, get the summarized stat:
  # onemodel.tidy.smoothTerms <- onemodel %>% broom::tidy(parametric = FALSE)
  # onemodel.tidy.parametricTerms <- onemodel %>% broom::tidy(parametric = TRUE)
  # onemodel.glance <- onemodel %>% broom::glance()
  # onemodel.summary <- onemodel %>% summary()
}

#' Generate GAM formula with factor-smooth interaction
#'
#' @description
#' Generates a formula in the format
#' \code{y ~ orderedFactor + s(x) + s(x, by = orderedFactor)},
#' where \code{y} is \code{response.var}, \code{x} is \code{smooth.var},
#' and \code{orderedFactor} is \code{factor.var}.
#' The formula generated can be further modified, e.g. by adding covariates.
#'
#' @details
#' This helper exists because setting up factor-smooth interactions in
#' \code{\link[mgcv]{gam}} requires an ordered factor and a specific formula
#' structure. If \code{factor.var} in \code{phenotypes} is not already an
#' ordered factor, this function creates one using \code{reference.group} as
#' the baseline level and adds it as a new column (named with
#' \code{prefix.ordered.factor} prepended to \code{factor.var}).
#'
#' The returned \code{phenotypes} data.frame must be used in the subsequent
#' \code{\link{ModelArray.gam}} call so that the ordered factor column is
#' available to the model.
#'
#' @param response.var Character. The variable name for the response
#'   (dependent variable), typically a scalar name like \code{"FD"}.
#' @param factor.var Character. The variable name for the factor. It should
#'   be an ordered factor in \code{phenotypes}. If not, an ordered factor
#'   will be generated as a new column, which requires \code{reference.group}.
#' @param smooth.var Character. The variable name for the smooth term main
#'   effect (e.g. \code{"age"}).
#' @param phenotypes A data.frame of the cohort with columns of independent
#'   variables, including \code{factor.var} and \code{smooth.var}.
#' @param reference.group Character. The reference (baseline) group for the
#'   ordered factor of \code{factor.var}. Required when \code{factor.var}
#'   in \code{phenotypes} is not already an ordered factor.
#' @param prefix.ordered.factor Character. Prefix for the ordered factor
#'   column name. Required when \code{factor.var} in \code{phenotypes} is
#'   not already an ordered factor. Default is \code{"o"}.
#' @param fx Logical. Passed to \code{\link[mgcv]{s}}. If \code{TRUE}
#'   (recommended), the smooth is treated as fixed degrees of freedom.
#'   Default is \code{TRUE}.
#' @param k Integer or \code{NULL}. Basis dimension passed to
#'   \code{\link[mgcv]{s}} for both the main smooth and interaction terms.
#'   If \code{NULL} (default), uses the default from \code{mgcv::s()}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{formula}{The generated \code{\link[stats]{formula}} object.}
#'     \item{phenotypes}{The (possibly updated) data.frame. If
#'       \code{factor.var} was not already an ordered factor, a new column
#'       named \code{paste0(prefix.ordered.factor, factor.var)} is added.
#'       Otherwise identical to the input.}
#'   }
#'
#' @seealso \code{\link{gen_gamFormula_contIx}} for continuous-by-continuous
#'   interactions, \code{\link{ModelArray.gam}} which accepts the generated
#'   formula.
#'
#' @examples
#' \dontrun{
#' phenotypes <- read.csv("cohort.csv")
#'
#' # factor.var is not yet ordered - function creates it
#' result <- gen_gamFormula_fxSmooth(
#'   response.var = "FD",
#'   factor.var = "sex",
#'   smooth.var = "age",
#'   phenotypes = phenotypes,
#'   reference.group = "female"
#' )
#' result$formula
#'
#' # Use the updated phenotypes (contains the ordered factor column)
#' results <- ModelArray.gam(
#'   result$formula,
#'   data = ma,
#'   phenotypes = result$phenotypes,
#'   scalar = "FD"
#' )
#' }
#'
#' @rdname gen_gamFormula_fxSmooth
#' @export
#'
gen_gamFormula_fxSmooth <- function(response.var, factor.var, smooth.var, phenotypes,
                                    reference.group = NULL, prefix.ordered.factor = "o",
                                    fx = TRUE, k = NULL) {
  class.factor.var <- class(phenotypes[[factor.var]])
  if (
    !(
      (length(class.factor.var) == 2) &&
        (class.factor.var[1] == "ordered") &&
        (class.factor.var[2] == "factor")
    )
  ) { # class is not c("ordered", "factor")

    message(
      "input `factor.var` is not an ordered factor; will generate",
      " one in data.frame `phenotypes` which will be returned"
    )
    if (is.null(reference.group)) {
      stop("requires a reference.group to generate the ordered factor")
    }

    # name of the ordered factor:
    unordered.factor.var <- factor.var
    factor.var <- paste0(prefix.ordered.factor, unordered.factor.var)

    message(paste0("the ordered factor will be named as: ", factor.var))

    # check if factor.var already exists:
    if (factor.var %in% colnames(phenotypes)) {
      stop(
        paste0(
          "a column with the same name '",
          factor.var,
          "' already exists in `phenotypes` data frame! Please change the `prefix.ordered.factor`!"
        )
      )
    }

    # add the column to phenotypes:
    list.groups <- unique(phenotypes[[unordered.factor.var]])
    list.groups <- list.groups[list.groups != reference.group] # temporarily drop the reference group
    list.groups <- c(reference.group, list.groups) # add as first
    phenotypes[[factor.var]] <- ordered(phenotypes[[unordered.factor.var]], levels = list.groups)
    # the first element in the list.groups would be the reference group
  }

  if (is.null(k)) {
    k <- invisible(eval(formals(mgcv::s)[["k"]]))
  }

  # generate the formula:
  formula <- paste0(response.var, "~", factor.var, "+")
  formula <- paste0(formula, "s(", smooth.var, ",k=", toString(k), ",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "s(", smooth.var, ",k=", toString(k), ",by=", factor.var, ",fx=", toString(fx), ")")

  # return
  formula <- stats::as.formula(formula) # when printing formula, it could be shown in more than one lines...
  toReturn <- list(
    formula = formula,
    phenotypes = phenotypes
  )
  return(toReturn)
}


#' Generate GAM formula with continuous-by-continuous interaction
#'
#' @description
#' Generates a formula in the format
#' \code{y ~ ti(x) + ti(z) + ti(x, z)}, where \code{y} is
#' \code{response.var}, \code{x} is \code{cont1.var}, and \code{z} is
#' \code{cont2.var}. The formula generated can be further modified, e.g.
#' by adding covariates.
#'
#' @details
#' This helper uses \code{\link[mgcv]{ti}} (tensor product interaction)
#' terms so that the interaction \code{ti(x, z)} captures only the
#' interaction effect, separate from the main effects \code{ti(x)} and
#' \code{ti(z)}. This decomposition is important for interpretability and
#' for requesting \code{changed.rsq.term.index} in
#' \code{\link{ModelArray.gam}}.
#'
#' @param response.var Character. The variable name for the response
#'   (dependent variable), typically a scalar name like \code{"FD"}.
#' @param cont1.var Character. The name of the first continuous variable.
#' @param cont2.var Character. The name of the second continuous variable.
#' @param fx Logical. Passed to \code{\link[mgcv]{ti}}. If \code{TRUE}
#'   (recommended), the smooth is treated as fixed degrees of freedom.
#'   Default is \code{TRUE}.
#' @param k Integer or \code{NULL}. Basis dimension passed to
#'   \code{\link[mgcv]{ti}} for all three terms (both main effects and the
#'   interaction). If \code{NULL} (default), uses the default from
#'   \code{mgcv::ti()}.
#'
#' @return A \code{\link[stats]{formula}} object.
#'
#' @seealso \code{\link{gen_gamFormula_fxSmooth}} for factor-smooth
#'   interactions, \code{\link{ModelArray.gam}} which accepts the generated
#'   formula.
#'
#' @examples
#' \dontrun{
#' formula <- gen_gamFormula_contIx(
#'   response.var = "FD",
#'   cont1.var = "age",
#'   cont2.var = "cognition"
#' )
#' formula
#'
#' # Use in ModelArray.gam with changed R-squared for the interaction
#' results <- ModelArray.gam(
#'   formula,
#'   data = ma,
#'   phenotypes = phenotypes,
#'   scalar = "FD",
#'   changed.rsq.term.index = list(3)
#' )
#' }
#'
#' @rdname gen_gamFormula_contIx
#' @export
gen_gamFormula_contIx <- function(response.var, cont1.var, cont2.var,
                                  fx = TRUE, k = NULL) {
  if (is.null(k)) {
    k <- invisible(eval(formals(mgcv::ti)[["k"]]))
  }

  formula <- paste0(response.var, "~")
  formula <- paste0(formula, "ti(", cont1.var, ",k=", toString(k), ",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "ti(", cont2.var, ",k=", toString(k), ",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "ti(", cont1.var, ",", cont2.var, ",k=", toString(k), ",fx=", toString(fx), ")")

  formula <- stats::as.formula(formula)
  return(formula)
}


#' Bind two tibbles together, considering one or both of them is empty tibble()
#' @details Without this function, bind_cols(A, tibble()) will becomes tibble()...
#' @param a A tibble, can be empty tibble()
#' @param b A tibble, can be empty tibble()
#' @return c, A tibble after binding a and b together
#' @noRd
bind_cols_check_emptyTibble <- function(a, b) {
  flag_a_empty <- !all(dim(a)) # if TRUE, a is empty
  flag_b_empty <- !all(dim(b)) # if TRUE, b is empty

  if (flag_a_empty && flag_b_empty) c <- tibble::tibble() # both are empty
  if (flag_a_empty && (!flag_b_empty)) c <- b # b is not empty ==> taking b
  if ((!flag_a_empty) && flag_b_empty) c <- a # a is not empty ==> taking a
  if ((!flag_a_empty) && (!flag_b_empty)) c <- dplyr::bind_cols(a, b) # neither of them is empty ==> simply combine

  c
}


#' Summarize ModelArray storage without loading a full ModelArray
#'
#' @description
#' Reads the storage structure and returns a summary of available scalars,
#' their dimensions, and any saved analyses. Useful for inspecting large files
#' without constructing a full \linkS4class{ModelArray} object.
#'
#' @details
#' For HDF5, this function opens the file read-only via
#' \code{\link[rhdf5]{h5ls}}. For TileDB, it inspects arrays under
#' \code{/scalars/} and \code{/results/}. It does not load full data matrices
#' into memory. The returned
#' object has a \code{print} method that displays a formatted summary.
#'
#' \code{h5summary()} is a backward-compatible alias for
#' \code{ModelArraySummary()}. Despite the historical name, it supports both
#' HDF5 files and TileDB stores.
#'
#' @param filepath Character. Path to an HDF5 (\code{.h5}) file or TileDB
#'   (\code{.tdb}) store.
#' @param backend Character. Storage backend: \code{"auto"} (default),
#'   \code{"hdf5"}, or \code{"tiledb"}.
#'
#' @return \code{ModelArraySummary()} returns an object of class
#'   \code{c("ModelArraySummary", "h5summary")}. \code{h5summary()} returns
#'   an object of class \code{"h5summary"} for backward compatibility. Both
#'   are lists with components:
#'   \describe{
#'     \item{scalars}{A data.frame with columns \code{name},
#'       \code{nElements}, and \code{nInputFiles}.}
#'     \item{analyses}{Character vector of analysis names found under
#'       \code{/results/}.}
#'     \item{filepath}{The input filepath.}
#'     \item{backend}{The resolved storage backend.}
#'   }
#'
#' @seealso \code{\link{ModelArray}} for loading the full object,
#'   \linkS4class{ModelArray} for the class definition.
#'
#' @examples
#' \dontrun{
#' ModelArraySummary("path/to/data.h5")
#' ModelArraySummary("path/to/store.tdb")
#'
#' # Inspect before deciding which scalars to load
#' info <- ModelArraySummary("path/to/data.h5")
#' info$scalars$name
#' ma <- ModelArray("path/to/data.h5", scalar_types = info$scalars$name)
#'
#' # Historical alias, still supported
#' h5summary("path/to/data.h5")
#' }
#'
#' @rdname ModelArraySummary
#' @export
ModelArraySummary <- function(filepath, backend = c("auto", "hdf5", "tiledb")) {
  structure(
    .modelarray_summary_impl(filepath, backend = backend),
    class = c("ModelArraySummary", "h5summary")
  )
}

#' @rdname ModelArraySummary
#' @export
h5summary <- function(filepath, backend = c("auto", "hdf5", "tiledb")) {
  structure(
    .modelarray_summary_impl(filepath, backend = backend),
    class = "h5summary"
  )
}

#' @noRd
.modelarray_summary_impl <- function(filepath, backend = c("auto", "hdf5", "tiledb")) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  backend <- .resolve_storage_backend(filepath, backend)
  if (identical(backend, "tiledb")) {
    .require_tiledb_support()
    scalars_dir <- .modelarray_path(filepath, "scalars")
    scalar_candidates <- if (dir.exists(scalars_dir)) {
      basename(list.dirs(scalars_dir, recursive = FALSE, full.names = TRUE))
    } else {
      character(0)
    }
    scalar_names <- scalar_candidates[vapply(
      scalar_candidates,
      function(scalar_name) {
        .tiledb_object_exists(
          .modelarray_path(filepath, "scalars", scalar_name, "values"),
          "ARRAY"
        )
      },
      logical(1)
    )]
    scalar_info <- do.call(rbind, lapply(scalar_names, function(scalar_name) {
      values_uri <- .modelarray_path(filepath, "scalars", scalar_name, "values")
      dims <- .tiledb_array_dims(values_uri)
      source_names <- .read_scalar_column_names(filepath, scalar_name, backend)
      n_input_files <- length(source_names)
      if (dims[2] == n_input_files) {
        n_elements <- dims[1]
      } else if (dims[1] == n_input_files) {
        n_elements <- dims[2]
      } else {
        stop(
          "TileDB scalar '", scalar_name, "' dimensions ",
          paste(dims, collapse = " x "),
          " do not match its ", n_input_files, " stored column names",
          call. = FALSE
        )
      }
      data.frame(
        name = scalar_name,
        nElements = as.integer(n_elements),
        nInputFiles = as.integer(n_input_files),
        stringsAsFactors = FALSE
      )
    }))
    if (is.null(scalar_info)) {
      scalar_info <- data.frame(
        name = character(0),
        nElements = integer(0),
        nInputFiles = integer(0),
        stringsAsFactors = FALSE
      )
    }

    results_dir <- .modelarray_path(filepath, "results")
    analysis_candidates <- if (dir.exists(results_dir)) {
      basename(list.dirs(results_dir, recursive = FALSE, full.names = TRUE))
    } else {
      character(0)
    }
    analyses <- analysis_candidates[vapply(
      analysis_candidates,
      function(analysis_name) {
        .tiledb_object_exists(
          .modelarray_path(filepath, "results", analysis_name, "results_matrix"),
          "ARRAY"
        )
      },
      logical(1)
    )]

    return(
      list(
        scalars = scalar_info,
        analyses = analyses,
        filepath = filepath,
        backend = backend
      )
    )
  }

  listing <- rhdf5::h5ls(filepath)
  rhdf5::h5closeAll()

  # Find scalar datasets: /scalars/<name>/values
  scalar_rows <- listing[
    listing$group != "/scalars" &
      grepl("^/scalars/", listing$group) &
      !grepl("^/scalars/scalars", listing$group) &
      listing$name == "values",
  ]
  scalar_info <- data.frame(
    name = sub("^/scalars/", "", scalar_rows$group),
    dim = scalar_rows$dim,
    stringsAsFactors = FALSE
  )

  # Parse dimensions
  if (nrow(scalar_info) > 0) {
    dims <- strsplit(scalar_info$dim, " x ")
    scalar_info$nElements <- as.integer(sapply(dims, `[`, 1))
    scalar_info$nInputFiles <- as.integer(sapply(dims, `[`, 2))
    scalar_info$dim <- NULL
  }

  # Find analyses: /results/<name>
  result_groups <- listing[
    listing$group == "/results" & listing$otype == "H5I_GROUP",
  ]
  analyses <- result_groups$name

  list(
    scalars = scalar_info,
    analyses = analyses,
    filepath = filepath,
    backend = backend
  )
}

#' @rdname ModelArraySummary
#'
#' @param x A \code{ModelArraySummary} or \code{h5summary} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisible \code{x}. Called for its side effect of printing a
#'   human-readable summary to the console.
#'
#' @method print ModelArraySummary
#' @export
print.ModelArraySummary <- function(x, ...) {
  .print_modelarray_summary(x, ...)
}

#' @rdname ModelArraySummary
#' @method print h5summary
#' @export
print.h5summary <- function(x, ...) {
  .print_modelarray_summary(x, ...)
}

#' @noRd
.print_modelarray_summary <- function(x, ...) {
  backend <- if (is.null(x$backend)) "hdf5" else x$backend
  label <- if (identical(backend, "tiledb")) "TileDB store" else "H5 file"
  cat(label, ": ", x$filepath, "\n\n", sep = "")
  if (nrow(x$scalars) > 0) {
    cat("Scalars:\n")
    for (i in seq_len(nrow(x$scalars))) {
      cat("  ", x$scalars$name[i], ": ",
        x$scalars$nElements[i], " elements x ",
        x$scalars$nInputFiles[i], " input files\n",
        sep = ""
      )
    }
  } else {
    cat("  (no scalars found)\n")
  }
  cat("\nAnalyses:", if (length(x$analyses) > 0) {
    paste(x$analyses, collapse = ", ")
  } else {
    "(none)"
  }, "\n")
  invisible(x)
}
