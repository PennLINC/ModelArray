skip_if_no_tiledb <- function() {
  testthat::skip_if_not_installed("tiledb")
  testthat::skip_if_not_installed("TileDBArray")
  testthat::skip_if_not_installed("jsonlite")
}

write_tiledb_group_metadata <- function(uri, key, value) {
  dir.create(uri, recursive = TRUE, showWarnings = FALSE)
  if (!identical(toupper(tiledb::tiledb_object_type(uri)), "GROUP")) {
    tiledb::tiledb_group_create(uri)
  }
  grp <- tiledb::tiledb_group(uri, type = "WRITE")
  on.exit(try(tiledb::tiledb_group_close(grp), silent = TRUE), add = TRUE)
  tiledb::tiledb_group_put_metadata(
    grp,
    key,
    as.character(jsonlite::toJSON(value, auto_unbox = TRUE))
  )
}

write_tiledb_scalar <- function(store, scalar, values, sources, metadata_location = c("array", "group")) {
  metadata_location <- match.arg(metadata_location)
  dir.create(file.path(store, "scalars", scalar), recursive = TRUE, showWarnings = FALSE)
  uri <- file.path(store, "scalars", scalar, "values")
  TileDBArray::writeTileDBArray(values, path = uri, attr = "values")

  dir.create(file.path(store, "scalars", "__group"), showWarnings = FALSE)
  dir.create(file.path(store, "scalars", "__meta"), showWarnings = FALSE)

  if (identical(metadata_location, "array")) {
    arr <- tiledb::tiledb_array(uri, query_type = "WRITE", keep_open = TRUE)
    on.exit(try(tiledb::tiledb_array_close(arr), silent = TRUE), add = TRUE)
    tiledb::tiledb_put_metadata(
      arr,
      "column_names",
      as.character(jsonlite::toJSON(sources, auto_unbox = TRUE))
    )
  } else {
    write_tiledb_group_metadata(file.path(store, "scalars", scalar), "column_names", sources)
  }
  invisible(uri)
}

write_tiledb_result_with_group_metadata <- function(store, analysis_name, matrix_values, column_names, luts = list()) {
  analysis_uri <- file.path(store, "results", analysis_name)
  matrix_uri <- file.path(analysis_uri, "results_matrix")
  dir.create(analysis_uri, recursive = TRUE, showWarnings = FALSE)
  TileDBArray::writeTileDBArray(matrix_values, path = matrix_uri, attr = "values")
  write_tiledb_group_metadata(analysis_uri, "column_names", column_names)
  for (lut_name in names(luts)) {
    write_tiledb_group_metadata(analysis_uri, lut_name, luts[[lut_name]])
  }
  invisible(matrix_uri)
}

test_that("ModelArray reads TileDB scalar arrays written subjects by elements", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  values <- matrix(seq_len(8), nrow = 2, byrow = TRUE)
  sources <- c("sub-1.nii.gz", "sub-2.nii.gz")
  write_tiledb_scalar(store, "FD", values, sources)

  ma <- ModelArray(store, scalar_types = "FD")

  expect_equal(sources(ma)[["FD"]], sources)
  loaded <- as.matrix(scalars(ma)[["FD"]])
  dimnames(loaded) <- NULL
  expect_equal(loaded, t(values))
})

test_that("ModelArraySummary summarizes TileDB scalar arrays", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  write_tiledb_scalar(
    store,
    "FA",
    matrix(seq_len(12), nrow = 3),
    paste0("sub-", 1:3)
  )

  info <- ModelArraySummary(store, backend = "tiledb")

  expect_s3_class(info, "ModelArraySummary")
  expect_s3_class(info, "h5summary")
  expect_equal(info$backend, "tiledb")
  expect_equal(info$scalars$name, "FA")
  expect_equal(info$scalars$nInputFiles, 3L)
  expect_equal(info$scalars$nElements, 4L)

  legacy <- h5summary(store, backend = "tiledb")
  expect_s3_class(legacy, "h5summary")
  expect_false(inherits(legacy, "ModelArraySummary"))
  expect_equal(legacy$scalars, info$scalars)
  expect_equal(legacy$backend, "tiledb")
})

test_that("ModelArraySummary agrees for fake HDF5 and TileDB stores", {
  skip_if_no_tiledb()

  h5_store <- tempfile(fileext = ".h5")
  tiledb_store <- tempfile(fileext = ".tdb")
  on.exit(unlink(h5_store), add = TRUE)
  on.exit(unlink(tiledb_store, recursive = TRUE), add = TRUE)

  values <- matrix(seq_len(6), nrow = 3, ncol = 2)
  source_names <- c("sub-1", "sub-2")
  for (store in c(h5_store, tiledb_store)) {
    writer <- .init_scalar_stream_writer(
      write_scalar_name = "FA",
      write_scalar_file = store,
      n_rows = nrow(values),
      column_names = source_names,
      flush_every = 2L
    )
    writer <- .scalar_stream_write_block(writer, values[1:2, , drop = FALSE])
    writer <- .scalar_stream_write_block(writer, values[3, , drop = FALSE])
    .finalize_scalar_stream_writer(writer)

    writeResults(
      store,
      data.frame(element_id = 0:2, estimate = c(0.1, 0.2, 0.3)),
      analysis_name = "fake_analysis"
    )
  }

  h5_info <- ModelArraySummary(h5_store, backend = "hdf5")
  tiledb_info <- ModelArraySummary(tiledb_store, backend = "tiledb")

  expect_equal(
    h5_info[c("scalars", "analyses")],
    tiledb_info[c("scalars", "analyses")]
  )
  expect_equal(h5_info$scalars$nElements, 3L)
  expect_equal(h5_info$scalars$nInputFiles, 2L)
})

test_that("writeResults can round-trip TileDB results metadata and LUTs", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  write_tiledb_scalar(
    store,
    "FD",
    matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
    c("s1", "s2")
  )

  df <- data.frame(
    element_id = 0:1,
    beta = c(0.25, 0.5),
    gamma = c(1.25, 1.5),
    label = c("low", "high")
  )
  writeResults(store, df, analysis_name = "lm_tile", backend = "tiledb")

  ma <- ModelArray(
    store,
    scalar_types = "FD",
    analysis_names = "lm_tile",
    backend = "tiledb"
  )
  result_matrix <- results(ma, "lm_tile")$results_matrix
  expect_s4_class(result_matrix, "DelayedMatrix")
  expect_false(is.matrix(result_matrix))
  loaded <- as.matrix(result_matrix)

  expect_equal(colnames(loaded), colnames(df))
  expect_equal(as.numeric(loaded[, "element_id"]), df$element_id)
  expect_equal(as.numeric(loaded[, "beta"]), df$beta)
  expect_equal(as.numeric(loaded[, "gamma"]), df$gamma)
  expect_equal(loaded[, "label"], df$label)
})

test_that("TileDB writers reject unsafe names without deleting stored data", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  dir.create(file.path(store, "scalars"), recursive = TRUE)
  dir.create(file.path(store, "results"), recursive = TRUE)
  scalar_marker <- file.path(store, "scalars", "keep.txt")
  result_marker <- file.path(store, "results", "keep.txt")
  writeLines("keep", scalar_marker)
  writeLines("keep", result_marker)
  result_data <- data.frame(element_id = 0, estimate = 0.25)

  unsafe_names <- c("", ".", "..", "../scalars", "nested/name", "nested\\name")
  for (unsafe_name in unsafe_names) {
    expect_error(
      writeResults(
        store,
        result_data,
        analysis_name = unsafe_name,
        backend = "tiledb"
      ),
      "analysis_name"
    )
    expect_true(file.exists(scalar_marker))
    expect_true(file.exists(result_marker))
  }

  expect_error(
    .init_results_stream_writer(
      write_results_name = "../scalars",
      write_results_file = store,
      n_rows = 1L,
      column_names = names(result_data)
    ),
    "write_results_name"
  )
  expect_true(file.exists(scalar_marker))
  expect_true(file.exists(result_marker))

  expect_error(
    .init_scalar_stream_writer(
      write_scalar_name = "../results",
      write_scalar_file = store,
      n_rows = 1L,
      column_names = "estimate"
    ),
    "write_scalar_name"
  )
  expect_true(file.exists(scalar_marker))
  expect_true(file.exists(result_marker))
})

test_that("ModelArray reads TileDB scalar and result metadata from groups", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  write_tiledb_scalar(
    store,
    "FD",
    matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
    c("s1", "s2"),
    metadata_location = "group"
  )
  write_tiledb_result_with_group_metadata(
    store,
    "group_meta",
    matrix(c(0, 1, 1, 2), nrow = 2, byrow = TRUE),
    c("element_id", "label"),
    luts = list(lut_forcol2 = c("low", "high"))
  )

  ma <- ModelArray(
    store,
    scalar_types = "FD",
    analysis_names = "group_meta",
    backend = "tiledb"
  )
  loaded <- as.matrix(results(ma, "group_meta")$results_matrix)

  expect_equal(sources(ma)[["FD"]], c("s1", "s2"))
  expect_equal(colnames(loaded), c("element_id", "label"))
  expect_equal(loaded[, "label"], c("low", "high"))
})

test_that("analysis functions can consume TileDB-backed scalar arrays", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  src <- paste0("s", 1:8)
  write_tiledb_scalar(
    store,
    "FD",
    matrix(c(
      1, 3, 2, 7, 5, 9, 8, 6,
      2, 1, 4, 3, 6, 5, 7, 8
    ), nrow = 2, byrow = TRUE),
    src
  )
  ma <- ModelArray(store, scalar_types = "FD", backend = "tiledb")
  phen <- data.frame(source_file = src, age = seq(10, 80, by = 10))

  lm_out <- ModelArray.lm(
    FD ~ age,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = as.integer(c(1, 2)),
    var.terms = c("estimate"),
    var.model = c("adj.r.squared"),
    correct.p.value.terms = c("none"),
    correct.p.value.model = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )
  gam_out <- ModelArray.gam(
    FD ~ age,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = as.integer(c(1, 2)),
    var.smoothTerms = c(),
    var.parametricTerms = c("estimate"),
    var.model = c("dev.expl"),
    correct.p.value.smoothTerms = c("none"),
    correct.p.value.parametricTerms = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )
  wrap_out <- ModelArray.wrap(
    FUN = function(data) data.frame(mean_fd = mean(data$FD)),
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = as.integer(c(1, 2)),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )

  expect_equal(nrow(lm_out), 2)
  expect_equal(nrow(gam_out), 2)
  expect_equal(wrap_out$mean_fd, rowMeans(as.matrix(scalars(ma, "FD"))))
})

test_that("TileDB result streaming writes each block without retaining it", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  writer <- .init_results_stream_writer(
    write_results_name = "direct_stream",
    write_results_file = store,
    n_rows = 3L,
    column_names = c("element_id", "estimate"),
    flush_every = 1L
  )

  expect_false("blocks" %in% names(writer))
  writer <- .results_stream_write_block(
    writer,
    data.frame(element_id = 0, estimate = 0.25)
  )
  expect_false("blocks" %in% names(writer))

  fragment_info <- tiledb::tiledb_fragment_info(writer$values_uri)
  expect_equal(tiledb::tiledb_fragment_info_get_num(fragment_info), 1L)

  writer <- .results_stream_write_block(
    writer,
    data.frame(
      element_id = c(1, 2),
      estimate = c(0.5, NaN)
    )
  )
  .finalize_results_stream_writer(writer)

  fragment_info <- tiledb::tiledb_fragment_info(writer$values_uri)
  expect_equal(tiledb::tiledb_fragment_info_get_num(fragment_info), 2L)
})

test_that("TileDB streamed outputs exactly match regular lm, gam, and wrap outputs", {
  skip_if_no_tiledb()

  src <- paste0("s", 1:8)
  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  write_tiledb_scalar(
    store,
    "FD",
    matrix(c(
      1, 3, 2, 7, 5, 9, 8, 6,
      2, 1, 4, 3, 6, 5, 7, 8,
      rep(NA_real_, 8)
    ), nrow = 3, byrow = TRUE),
    src
  )
  ma <- ModelArray(store, scalar_types = "FD", backend = "tiledb")
  phen <- data.frame(source_file = src, age = seq(10, 80, by = 10))
  element_subset <- as.integer(c(3, 1, 2))

  expected_lm <- ModelArray.lm(
    FD ~ age,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    var.terms = c("estimate", "std.error", "statistic", "p.value"),
    var.model = c("adj.r.squared", "AIC"),
    correct.p.value.terms = c("none"),
    correct.p.value.model = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )

  out_lm <- ModelArray.lm(
    FD ~ age,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    var.terms = c("estimate", "std.error", "statistic", "p.value"),
    var.model = c("adj.r.squared", "AIC"),
    correct.p.value.terms = c("none"),
    correct.p.value.model = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE,
    write_results_name = "lm_stream",
    write_results_file = store,
    write_results_flush_every = 1L,
    return_output = FALSE
  )

  expected_gam <- ModelArray.gam(
    FD ~ age,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    var.smoothTerms = c(),
    var.parametricTerms = c("estimate", "std.error", "statistic", "p.value"),
    var.model = c("adj.r.squared", "dev.expl", "AIC"),
    correct.p.value.smoothTerms = c("none"),
    correct.p.value.parametricTerms = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )
  out_gam <- ModelArray.gam(
    FD ~ age,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    var.smoothTerms = c(),
    var.parametricTerms = c("estimate", "std.error", "statistic", "p.value"),
    var.model = c("adj.r.squared", "dev.expl", "AIC"),
    correct.p.value.smoothTerms = c("none"),
    correct.p.value.parametricTerms = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE,
    write_results_name = "gam_stream",
    write_results_file = store,
    write_results_flush_every = 1L,
    return_output = FALSE
  )

  wrap_fun <- function(data) {
    data.frame(
      mean_fd = mean(data$FD),
      sd_fd = stats::sd(data$FD)
    )
  }
  expected_wrap <- ModelArray.wrap(
    FUN = wrap_fun,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )
  out_wrap <- ModelArray.wrap(
    FUN = wrap_fun,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE,
    write_results_name = "wrap_stream",
    write_results_file = store,
    write_results_flush_every = 1L,
    return_output = FALSE
  )

  ma_with_results <- ModelArray(
    store,
    scalar_types = "FD",
    analysis_names = c("lm_stream", "gam_stream", "wrap_stream"),
    backend = "tiledb"
  )
  expected_outputs <- list(
    lm_stream = expected_lm,
    gam_stream = expected_gam,
    wrap_stream = expected_wrap
  )
  for (analysis_name in names(expected_outputs)) {
    saved <- as.matrix(results(ma_with_results, analysis_name)$results_matrix)
    expected <- as.matrix(expected_outputs[[analysis_name]])
    expect_identical(colnames(saved), colnames(expected))
    expect_identical(unname(saved), unname(expected))

    matrix_uri <- file.path(
      store,
      "results",
      analysis_name,
      "results_matrix"
    )
    fragment_info <- tiledb::tiledb_fragment_info(matrix_uri)
    expect_equal(
      tiledb::tiledb_fragment_info_get_num(fragment_info),
      length(element_subset)
    )
  }
  expect_null(out_lm)
  expect_null(out_gam)
  expect_null(out_wrap)
})

test_that("HDF5 and TileDB outputs are identical for regular and streamed writes", {
  skip_if_no_tiledb()

  src <- paste0("s", 1:8)
  values <- matrix(c(
    1, 3, 2, 7, 5, 9, 8, 6,
    2, 1, 4, 3, 6, 5, 7, 8,
    rep(NA_real_, 8)
  ), nrow = 3, byrow = TRUE)
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src),
    scalars = list(FD = values),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = seq(10, 80, by = 10))
  element_subset <- as.integer(c(3, 1, 2))

  regular <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    var.terms = c("estimate", "std.error", "statistic", "p.value"),
    var.model = c("adj.r.squared", "AIC"),
    correct.p.value.terms = c("none"),
    correct.p.value.model = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )

  h5_regular_file <- tempfile(fileext = ".h5")
  tiledb_regular_file <- tempfile(fileext = ".tdb")
  h5_stream_file <- tempfile(fileext = ".h5")
  tiledb_stream_file <- tempfile(fileext = ".tdb")
  output_paths <- c(
    h5_regular_file,
    tiledb_regular_file,
    h5_stream_file,
    tiledb_stream_file
  )
  on.exit(unlink(output_paths, recursive = TRUE), add = TRUE)

  writeResults(
    h5_regular_file,
    regular,
    analysis_name = "regular",
    backend = "hdf5"
  )
  writeResults(
    tiledb_regular_file,
    regular,
    analysis_name = "regular",
    backend = "tiledb"
  )

  stream_lm <- function(output_file) {
    ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phen,
      scalar = "FD",
      element.subset = element_subset,
      var.terms = c("estimate", "std.error", "statistic", "p.value"),
      var.model = c("adj.r.squared", "AIC"),
      correct.p.value.terms = c("none"),
      correct.p.value.model = c("none"),
      num.subj.lthr.abs = 0,
      n_cores = 1,
      pbar = FALSE,
      verbose = FALSE,
      write_results_name = "stream",
      write_results_file = output_file,
      write_results_flush_every = 1L,
      return_output = FALSE
    )
  }
  expect_null(stream_lm(h5_stream_file))
  expect_null(stream_lm(tiledb_stream_file))

  read_results_matrix <- function(output_file, analysis_name, backend) {
    saved <- as.matrix(DelayedArray::DelayedArray(ModelArraySeed(
      output_file,
      name = file.path("results", analysis_name, "results_matrix"),
      backend = backend
    )))
    column_names <- .read_result_column_names(
      output_file,
      analysis_name,
      backend
    )
    if (ncol(saved) != length(column_names)) {
      saved <- t(saved)
    }
    colnames(saved) <- column_names
    saved
  }

  h5_regular <- read_results_matrix(h5_regular_file, "regular", "hdf5")
  tiledb_regular <- read_results_matrix(
    tiledb_regular_file,
    "regular",
    "tiledb"
  )
  h5_stream <- read_results_matrix(h5_stream_file, "stream", "hdf5")
  tiledb_stream <- read_results_matrix(
    tiledb_stream_file,
    "stream",
    "tiledb"
  )

  expect_identical(h5_regular, tiledb_regular)
  expect_identical(h5_stream, tiledb_stream)
  expect_identical(h5_regular, h5_stream)
  expect_identical(tiledb_regular, tiledb_stream)
})

test_that("write_scalar streaming saves compressed scalars to TileDB", {
  skip_if_no_tiledb()

  src <- paste0("s", 1:8)
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src),
    scalars = list(FD = matrix(c(
      1, 3, 2, 7, 5, 9, 8, 6,
      2, 1, 4, 3, 6, 5, 7, 8
    ), nrow = 2, byrow = TRUE)),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = seq(10, 80, by = 10))
  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)

  out <- ModelArray.wrap(
    FUN = function(data) data.frame(mean_fd = mean(data$FD)),
    data = modelarray,
    phenotypes = phen,
    scalar = "FD",
    element.subset = as.integer(c(1, 2)),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE,
    write_scalar_name = "FD_harmonized",
    write_scalar_file = store,
    write_scalar_flush_every = 1L,
    write_scalar_compression_level = 6L,
    return_output = FALSE
  )

  saved_modelarray <- ModelArray(
    store,
    scalar_types = "FD_harmonized",
    backend = "tiledb"
  )
  saved <- as.matrix(scalars(saved_modelarray, "FD_harmonized"))
  values_uri <- file.path(store, "scalars", "FD_harmonized", "values")
  array <- tiledb::tiledb_array(values_uri, query_type = "READ", keep_open = TRUE)
  on.exit(try(tiledb::tiledb_array_close(array), silent = TRUE), add = TRUE)
  value_filters <- tiledb::filter_list(
    tiledb::attrs(tiledb::schema(array))[["values"]]
  )
  gzip_filter <- value_filters[0]
  fragment_info <- tiledb::tiledb_fragment_info(values_uri)

  expect_null(out)
  expect_equal(colnames(saved), "mean_fd")
  expect_equal(as.vector(saved), rowMeans(modelarray@scalars$FD))
  expect_equal(tiledb::tiledb_fragment_info_get_num(fragment_info), 2L)
  expect_equal(tiledb::nfilters(value_filters), 1L)
  expect_equal(tiledb::tiledb_filter_type(gzip_filter), "GZIP")
  expect_equal(
    tiledb::tiledb_filter_get_option(gzip_filter, "COMPRESSION_LEVEL"),
    6L
  )
})

test_that("elementMetadata reads TileDB top-level metadata arrays", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  dir.create(store, recursive = TRUE)
  metadata <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 3, byrow = TRUE)
  TileDBArray::writeTileDBArray(
    metadata,
    path = file.path(store, "fixels"),
    attr = "values"
  )
  ma <- methods::new(
    "ModelArray",
    sources = list(FD = c("s1", "s2")),
    scalars = list(FD = matrix(1, nrow = 3, ncol = 2)),
    results = list(),
    path = store
  )

  expect_equal(elementMetadata(ma), metadata)
})

test_that("elementMetadata returns NULL when no TileDB metadata arrays exist", {
  skip_if_no_tiledb()

  store <- tempfile(fileext = ".tdb")
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  dir.create(store, recursive = TRUE)
  ma <- methods::new(
    "ModelArray",
    sources = list(FD = c("s1", "s2")),
    scalars = list(FD = matrix(1, nrow = 2, ncol = 2)),
    results = list(),
    path = store
  )

  expect_null(elementMetadata(ma))
})
