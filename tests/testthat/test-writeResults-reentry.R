test_that("writeResults reentry mode writes scalar-style output (single-column)", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
  phenotypes <- read.csv(csv_path)

  h5_tmp <- tempfile(fileext = ".h5")
  file.copy(h5_path, h5_tmp, overwrite = TRUE)
  on.exit(unlink(h5_tmp), add = TRUE)

  modelarray <- ModelArray(h5_tmp, scalar_types = c("FD"))
  element.subset <- 1:10
  mylm <- ModelArray.lm(FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = "FD",
    element.subset = element.subset,
    n_cores = 1,
    pbar = FALSE
  )

  writeResults(
    h5_tmp,
    df.output = mylm,
    analysis_name = "result_lm",
    overwrite = TRUE,
    reentry = TRUE,
    reentry_scalar_name = "FD_reentry_single",
    reentry_col = "age.estimate",
    reentry_overwrite = TRUE
  )

  modelarray_new <- ModelArray(
    h5_tmp,
    scalar_types = c("FD", "FD_reentry_single"),
    analysis_names = c("result_lm")
  )

  fd_reentry <- scalars(modelarray_new)[["FD_reentry_single"]]
  expect_equal(nrow(fd_reentry), nrow(scalars(modelarray_new)[["FD"]]))
  expect_equal(ncol(fd_reentry), 1)

  idx <- mylm$element_id + 1L
  reentry_vals <- as.numeric(fd_reentry[idx, 1])
  expect_equal(reentry_vals, mylm$age.estimate, tolerance = 1e-12)

  # check one row outside subset remains NaN
  outside_idx <- max(idx) + 1L
  expect_true(is.nan(as.numeric(fd_reentry[outside_idx, 1])))
})

test_that("writeResults reentry mode writes scalar-style output (matrix-style columns)", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  h5_tmp <- tempfile(fileext = ".h5")
  file.copy(h5_path, h5_tmp, overwrite = TRUE)
  on.exit(unlink(h5_tmp), add = TRUE)

  modelarray <- ModelArray(h5_tmp, scalar_types = c("FD"))
  src_names <- sources(modelarray)[["FD"]]

  element.subset <- 1:10
  element_id <- element.subset - 1L
  set.seed(20260209)
  mock_mat <- matrix(
    rnorm(length(element.subset) * length(src_names)),
    nrow = length(element.subset),
    ncol = length(src_names)
  )
  colnames(mock_mat) <- src_names
  df_mock <- cbind(data.frame(element_id = element_id), as.data.frame(mock_mat))

  writeResults(
    h5_tmp,
    df.output = df_mock,
    analysis_name = "result_mock",
    overwrite = TRUE,
    reentry = TRUE,
    reentry_scalar_name = "FD_reentry_matrix",
    reentry_overwrite = TRUE
  )

  modelarray_new <- ModelArray(h5_tmp, scalar_types = c("FD_reentry_matrix"))
  fd_reentry <- scalars(modelarray_new)[["FD_reentry_matrix"]]
  expect_equal(dim(fd_reentry), c(nrow(scalars(modelarray)[["FD"]]), length(src_names)))
  expect_equal(colnames(fd_reentry), src_names)

  idx <- element_id + 1L
  loaded <- as.matrix(fd_reentry[idx, ])
  expect_equal(loaded, mock_mat, tolerance = 1e-12)
})

test_that("ModelArray scalar loader accepts legacy 'colnames' attribute", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  h5_tmp <- tempfile(fileext = ".h5")
  file.copy(h5_path, h5_tmp, overwrite = TRUE)
  on.exit(unlink(h5_tmp), add = TRUE)

  modelarray <- ModelArray(h5_tmp, scalar_types = c("FD"))
  n_elements <- nrow(scalars(modelarray)[["FD"]])
  legacy_vals <- matrix(NaN, nrow = n_elements, ncol = 1)
  legacy_vals[1:10, 1] <- seq_len(10)

  rhdf5::h5createGroup(h5_tmp, "scalars/FD_reentry_legacy")
  rhdf5::h5write(legacy_vals, h5_tmp, "scalars/FD_reentry_legacy/values")
  rhdf5::h5writeAttribute(
    attr = "legacy_source",
    h5obj = h5_tmp,
    h5loc = "/scalars/FD_reentry_legacy/values",
    name = "colnames"
  )

  legacy_loaded <- ModelArray(h5_tmp, scalar_types = c("FD_reentry_legacy"))
  expect_equal(dim(scalars(legacy_loaded)[["FD_reentry_legacy"]]), c(n_elements, 1))
  expect_equal(sources(legacy_loaded)[["FD_reentry_legacy"]], "legacy_source")
})

test_that("writeResults reentry rejects duplicate element_id values", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  h5_tmp <- tempfile(fileext = ".h5")
  file.copy(h5_path, h5_tmp, overwrite = TRUE)
  on.exit(unlink(h5_tmp), add = TRUE)

  df_bad <- data.frame(
    element_id = c(0, 0, 1),
    age.estimate = c(1.1, 2.2, 3.3)
  )

  expect_error(
    writeResults(
      h5_tmp,
      df.output = df_bad,
      analysis_name = "result_dup_element_id",
      overwrite = TRUE,
      reentry = TRUE,
      reentry_scalar_name = "FD_reentry_dup_element_id"
    ),
    "duplicate element_id"
  )
})

test_that("writeResults reentry rejects duplicate non-element_id column names", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  h5_tmp <- tempfile(fileext = ".h5")
  file.copy(h5_path, h5_tmp, overwrite = TRUE)
  on.exit(unlink(h5_tmp), add = TRUE)

  df_bad <- data.frame(
    element_id = c(0, 1, 2),
    x = c(10, 20, 30),
    x = c(40, 50, 60),
    check.names = FALSE
  )

  expect_error(
    writeResults(
      h5_tmp,
      df.output = df_bad,
      analysis_name = "result_dup_colnames",
      overwrite = TRUE,
      reentry = TRUE,
      reentry_scalar_name = "FD_reentry_dup_colnames"
    ),
    "duplicate column names"
  )
})

test_that("writeResults reentry uses explicit reentry_ref_scalar", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  h5_tmp <- tempfile(fileext = ".h5")
  file.copy(h5_path, h5_tmp, overwrite = TRUE)
  on.exit(unlink(h5_tmp), add = TRUE)

  # add a tiny scalar so default-first behavior could be wrong/non-deterministic
  tiny_vals <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  tiny_src <- c("tiny_a", "tiny_b")
  rhdf5::h5createGroup(h5_tmp, "scalars/AA_tiny")
  rhdf5::h5write(tiny_vals, h5_tmp, "scalars/AA_tiny/values")
  rhdf5::h5write(tiny_src, h5_tmp, "scalars/AA_tiny/column_names")
  rhdf5::h5writeAttribute(
    attr = tiny_src,
    h5obj = h5_tmp,
    h5loc = "/scalars/AA_tiny/values",
    name = "column_names"
  )

  modelarray <- ModelArray(h5_tmp, scalar_types = c("FD"))
  fd_src <- sources(modelarray)[["FD"]]
  element_id <- 0:4
  payload <- matrix(
    seq_len(length(element_id) * length(fd_src)),
    nrow = length(element_id),
    ncol = length(fd_src)
  )
  colnames(payload) <- fd_src
  df_ok <- cbind(data.frame(element_id = element_id), as.data.frame(payload))

  writeResults(
    h5_tmp,
    df.output = df_ok,
    analysis_name = "result_ref_scalar_ok",
    overwrite = TRUE,
    reentry = TRUE,
    reentry_scalar_name = "FD_reentry_ref_scalar",
    reentry_ref_scalar = "FD",
    reentry_overwrite = TRUE
  )

  ma_new <- ModelArray(h5_tmp, scalar_types = c("FD_reentry_ref_scalar"))
  expect_equal(colnames(scalars(ma_new)[["FD_reentry_ref_scalar"]]), fd_src)
  expect_equal(
    as.matrix(scalars(ma_new)[["FD_reentry_ref_scalar"]])[element_id + 1L, ],
    payload
  )

  expect_error(
    writeResults(
      h5_tmp,
      df.output = df_ok,
      analysis_name = "result_ref_scalar_bad",
      overwrite = TRUE,
      reentry = TRUE,
      reentry_scalar_name = "FD_reentry_ref_scalar_bad",
      reentry_ref_scalar = "DOES_NOT_EXIST",
      reentry_overwrite = TRUE
    ),
    "reentry_ref_scalar not found"
  )
})
