test_that("ModelArray.wrap can stream selected columns to scalars", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src),
    scalars = list(FD = matrix(c(1, 2, 3, 4, 2, 3, 4, 5), nrow = 2, byrow = TRUE)),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))
  h5_out <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_out), add = TRUE)

  fun_harmonize <- function(data) {
    vals <- data$FD + 10
    names(vals) <- data$source_file
    vals
  }

  out <- ModelArray.wrap(
    FUN = fun_harmonize,
    data = modelarray,
    phenotypes = phen,
    scalar = "FD",
    element.subset = as.integer(c(1, 2)),
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE,
    num.subj.lthr.abs = 0,
    write_scalar_name = "FD_harmonized",
    write_scalar_file = h5_out,
    write_scalar_flush_every = 1L
  )

  expect_equal(dim(out), c(2, 5))
  expect_true(all(c("element_id", src) %in% colnames(out)))

  h5_vals <- rhdf5::h5read(h5_out, "scalars/FD_harmonized/values")
  expect_equal(
    h5_vals,
    matrix(c(11, 12, 13, 14, 12, 13, 14, 15), nrow = 2, byrow = TRUE)
  )

  ma_out <- ModelArray(h5_out, scalar_types = c("FD_harmonized"))
  loaded <- as.matrix(scalars(ma_out)[["FD_harmonized"]])
  dimnames(loaded) <- NULL
  expect_equal(loaded, h5_vals)
  expect_equal(sources(ma_out)[["FD_harmonized"]], src)
})

test_that("ModelArray.wrap write-scalar mode validates column metadata", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src),
    scalars = list(FD = matrix(c(1, 2, 3, 4, 2, 3, 4, 5), nrow = 2, byrow = TRUE)),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))
  h5_out <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_out), add = TRUE)

  simple_fun <- function(data) {
    vals <- data$FD
    names(vals) <- data$source_file
    vals
  }

  expect_error(
    ModelArray.wrap(
      FUN = simple_fun,
      data = modelarray,
      phenotypes = phen,
      scalar = "FD",
      element.subset = as.integer(c(1, 2)),
      n_cores = 1,
      pbar = FALSE,
      verbose = FALSE,
      num.subj.lthr.abs = 0,
      write_scalar_name = "FD_harmonized",
      write_scalar_file = h5_out,
      write_scalar_column_names = c("only_one_name")
    ),
    "length\\(write_scalar_column_names\\) must equal number of selected write_scalar_columns"
  )
})
