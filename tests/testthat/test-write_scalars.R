test_that("writeScalars writes block-flushed scalars and supports overwrite controls", {
  h5_path <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_path), add = TRUE)

  mat1 <- matrix(c(
    1, 2,
    3, 4,
    5, 6
  ), nrow = 3, byrow = TRUE)
  sources1 <- c("sub1", "sub2")

  writeScalars(
    fn.output = h5_path,
    scalar_matrix = mat1,
    scalar_name = "HARM",
    column_names = sources1,
    overwrite = TRUE,
    flush_every = 2L
  )

  h5 <- hdf5r::H5File$new(h5_path, mode = "r")
  expect_true(h5$exists("scalars/HARM/values"))
  expect_true(h5$exists("scalars/HARM/column_names"))
  h5$close_all()

  vals1 <- rhdf5::h5read(h5_path, "scalars/HARM/values")
  attrs <- rhdf5::h5readAttributes(h5_path, "scalars/HARM/values")
  cn1 <- rhdf5::h5read(h5_path, "scalars/HARM/column_names")
  expect_equal(vals1, mat1)
  expect_equal(as.character(attrs$column_names), sources1)
  expect_equal(as.character(cn1), sources1)

  modelarray <- ModelArray(h5_path, scalar_types = c("HARM"))
  vals_loaded <- as.matrix(scalars(modelarray)[["HARM"]])
  dimnames(vals_loaded) <- NULL
  expect_equal(vals_loaded, mat1)
  expect_equal(sources(modelarray)[["HARM"]], sources1)

  mat2 <- matrix(c(
    10, 20,
    30, 40
  ), nrow = 2, byrow = TRUE)
  expect_warning(
    writeScalars(
      fn.output = h5_path,
      scalar_matrix = mat2,
      scalar_name = "HARM",
      column_names = sources1,
      overwrite = FALSE
    ),
    "exists but not to overwrite"
  )
  vals2 <- rhdf5::h5read(h5_path, "scalars/HARM/values")
  expect_equal(dim(vals2), dim(mat1))

  writeScalars(
    fn.output = h5_path,
    scalar_matrix = mat2,
    scalar_name = "HARM",
    column_names = sources1,
    overwrite = TRUE,
    flush_every = 1L
  )
  vals3 <- rhdf5::h5read(h5_path, "scalars/HARM/values")
  expect_equal(vals3, mat2)

  expect_error(
    writeScalars(
      fn.output = h5_path,
      scalar_matrix = mat2,
      scalar_name = "HARM_BAD",
      column_names = "sub1"
    ),
    "length\\(column_names\\) must equal ncol\\(scalar_matrix\\)"
  )
})
