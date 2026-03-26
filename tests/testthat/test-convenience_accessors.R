test_that("nElements returns correct count", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  expect_equal(nElements(ma), 182581L)
  expect_equal(nElements(ma, "FD"), 182581L)
})

test_that("nInputFiles returns correct count", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  expect_equal(nInputFiles(ma), 50L)
  expect_equal(nInputFiles(ma, "FD"), 50L)
})

test_that("scalarNames returns scalar names", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  expect_equal(scalarNames(ma), "FD")
})

test_that("analysisNames returns empty for no analyses", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  expect_true(length(analysisNames(ma)) == 0)
})

test_that("elementMetadata returns fixel metadata", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  em <- elementMetadata(ma)
  expect_true(!is.null(em))
  expect_equal(nrow(em), 182581L)
})

test_that("h5summary returns correct structure", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  s <- h5summary(h5_path)

  expect_s3_class(s, "h5summary")
  expect_equal(s$scalars$name, "FD")
  expect_equal(s$scalars$nElements, 182581L)
  expect_equal(s$scalars$nInputFiles, 50L)
  expect_equal(s$filepath, h5_path)

  # Print method works
  output <- capture.output(print(s))
  expect_true(any(grepl("FD", output)))
  expect_true(any(grepl("182581", output)))
})

test_that("h5summary errors on missing file", {
  expect_error(h5summary("/nonexistent/file.h5"), "not found")
})
