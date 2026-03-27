test_that("ModelArray.lm can stream results to HDF5", {
  src <- paste0("s", 1:8)
  mat <- matrix(c(
    1, 3, 2, 7, 5, 9, 8, 6,
    2, 1, 4, 3, 6, 5, 7, 8,
    5, 2, 6, 1, 7, 3, 8, 4
  ), nrow = 3, byrow = TRUE)
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src),
    scalars = list(FD = mat),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = seq(10, 80, by = 10))
  h5_out <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_out), add = TRUE)

  out <- ModelArray.lm(
    formula = FD ~ age,
    data = modelarray,
    phenotypes = phen,
    scalar = "FD",
    element.subset = as.integer(c(1, 2, 3)),
    var.terms = c("estimate"),
    var.model = c("adj.r.squared"),
    correct.p.value.terms = c("none"),
    correct.p.value.model = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE,
    write_results_name = "lm_stream",
    write_results_file = h5_out,
    write_results_flush_every = 2L,
    return_output = TRUE
  )

  mat_h5 <- rhdf5::h5read(h5_out, "results/lm_stream/results_matrix")
  expect_equal(dim(mat_h5), dim(out))
  expect_equal(unname(mat_h5), unname(as.matrix(out)))
})

test_that("ModelArray.gam can stream results to HDF5", {
  src <- paste0("s", 1:8)
  mat <- matrix(c(
    1, 3, 2, 7, 5, 9, 8, 6,
    2, 1, 4, 3, 6, 5, 7, 8,
    5, 2, 6, 1, 7, 3, 8, 4
  ), nrow = 3, byrow = TRUE)
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src),
    scalars = list(FD = mat),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = seq(10, 80, by = 10))
  h5_out <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_out), add = TRUE)

  out <- ModelArray.gam(
    formula = FD ~ age,
    data = modelarray,
    phenotypes = phen,
    scalar = "FD",
    element.subset = as.integer(c(1, 2, 3)),
    var.smoothTerms = c(),
    var.parametricTerms = c("estimate"),
    var.model = c("dev.expl"),
    correct.p.value.smoothTerms = c("none"),
    correct.p.value.parametricTerms = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE,
    write_results_name = "gam_stream",
    write_results_file = h5_out,
    write_results_flush_every = 2L,
    return_output = TRUE
  )

  mat_h5 <- rhdf5::h5read(h5_out, "results/gam_stream/results_matrix")
  expect_equal(dim(mat_h5), dim(out))
  expect_equal(unname(mat_h5), unname(as.matrix(out)))
})

test_that("ModelArray.wrap can stream results to HDF5 without returning output", {
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

  fun_wrap <- function(data) {
    data.frame(m = mean(data$FD), a = data$age[1])
  }

  out <- ModelArray.wrap(
    FUN = fun_wrap,
    data = modelarray,
    phenotypes = phen,
    scalar = "FD",
    element.subset = as.integer(c(1, 2)),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE,
    write_results_name = "wrap_stream",
    write_results_file = h5_out,
    write_results_flush_every = 1L,
    return_output = FALSE
  )

  expect_null(out)
  mat_h5 <- rhdf5::h5read(h5_out, "results/wrap_stream/results_matrix")
  expect_equal(dim(mat_h5), c(2, 3))
})
