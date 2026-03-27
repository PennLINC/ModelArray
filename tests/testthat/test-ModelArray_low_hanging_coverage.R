test_that("S4 accessors and helpers cover low-hanging branches", {
  fd <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
  fa <- matrix(c(6, 5, 4, 3, 2, 1), nrow = 3, byrow = TRUE)
  src <- c("sub1", "sub2")

  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src, FA = src),
    scalars = list(FD = fd, FA = fa),
    results = list(my_analysis = list(results_matrix = matrix(1:3, ncol = 1))),
    path = tempfile(fileext = ".h5")
  )

  # accessors with and without optional selector argument
  expect_identical(scalars(modelarray), modelarray@scalars)
  expect_identical(scalars(modelarray, "FD"), fd)
  expect_identical(results(modelarray), modelarray@results)
  expect_identical(results(modelarray, "my_analysis"), modelarray@results$my_analysis)
  expect_identical(sources(modelarray), modelarray@sources)

  # show() formatting branch
  shown <- paste(capture.output(show(modelarray)), collapse = "\n")
  expect_match(shown, "elements x")
  expect_match(shown, "input files")
  expect_match(shown, "Analyses:")

  # helper success + error branches
  phenotypes <- data.frame(source_file = src, age = c(20, 30))
  out <- exampleElementData(modelarray, scalar = "FD", i_element = 2L, phenotypes = phenotypes)
  expect_true(is.data.frame(out))
  expect_identical(out$FD, as.numeric(fd[2, ]))

  expect_error(
    exampleElementData(modelarray, scalar = "FD", i_element = 1L, phenotypes = c(1, 2)),
    "phenotypes must be a data.frame"
  )
  expect_error(
    exampleElementData(modelarray, scalar = "NOT_A_SCALAR", i_element = 1L, phenotypes = phenotypes),
    "scalar not found in modelarray"
  )
  expect_error(
    exampleElementData(modelarray, scalar = "FD", i_element = 0L, phenotypes = phenotypes),
    "i_element is out of range"
  )
  expect_error(
    exampleElementData(modelarray, scalar = "FD", i_element = NA_integer_, phenotypes = phenotypes),
    "i_element is out of range"
  )

  expect_error(
    numElementsTotal(modelarray, scalar_name = "UNKNOWN"),
    "scalar_name requested in not in modelarray"
  )
})

test_that("ModelArray constructor falls back to dataset column names and transposes values", {
  h5_path <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_path), add = TRUE)

  h5 <- hdf5r::H5File$new(h5_path, mode = "w")
  scalars_grp <- h5$create_group("scalars")
  fd_grp <- scalars_grp$create_group("FD")
  # Store as subjects x elements so constructor needs to transpose.
  fd_grp[["values"]] <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  fd_grp[["column_names"]] <- c("subA", "subB")
  h5$close_all()

  modelarray <- ModelArray(
    h5_path,
    scalar_types = c("FD"),
    analysis_names = c("analysis_not_present")
  )

  expect_identical(sources(modelarray)$FD, c("subA", "subB"))
  expect_identical(dim(scalars(modelarray)$FD), c(3L, 2L))
  expect_identical(colnames(scalars(modelarray)$FD), c("subA", "subB"))
  # No /results group in file, so requesting analyses should return empty list.
  expect_identical(results(modelarray), list())
})

test_that("ModelArray constructor errors when scalar column names are unavailable", {
  h5_path <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_path), add = TRUE)

  h5 <- hdf5r::H5File$new(h5_path, mode = "w")
  scalars_grp <- h5$create_group("scalars")
  fd_grp <- scalars_grp$create_group("FD")
  fd_grp[["values"]] <- matrix(1:6, nrow = 3, ncol = 2)
  h5$close_all()

  expect_error(
    ModelArray(h5_path, scalar_types = c("FD")),
    "Neither attribute 'column_names' nor a dataset with column names found"
  )
})

test_that("ModelArray constructor supports backward-compatible column-name attributes", {
  h5_path <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_path), add = TRUE)

  rhdf5::h5createFile(h5_path)
  rhdf5::h5createGroup(h5_path, "scalars")
  rhdf5::h5createGroup(h5_path, "scalars/FD")
  rhdf5::h5write(matrix(1:6, nrow = 3, ncol = 2), h5_path, "scalars/FD/values")
  rhdf5::h5writeAttribute(
    attr = c("attr_sub1", "attr_sub2"),
    h5obj = h5_path,
    name = "column_names",
    h5loc = "scalars/FD/values"
  )

  rhdf5::h5createGroup(h5_path, "results")
  rhdf5::h5createGroup(h5_path, "results/my_analysis")
  rhdf5::h5write(matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2), h5_path, "results/my_analysis/results_matrix")
  rhdf5::h5writeAttribute(
    attr = c("beta", "p.value"),
    h5obj = h5_path,
    name = "colnames",
    h5loc = "results/my_analysis/results_matrix"
  )

  modelarray <- ModelArray(
    h5_path,
    scalar_types = c("FD"),
    analysis_names = c("my_analysis")
  )

  expect_identical(sources(modelarray)$FD, c("attr_sub1", "attr_sub2"))
  expect_identical(colnames(scalars(modelarray)$FD), c("attr_sub1", "attr_sub2"))
  expect_identical(
    colnames(results(modelarray)$my_analysis$results_matrix),
    c("beta", "p.value")
  )
})

test_that("ModelArray constructor prefers attributes over dataset column-name fallbacks", {
  h5_path <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_path), add = TRUE)

  rhdf5::h5createFile(h5_path)
  rhdf5::h5createGroup(h5_path, "scalars")
  rhdf5::h5createGroup(h5_path, "scalars/FD")
  rhdf5::h5write(matrix(1:6, nrow = 3, ncol = 2), h5_path, "scalars/FD/values")
  rhdf5::h5writeAttribute(
    attr = c("attr_sub1", "attr_sub2"),
    h5obj = h5_path,
    name = "column_names",
    h5loc = "scalars/FD/values"
  )
  # Conflicting fallback path should be ignored when attribute exists.
  rhdf5::h5write(c("fallback_sub1", "fallback_sub2"), h5_path, "scalars/FD/column_names")

  rhdf5::h5createGroup(h5_path, "results")
  rhdf5::h5createGroup(h5_path, "results/my_analysis")
  rhdf5::h5write(matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2), h5_path, "results/my_analysis/results_matrix")
  rhdf5::h5writeAttribute(
    attr = c("attr_beta", "attr_p.value"),
    h5obj = h5_path,
    name = "colnames",
    h5loc = "results/my_analysis/results_matrix"
  )
  # Conflicting fallback path should be ignored when attribute exists.
  rhdf5::h5write(c("fallback_beta", "fallback_p"), h5_path, "results/my_analysis/column_names")

  modelarray <- ModelArray(
    h5_path,
    scalar_types = c("FD"),
    analysis_names = c("my_analysis")
  )

  expect_identical(sources(modelarray)$FD, c("attr_sub1", "attr_sub2"))
  expect_identical(
    colnames(results(modelarray)$my_analysis$results_matrix),
    c("attr_beta", "attr_p.value")
  )
})

test_that("analyseOneElement.lm covers predictor/collision/error branches", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src, FA = src),
    scalars = list(
      FD = matrix(c(1, 2, 3, 4, 2, 3, 4, 5), nrow = 2, byrow = TRUE),
      FA = matrix(c(1, 2, NaN, 4, 5, 6, 7, 8), nrow = 2, byrow = TRUE)
    ),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))

  expect_error(
    analyseOneElement.lm(
      i_element = 1, formula = FD ~ age + FA, modelarray = modelarray,
      phenotypes = transform(phen, FA = 1:4), scalar = "FD",
      var.terms = c("estimate"), var.model = c("adj.r.squared"),
      num.subj.lthr = 0, flag_initiate = TRUE
    ),
    "Column name collision"
  )

  modelarray_bad_src <- methods::new(
    "ModelArray",
    sources = list(FD = src, FA = c("x1", "x2", "x3", "x4")),
    scalars = modelarray@scalars,
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  expect_error(
    analyseOneElement.lm(
      i_element = 1, formula = FD ~ age + FA, modelarray = modelarray_bad_src,
      phenotypes = phen, scalar = "FD",
      var.terms = c("estimate"), var.model = c("adj.r.squared"),
      num.subj.lthr = 0, flag_initiate = TRUE
    ),
    "do not match phenotypes\\$source_file"
  )

  # intersection threshold branch (finite FD but not enough finite FA)
  res_insuf <- analyseOneElement.lm(
    i_element = 1, formula = FD ~ age + FA, modelarray = modelarray,
    phenotypes = phen, scalar = "FD",
    var.terms = c("estimate"), var.model = c("adj.r.squared"),
    num.subj.lthr = 3, num.stat.output = 3, flag_initiate = FALSE
  )
  expect_identical(res_insuf[1], 0)
  expect_true(all(is.nan(res_insuf[-1])))

  expect_warning(
    res_skip <- analyseOneElement.lm(
      i_element = 1, formula = FD ~ missingVar, modelarray = modelarray,
      phenotypes = phen, scalar = "FD",
      var.terms = c("estimate"), var.model = c("adj.r.squared"),
      num.subj.lthr = 0, num.stat.output = 3, flag_initiate = FALSE,
      on_error = "skip"
    )
  )
  expect_identical(res_skip[1], 0)
  expect_true(all(is.nan(res_skip[-1])))

  expect_warning(
    init_skip <- analyseOneElement.lm(
      i_element = 1, formula = FD ~ missingVar, modelarray = modelarray,
      phenotypes = phen, scalar = "FD",
      var.terms = c("estimate"), var.model = c("adj.r.squared"),
      num.subj.lthr = 0, flag_initiate = TRUE, on_error = "skip"
    )
  )
  expect_true(is.nan(init_skip$column_names)[1])
})

test_that("analyseOneElement.gam covers predictor/collision/error branches", {
  src <- c("s1", "s2", "s3", "s4", "s5")
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src, FA = src),
    scalars = list(
      FD = matrix(c(1, 2, 3, 4, 5, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE),
      FA = matrix(c(1, 2, NaN, 4, 5, 5, 6, 7, 8, 9), nrow = 2, byrow = TRUE)
    ),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50, 60))

  expect_error(
    analyseOneElement.gam(
      i_element = 1, formula = FD ~ age + FA, modelarray = modelarray,
      phenotypes = transform(phen, FA = 1:5), scalar = "FD",
      var.smoothTerms = c("p.value"), var.parametricTerms = c("estimate"),
      var.model = c("dev.expl"),
      num.subj.lthr = 0, flag_initiate = TRUE
    ),
    "Column name collision"
  )

  modelarray_bad_src <- methods::new(
    "ModelArray",
    sources = list(FD = src, FA = c("x1", "x2", "x3", "x4", "x5")),
    scalars = modelarray@scalars,
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  expect_error(
    analyseOneElement.gam(
      i_element = 1, formula = FD ~ age + FA, modelarray = modelarray_bad_src,
      phenotypes = phen, scalar = "FD",
      var.smoothTerms = c("p.value"), var.parametricTerms = c("estimate"),
      var.model = c("dev.expl"),
      num.subj.lthr = 0, flag_initiate = TRUE
    ),
    "do not match phenotypes\\$source_file"
  )

  res_insuf <- analyseOneElement.gam(
    i_element = 1, formula = FD ~ age + FA, modelarray = modelarray,
    phenotypes = phen, scalar = "FD",
    var.smoothTerms = c("p.value"), var.parametricTerms = c("estimate"),
    var.model = c("dev.expl"),
    num.subj.lthr = 4, num.stat.output = 3, flag_initiate = FALSE
  )
  expect_identical(res_insuf[1], 0)
  expect_true(all(is.nan(res_insuf[-1])))

  expect_warning(
    res_skip <- analyseOneElement.gam(
      i_element = 1, formula = FD ~ missingVar, modelarray = modelarray,
      phenotypes = phen, scalar = "FD",
      var.smoothTerms = c("p.value"), var.parametricTerms = c("estimate"),
      var.model = c("dev.expl"),
      num.subj.lthr = 0, num.stat.output = 3, flag_initiate = FALSE,
      on_error = "skip"
    )
  )
  expect_identical(res_skip[1], 0)
  expect_true(all(is.nan(res_skip[-1])))
})

test_that("analyseOneElement.wrap covers coercion and error branches", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src, FA = src),
    scalars = list(
      FD = matrix(c(1, 2, 3, 4, 2, 3, 4, 5), nrow = 2, byrow = TRUE),
      FA = matrix(c(11, 12, 13, 14, 21, 22, 23, 24), nrow = 2, byrow = TRUE)
    ),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))

  list_fun <- function(data) list(m = mean(data$FD), a = mean(data$FA))
  init <- analyseOneElement.wrap(
    i_element = 1, user_fun = list_fun, modelarray = modelarray,
    phenotypes = phen, scalar = "FD",
    num.subj.lthr = 0, flag_initiate = TRUE
  )
  expect_true(all(c("element_id", "m", "a") %in% init$column_names))

  atomic_fun <- function(data) c(mean(data$FD), mean(data$FA))
  out_atomic <- analyseOneElement.wrap(
    i_element = 1, user_fun = atomic_fun, modelarray = modelarray,
    phenotypes = phen, scalar = "FD",
    num.subj.lthr = 0, num.stat.output = 3, flag_initiate = FALSE
  )
  expect_equal(length(out_atomic), 3)

  expect_error(
    analyseOneElement.wrap(
      i_element = 1, user_fun = function(data) data.frame(a = 1:2), modelarray = modelarray,
      phenotypes = phen, scalar = "FD",
      num.subj.lthr = 0, flag_initiate = FALSE, num.stat.output = 2
    ),
    "must return a one-row data.frame/tibble"
  )

  expect_error(
    analyseOneElement.wrap(
      i_element = 1, user_fun = function(data) new.env(), modelarray = modelarray,
      phenotypes = phen, scalar = "FD",
      num.subj.lthr = 0, flag_initiate = FALSE, num.stat.output = 2
    ),
    "Unsupported return type"
  )

  expect_warning(
    out_skip <- analyseOneElement.wrap(
      i_element = 1, user_fun = function(data) stop("boom"), modelarray = modelarray,
      phenotypes = phen, scalar = "FD",
      num.subj.lthr = 0, flag_initiate = FALSE, num.stat.output = 3, on_error = "skip"
    )
  )
  expect_identical(out_skip[1], 0)
  expect_true(all(is.nan(out_skip[-1])))

  expect_error(
    analyseOneElement.wrap(
      i_element = 1, user_fun = list_fun, modelarray = modelarray,
      phenotypes = transform(phen, FD = 1:4), scalar = "FD",
      num.subj.lthr = 0, flag_initiate = TRUE
    ),
    "Column name collision"
  )
})

test_that("writeResults handles validation, non-numeric LUT, and overwrite=false", {
  h5_path <- tempfile(fileext = ".h5")
  on.exit(unlink(h5_path), add = TRUE)

  expect_error(
    writeResults(h5_path, df.output = c(1, 2, 3), analysis_name = "bad"),
    "must be data of type `data.frame`"
  )

  df1 <- data.frame(
    element_id = 0:1,
    score = c(1.1, 2.2),
    group = c("A", "B")
  )
  writeResults(h5_path, df.output = df1, analysis_name = "my_res", overwrite = TRUE)

  h5 <- hdf5r::H5File$new(h5_path, mode = "r")
  expect_true(h5[["results/my_res"]]$exists("results_matrix"))
  expect_true(h5[["results/my_res"]]$exists("column_names"))
  expect_true(h5[["results/my_res"]]$exists("lut_forcol3"))
  mat1 <- h5[["results/my_res/results_matrix"]]$read()
  h5$close_all()

  df2 <- data.frame(element_id = 0:2, score = c(10, 20, 30))
  expect_warning(
    writeResults(h5_path, df.output = df2, analysis_name = "my_res", overwrite = FALSE),
    "exists but not to overwrite"
  )

  h5b <- hdf5r::H5File$new(h5_path, mode = "r")
  mat2 <- h5b[["results/my_res/results_matrix"]]$read()
  h5b$close_all()
  expect_identical(dim(mat2), dim(mat1))
})

test_that("ModelArray.wrap validation branches are exercised", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new(
    "ModelArray",
    sources = list(FD = src),
    scalars = list(FD = matrix(c(1, 2, 3, 4, 2, 3, 4, 5), nrow = 2, byrow = TRUE)),
    results = list(),
    path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))
  simple_fun <- function(data) data.frame(m = mean(data$FD))

  expect_error(
    ModelArray.wrap(simple_fun, data = list(), phenotypes = phen, scalar = "FD", element.subset = as.integer(1)),
    "data's class is not ModelArray"
  )
  expect_error(
    ModelArray.wrap(simple_fun, data = modelarray, phenotypes = phen, scalar = "FD", element.subset = c(1, 2)),
    "Please enter integers for element.subset"
  )
  expect_error(
    ModelArray.wrap(
      simple_fun, data = modelarray,
      phenotypes = subset(phen, select = -source_file),
      scalar = "FD", element.subset = as.integer(1)
    ),
    "Did not find column 'source_file'"
  )
  expect_error(
    ModelArray.wrap(
      simple_fun, data = modelarray,
      phenotypes = rbind(phen, phen[1, ]),
      scalar = "FD", element.subset = as.integer(1)
    ),
    "not the same as the length of the source file list in ModelArray"
  )

  # non-identical but matchable source_file order should be handled
  phen_swapped <- phen[c(2, 1, 3, 4), ]
  out <- ModelArray.wrap(
    simple_fun, data = modelarray, phenotypes = phen_swapped, scalar = "FD",
    element.subset = as.integer(c(1, 2)), n_cores = 1, pbar = FALSE, verbose = FALSE,
    num.subj.lthr.abs = 0
  )
  expect_equal(out$element_id, c(0, 1))
})
