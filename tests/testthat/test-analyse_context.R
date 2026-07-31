# tests/testthat/test-analyse_context.R
#
# Tests for context builders (.build_*_context), shared element
# assembly (.assemble_element_data), per-element ctx vs. legacy
# equivalence, .find_initiator_element, and top-level integration.
#
# Organized to follow the two-tier hierarchy:
#   SECTION 1: Context builders (base → lm → gam → wrap)
#   SECTION 2: Shared per-element helper (.assemble_element_data)
#   SECTION 3: Per-element function ctx vs. legacy (lm → gam → wrap)
#   SECTION 4: Initiator search (.find_initiator_element)
#   SECTION 5: Top-level integration (ModelArray.gam, ModelArray.wrap)

# ---- Shared fixtures ----
# NOTE: Uses plain matrices in the scalars slot, not DelayedArray-backed
# HDF5 data. This works because the ModelArray S4 class accepts any
# matrix-like object for scalars [4]. Integration tests in Section 5
# use real HDF5-backed ModelArrays via system.file().

make_test_modelarray <- function(n_elements = 3, n_subjects = 4,
                                 scalars = c("FD"),
                                 inject_nan = FALSE) {
  src <- paste0("sub", seq_len(n_subjects))
  scalar_list <- list()
  source_list <- list()

  for (sn in scalars) {
    set.seed(match(sn, scalars))
    mat <- matrix(rnorm(n_elements * n_subjects),
                  nrow = n_elements, ncol = n_subjects)
    if (inject_nan && sn == scalars[1]) {
      mat[1, 1] <- NaN
      mat[1, 2] <- Inf
    }
    colnames(mat) <- src
    scalar_list[[sn]] <- mat
    source_list[[sn]] <- src
  }

  ma <- methods::new("ModelArray",
                     sources = source_list,
                     scalars = scalar_list,
                     results = list(),
                     path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 10, length.out = n_subjects),
    sex = rep(c("M", "F"), length.out = n_subjects),
    stringsAsFactors = FALSE
  )
  list(ma = ma, phen = phen)
}


# ==================================================================
# SECTION 1: Context builders (unit tests)
#   Ordered to match the two-tier hierarchy:
#   base → lm → gam → wrap
# ==================================================================

# ---- .build_base_context ----

test_that(".build_base_context attaches all scalars when scalar_subset is NULL", {
  fix <- make_test_modelarray(scalars = c("FD", "FC"))
  ctx <- .build_base_context(fix$ma, fix$phen, scalar = "FD",
                             scalar_subset = NULL)

  expect_identical(ctx$scalar, "FD")
  expect_identical(sort(ctx$attached_scalars), sort(c("FD", "FC")))
  expect_identical(sort(ctx$all_scalar_names), sort(c("FD", "FC")))
})

test_that(".build_base_context attaches only requested scalars when scalar_subset is provided", {
  fix <- make_test_modelarray(scalars = c("FD", "FC"))
  ctx <- .build_base_context(fix$ma, fix$phen, scalar = "FD",
                             scalar_subset = c("FD"))

  expect_identical(ctx$attached_scalars, "FD")
  # FC is in the ModelArray but not attached
  expect_true("FC" %in% ctx$all_scalar_names)
  expect_false("FC" %in% ctx$attached_scalars)
})

test_that(".build_base_context detects collision between scalar names and phenotype columns", {
  fix <- make_test_modelarray(scalars = c("FD"))
  phen_collision <- fix$phen
  phen_collision$FD <- seq_len(nrow(phen_collision))

  expect_error(
    .build_base_context(fix$ma, phen_collision, scalar = "FD",
                        scalar_subset = c("FD")),
    "Column name collision"
  )
})

test_that(".build_base_context sets predictor_reorder to NULL for aligned sources", {
  fix <- make_test_modelarray(scalars = c("FD"))
  ctx <- .build_base_context(fix$ma, fix$phen, scalar = "FD",
                             scalar_subset = c("FD"))

  expect_null(ctx$predictor_reorder[["FD"]])
})

test_that(".build_base_context computes reorder indices for misaligned sources", {
  fix <- make_test_modelarray(scalars = c("FD", "FC"))

  # Manually scramble FC sources
  fc_src <- rev(sources(fix$ma)[["FC"]])
  fix$ma@sources[["FC"]] <- fc_src

  ctx <- .build_base_context(fix$ma, fix$phen, scalar = "FD",
                             scalar_subset = c("FD", "FC"))

  expect_null(ctx$predictor_reorder[["FD"]])
  expect_false(is.null(ctx$predictor_reorder[["FC"]]))
  # Applying the reorder should recover the phenotype source order
  expect_identical(
    fc_src[ctx$predictor_reorder[["FC"]]],
    fix$phen$source_file
  )
})

test_that(".build_base_context errors when scalar sources don't match phenotypes", {
  fix <- make_test_modelarray(scalars = c("FD", "FC"))
  fix$ma@sources[["FC"]] <- paste0("wrong", seq_along(sources(fix$ma)[["FC"]]))

  expect_error(
    .build_base_context(fix$ma, fix$phen, scalar = "FD",
                        scalar_subset = c("FD", "FC")),
    "do not match phenotypes"
  )
})

# ---- .build_lm_context ----

test_that(".build_lm_context extracts formula components correctly", {
  fix <- make_test_modelarray(scalars = c("FD"))
  ctx <- .build_lm_context(FD ~ age + sex, fix$ma, fix$phen, scalar = "FD")

  expect_identical(ctx$formula, FD ~ age + sex)
  expect_identical(ctx$lhs_name, "FD")
  expect_true("age" %in% ctx$rhs_vars)
  expect_true("sex" %in% ctx$rhs_vars)
  expect_identical(ctx$scalar_predictors, character(0))
})

test_that(".build_lm_context detects scalar predictors in cross-scalar formulas", {
  fix <- make_test_modelarray(scalars = c("FD", "FC"))
  ctx <- .build_lm_context(FD ~ age + FC, fix$ma, fix$phen, scalar = "FD")

  expect_identical(ctx$scalar_predictors, "FC")
  expect_true("FC" %in% ctx$attached_scalars)
  expect_true("FD" %in% ctx$attached_scalars)
})

test_that(".build_lm_context does not attach unrelated scalars", {
  fix <- make_test_modelarray(scalars = c("FD", "FC", "FA"))
  ctx <- .build_lm_context(FD ~ age + FC, fix$ma, fix$phen, scalar = "FD")

  expect_true("FA" %in% ctx$all_scalar_names)
  expect_false("FA" %in% ctx$attached_scalars)
})

# ---- .build_gam_context ----

test_that(".build_gam_context caches interpret.gam breakdown", {
  fix <- make_test_modelarray(scalars = c("FD"))
  ctx <- .build_gam_context(FD ~ s(age) + sex, fix$ma, fix$phen, scalar = "FD")

  expect_true(!is.null(ctx$gam_formula_breakdown))
  expect_true("smooth.spec" %in% names(ctx$gam_formula_breakdown))
})

test_that(".build_gam_context errors on invalid GAM formula", {
  fix <- make_test_modelarray(scalars = c("FD"))

  expect_error(
    .build_gam_context(FD ~ s(age, d = 999), fix$ma, fix$phen, scalar = "FD"),
    "not valid for mgcv::gam"
  )
})

test_that(".build_gam_context inherits lm context fields", {
  fix <- make_test_modelarray(scalars = c("FD"))
  ctx <- .build_gam_context(FD ~ s(age) + sex, fix$ma, fix$phen, scalar = "FD")

  # Should have all the lm context fields
  expect_identical(ctx$lhs_name, "FD")
  expect_true("age" %in% ctx$rhs_vars)
  expect_true("sex" %in% ctx$rhs_vars)
  expect_identical(ctx$scalar_predictors, character(0))
})

test_that(".build_gam_context detects scalar predictors same as lm context", {
  fix <- make_test_modelarray(scalars = c("FD", "FC"))
  ctx <- .build_gam_context(FD ~ s(age) + FC, fix$ma, fix$phen, scalar = "FD")

  expect_identical(ctx$scalar_predictors, "FC")
  expect_true("FC" %in% ctx$attached_scalars)
  expect_true("FD" %in% ctx$attached_scalars)
})

# ---- .build_wrap_context ----

test_that(".build_wrap_context attaches all scalars", {
  fix <- make_test_modelarray(scalars = c("FD", "FC"))
  ctx <- .build_wrap_context(fix$ma, fix$phen, scalar = "FD")

  expect_identical(sort(ctx$attached_scalars), sort(c("FD", "FC")))
})

test_that(".build_wrap_context has no formula fields", {
  fix <- make_test_modelarray(scalars = c("FD"))
  ctx <- .build_wrap_context(fix$ma, fix$phen, scalar = "FD")

  expect_null(ctx$formula)
  expect_null(ctx$lhs_name)
})


# ==================================================================
# SECTION 2: Shared per-element helper
# ==================================================================

# ---- .assemble_element_data ----

test_that(".assemble_element_data returns filtered data.frame for valid element", {
  fix <- make_test_modelarray(scalars = c("FD"))
  ctx <- .build_lm_context(FD ~ age, fix$ma, fix$phen, scalar = "FD")

  elem <- .assemble_element_data(1, ctx, num.subj.lthr = 0)

  expect_true(elem$sufficient)
  expect_true(is.data.frame(elem$dat))
  expect_true("FD" %in% colnames(elem$dat))
  expect_true("age" %in% colnames(elem$dat))
  expect_true("source_file" %in% colnames(elem$dat))
})

test_that(".assemble_element_data returns insufficient when below threshold", {
  fix <- make_test_modelarray(n_subjects = 4, scalars = c("FD"))
  ctx <- .build_lm_context(FD ~ age, fix$ma, fix$phen, scalar = "FD")

  elem <- .assemble_element_data(1, ctx, num.subj.lthr = 100)

  expect_false(elem$sufficient)
  expect_null(elem$dat)
})

test_that(".assemble_element_data excludes NaN/Inf subjects from mask", {
  fix <- make_test_modelarray(n_subjects = 4, scalars = c("FD"),
                              inject_nan = TRUE)
  ctx <- .build_lm_context(FD ~ age, fix$ma, fix$phen, scalar = "FD")

  elem <- .assemble_element_data(1, ctx, num.subj.lthr = 0)

  # Element 1 has NaN in subject 1 and Inf in subject 2
  expect_true(elem$sufficient)
  expect_equal(elem$num_valid, 2)  # only subjects 3 and 4
  expect_equal(nrow(elem$dat), 2)
})

test_that(".assemble_element_data applies precomputed reorder for cross-scalar", {
  fix <- make_test_modelarray(n_elements = 2, n_subjects = 4,
                              scalars = c("FD", "FC"))

  # Scramble FC sources
  fc_src_scrambled <- rev(sources(fix$ma)[["FC"]])
  fix$ma@sources[["FC"]] <- fc_src_scrambled

  ctx <- .build_lm_context(FD ~ age + FC, fix$ma, fix$phen, scalar = "FD")
  elem <- .assemble_element_data(1, ctx, num.subj.lthr = 0)

  expect_true(elem$sufficient)
  expect_true("FC" %in% colnames(elem$dat))

  # The FC values should be reordered to match phenotype order
  fc_raw <- scalars(fix$ma)[["FC"]][1, ]
  reorder_idx <- ctx$predictor_reorder[["FC"]]
  fc_reordered <- fc_raw[reorder_idx]
  expect_equal(as.numeric(elem$dat$FC), as.numeric(fc_reordered))
})

test_that(".assemble_element_data intersection mask works across scalars", {
  n_subj <- 6
  src <- paste0("sub", 1:n_subj)

  fd_mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 1)
  fc_mat <- matrix(c(10, NaN, 30, NaN, 50, 60), nrow = 1)
  colnames(fd_mat) <- src
  colnames(fc_mat) <- src

  ma <- methods::new("ModelArray",
                     sources = list(FD = src, FC = src),
                     scalars = list(FD = fd_mat, FC = fc_mat),
                     results = list(),
                     path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 10, length.out = n_subj)
  )

  ctx <- .build_lm_context(FD ~ age + FC, ma, phen, scalar = "FD")
  elem <- .assemble_element_data(1, ctx, num.subj.lthr = 0)

  # Subjects 2 and 4 have NaN in FC, so should be excluded
  expect_equal(elem$num_valid, 4)
  expect_equal(nrow(elem$dat), 4)
  expect_false(any(is.nan(elem$dat$FC)))
})

test_that(".with_scalar_row_cache serves cached rows for TileDB-style stores", {
  store <- tempfile(fileext = ".tdb")
  dir.create(store)
  on.exit(unlink(store, recursive = TRUE), add = TRUE)

  old_block_size <- getOption("ModelArray.tiledb_read_block_size")
  on.exit(options(ModelArray.tiledb_read_block_size = old_block_size), add = TRUE)
  options(ModelArray.tiledb_read_block_size = 2L)

  src <- paste0("sub", 1:4)
  values <- matrix(seq_len(12), nrow = 3, ncol = 4)
  colnames(values) <- src
  ma <- methods::new("ModelArray",
                     sources = list(FD = src),
                     scalars = list(FD = values),
                     results = list(),
                     path = store
  )
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 10, length.out = length(src))
  )

  ctx <- .build_lm_context(FD ~ age, ma, phen, scalar = "FD")
  cached_ctx <- .with_scalar_row_cache(ctx, as.integer(c(3, 1)))

  expect_named(cached_ctx$scalar_row_cache, "FD")
  expect_equal(cached_ctx$scalar_row_cache[["FD"]], values[c(3, 1), , drop = FALSE])

  cached_ctx$modelarray@scalars[["FD"]][, ] <- 999
  elem <- .assemble_element_data(3L, cached_ctx, num.subj.lthr = 0)

  expect_true(elem$sufficient)
  expect_equal(as.numeric(elem$dat$FD), as.numeric(values[3, ]))
})

test_that("TileDB row caching honors the backend retained by ModelArray", {
  path_with_tiledb_suffix <- tempfile(fileext = ".tdb")
  file.create(path_with_tiledb_suffix)
  on.exit(unlink(path_with_tiledb_suffix), add = TRUE)

  src <- c("sub1", "sub2")
  values <- matrix(seq_len(4), nrow = 2, ncol = 2)
  colnames(values) <- src
  ma <- methods::new(
    "ModelArray",
    sources = list(FD = src),
    scalars = list(FD = values),
    results = list(),
    path = path_with_tiledb_suffix,
    backend = "hdf5"
  )
  phen <- data.frame(source_file = src, age = c(20, 30))
  ctx <- .build_lm_context(FD ~ age, ma, phen, scalar = "FD")

  expect_identical(.modelarray_backend(ma, path_with_tiledb_suffix), "hdf5")
  expect_false(.has_tiledb_attached_scalar(ctx))
  expect_null(.scalar_read_block_size(ctx))
})

test_that("TileDB read block sizing is memory bounded without a hard row cap", {
  store <- tempfile(fileext = ".tdb")
  dir.create(store)
  on.exit(unlink(store, recursive = TRUE), add = TRUE)

  old_block_size <- getOption("ModelArray.tiledb_read_block_size")
  old_block_mb <- getOption("ModelArray.tiledb_read_block_mb")
  on.exit(options(
    ModelArray.tiledb_read_block_size = old_block_size,
    ModelArray.tiledb_read_block_mb = old_block_mb
  ), add = TRUE)
  options(
    ModelArray.tiledb_read_block_size = NULL,
    ModelArray.tiledb_read_block_mb = 1
  )

  src <- paste0("sub", 1:2)
  values <- matrix(seq_len(4), nrow = 2, ncol = 2)
  colnames(values) <- src
  ma <- methods::new("ModelArray",
                     sources = list(FD = src),
                     scalars = list(FD = values),
                     results = list(),
                     path = store
  )
  phen <- data.frame(source_file = src, age = c(20, 30))
  ctx <- .build_lm_context(FD ~ age, ma, phen, scalar = "FD")

  expect_equal(.scalar_read_block_size(ctx), 65536L)

  options(ModelArray.tiledb_read_block_mb = -1)
  expect_error(
    .scalar_read_block_size(ctx),
    "ModelArray.tiledb_read_block_mb must be a positive single number"
  )
})

test_that(".report_iteration_progress labels chunked progress bars", {
  out <- capture.output(
    .report_iteration_progress(
      iteration = 1L,
      total_iterations = 3L,
      chunk_start = 1L,
      chunk_end = 10L,
      total_elements = 25L,
      pbar = TRUE
    ),
    type = "message"
  )
  expect_match(paste(out, collapse = "\n"), "iteration 1/3")
  expect_match(paste(out, collapse = "\n"), "elements 1-10 of 25")

  quiet <- capture.output(
    .report_iteration_progress(1L, 1L, 1L, 10L, 10L, TRUE),
    type = "message"
  )
  expect_length(quiet, 0L)
})

test_that("ModelArray.lm returns same output with TileDB-style row cache enabled", {
  store <- tempfile(fileext = ".tdb")
  dir.create(store)
  on.exit(unlink(store, recursive = TRUE), add = TRUE)

  old_block_size <- getOption("ModelArray.tiledb_read_block_size")
  on.exit(options(ModelArray.tiledb_read_block_size = old_block_size), add = TRUE)

  src <- paste0("sub", 1:8)
  values <- matrix(c(
    1, 3, 2, 7, 5, 9, 8, 6,
    2, 1, 4, 3, 6, 5, 7, 8,
    3, 2, 5, 4, 7, 6, 8, 9,
    4, 6, 5, 8, 9, 7, 10, 11
  ), nrow = 4, byrow = TRUE)
  colnames(values) <- src
  ma <- methods::new("ModelArray",
                     sources = list(FD = src),
                     scalars = list(FD = values),
                     results = list(),
                     path = store
  )
  phen <- data.frame(source_file = src, age = seq(10, 80, by = 10))
  element_subset <- as.integer(c(4, 1, 3, 2))

  options(ModelArray.tiledb_read_block_size = 2L)
  cached <- ModelArray.lm(
    FD ~ age,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    var.terms = c("estimate"),
    var.model = c("adj.r.squared"),
    correct.p.value.terms = c("none"),
    correct.p.value.model = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )

  options(ModelArray.tiledb_read_block_size = 0L)
  uncached <- ModelArray.lm(
    FD ~ age,
    data = ma,
    phenotypes = phen,
    scalar = "FD",
    element.subset = element_subset,
    var.terms = c("estimate"),
    var.model = c("adj.r.squared"),
    correct.p.value.terms = c("none"),
    correct.p.value.model = c("none"),
    num.subj.lthr.abs = 0,
    n_cores = 1,
    pbar = FALSE,
    verbose = FALSE
  )

  expect_equal(cached, uncached)
})


# ==================================================================
# SECTION 3: Per-element function ctx vs. legacy equivalence
#   Ordered: lm → gam → wrap
# ==================================================================

# ---- analyseOneElement.lm ----

test_that("analyseOneElement.lm produces identical results with ctx vs. without", {
  # Create fixture with 8 subjects (enough for a meaningful lm fit)
  fix <- make_test_modelarray(n_elements = 2, n_subjects = 8, scalars = c("FD"))

  # The fixture's phen has age and sex; rebuild with only age for a
  # clean single-predictor test. The fixture already has 8 subjects.
  src <- sources(fix$ma)[["FD"]]
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 5, length.out = length(src)),
    stringsAsFactors = FALSE
  )

  ctx <- .build_lm_context(FD ~ age, fix$ma, phen, scalar = "FD")

  # With context — initiation
  result_ctx <- analyseOneElement.lm(
    i_element = 1,
    ctx = ctx,
    var.terms = c("estimate", "p.value"),
    var.model = c("adj.r.squared"),
    num.subj.lthr = 0,
    num.stat.output = NULL,
    flag_initiate = TRUE
  )

  # Without context — legacy path
  result_legacy <- analyseOneElement.lm(
    i_element = 1,
    formula = FD ~ age,
    modelarray = fix$ma,
    phenotypes = phen,
    scalar = "FD",
    var.terms = c("estimate", "p.value"),
    var.model = c("adj.r.squared"),
    num.subj.lthr = 0,
    num.stat.output = NULL,
    flag_initiate = TRUE
  )

  expect_identical(result_ctx$column_names, result_legacy$column_names)
  expect_identical(result_ctx$list.terms, result_legacy$list.terms)

  # Also check runtime (non-initiate) path
  num_stat <- length(result_ctx$column_names)
  result_ctx_run <- analyseOneElement.lm(
    i_element = 1,
    ctx = ctx,
    var.terms = c("estimate", "p.value"),
    var.model = c("adj.r.squared"),
    num.subj.lthr = 0,
    num.stat.output = num_stat,
    flag_initiate = FALSE
  )
  result_legacy_run <- analyseOneElement.lm(
    i_element = 1,
    formula = FD ~ age,
    modelarray = fix$ma,
    phenotypes = phen,
    scalar = "FD",
    var.terms = c("estimate", "p.value"),
    var.model = c("adj.r.squared"),
    num.subj.lthr = 0,
    num.stat.output = num_stat,
    flag_initiate = FALSE
  )

  expect_equal(result_ctx_run, result_legacy_run)
})

# ---- analyseOneElement.gam ----

test_that("analyseOneElement.gam produces identical results with ctx vs. without", {
  src <- paste0("s", 1:8)
  set.seed(42)
  mat <- matrix(rnorm(3 * 8), nrow = 3)
  colnames(mat) <- src

  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src),
                             scalars = list(FD = mat),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 5, length.out = 8),
    stringsAsFactors = FALSE
  )

  ctx <- .build_gam_context(FD ~ s(age, k = 4), modelarray, phen, scalar = "FD")

  # With context — initiation
  result_ctx <- analyseOneElement.gam(
    i_element = 1,
    ctx = ctx,
    var.smoothTerms = c("statistic", "p.value"),
    var.parametricTerms = c("estimate", "p.value"),
    var.model = c("dev.expl"),
    num.subj.lthr = 0,
    num.stat.output = NULL,
    flag_initiate = TRUE
  )

  # Without context — legacy path
  result_legacy <- analyseOneElement.gam(
    i_element = 1,
    formula = FD ~ s(age, k = 4),
    modelarray = modelarray,
    phenotypes = phen,
    scalar = "FD",
    var.smoothTerms = c("statistic", "p.value"),
    var.parametricTerms = c("estimate", "p.value"),
    var.model = c("dev.expl"),
    num.subj.lthr = 0,
    num.stat.output = NULL,
    flag_initiate = TRUE
  )

  expect_identical(result_ctx$column_names, result_legacy$column_names)
  expect_identical(result_ctx$list.smoothTerms, result_legacy$list.smoothTerms)
  expect_identical(result_ctx$list.parametricTerms, result_legacy$list.parametricTerms)

  # Runtime path
  num_stat <- length(result_ctx$column_names)

  result_ctx_run <- analyseOneElement.gam(
    i_element = 1,
    ctx = ctx,
    var.smoothTerms = c("statistic", "p.value"),
    var.parametricTerms = c("estimate", "p.value"),
    var.model = c("dev.expl"),
    num.subj.lthr = 0,
    num.stat.output = num_stat,
    flag_initiate = FALSE
  )

  result_legacy_run <- analyseOneElement.gam(
    i_element = 1,
    formula = FD ~ s(age, k = 4),
    modelarray = modelarray,
    phenotypes = phen,
    scalar = "FD",
    var.smoothTerms = c("statistic", "p.value"),
    var.parametricTerms = c("estimate", "p.value"),
    var.model = c("dev.expl"),
    num.subj.lthr = 0,
    num.stat.output = num_stat,
    flag_initiate = FALSE
  )

  expect_equal(result_ctx_run, result_legacy_run)
})

test_that("analyseOneElement.gam with ctx normalizes (Intercept) to Intercept", {
  src <- paste0("s", 1:8)
  set.seed(42)
  mat <- matrix(rnorm(2 * 8), nrow = 2)
  colnames(mat) <- src

  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src),
                             scalars = list(FD = mat),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 5, length.out = 8),
    stringsAsFactors = FALSE
  )

  ctx <- .build_gam_context(FD ~ s(age, k = 4) + 1, modelarray, phen, scalar = "FD")

  result <- analyseOneElement.gam(
    i_element = 1,
    ctx = ctx,
    var.smoothTerms = c("p.value"),
    var.parametricTerms = c("estimate"),
    var.model = c("dev.expl"),
    num.subj.lthr = 0,
    num.stat.output = NULL,
    flag_initiate = TRUE
  )

  # Should contain "Intercept", not "(Intercept)"
  param_cols <- grep("estimate", result$column_names, value = TRUE)
  expect_true(any(grepl("^Intercept\\.", param_cols)))
  expect_false(any(grepl("\\(Intercept\\)", result$column_names)))
})

test_that("analyseOneElement.gam with ctx returns NaN row for insufficient subjects", {
  src <- paste0("s", 1:4)
  mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2)
  colnames(mat) <- src

  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src),
                             scalars = list(FD = mat),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(
    source_file = src,
    age = c(20, 30, 40, 50)
  )

  ctx <- .build_gam_context(FD ~ age, modelarray, phen, scalar = "FD")

  # Initiation with impossible threshold
  result_init <- analyseOneElement.gam(
    i_element = 1,
    ctx = ctx,
    var.smoothTerms = c(),
    var.parametricTerms = c("estimate"),
    var.model = c("dev.expl"),
    num.subj.lthr = 1000,
    flag_initiate = TRUE
  )
  expect_true(is.nan(result_init$column_names[1]))

  # Runtime with impossible threshold
  result_run <- analyseOneElement.gam(
    i_element = 1,
    ctx = ctx,
    var.smoothTerms = c(),
    var.parametricTerms = c("estimate"),
    var.model = c("dev.expl"),
    num.subj.lthr = 1000,
    num.stat.output = 5,
    flag_initiate = FALSE
  )
  expect_equal(result_run[1], 0)
  expect_true(all(is.nan(result_run[-1])))
})

test_that("analyseOneElement.gam with ctx handles on_error = 'skip'", {
  src <- paste0("s", 1:8)
  set.seed(42)
  mat <- matrix(rnorm(2 * 8), nrow = 2)
  colnames(mat) <- src

  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src),
                             scalars = list(FD = mat),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 5, length.out = 8)
  )

  ctx <- .build_gam_context(FD ~ nonexistent_var, modelarray, phen, scalar = "FD")

  expect_warning(
    result <- analyseOneElement.gam(
      i_element = 1,
      ctx = ctx,
      var.smoothTerms = c(),
      var.parametricTerms = c("estimate"),
      var.model = c("dev.expl"),
      num.subj.lthr = 0,
      num.stat.output = 4,
      flag_initiate = FALSE,
      on_error = "skip"
    )
  )
  expect_equal(result[1], 0)
  expect_true(all(is.nan(result[-1])))
})

# ---- analyseOneElement.wrap ----

test_that("analyseOneElement.wrap produces identical results with ctx vs. without", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src, FA = src),
                             scalars = list(
                               FD = matrix(c(1, 2, 3, 4, 2, 3, 4, 5), nrow = 2, byrow = TRUE),
                               FA = matrix(c(10, 20, 30, 40, 50, 60, 70, 80), nrow = 2, byrow = TRUE)
                             ),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(
    source_file = src,
    age = c(20, 30, 40, 50)
  )

  my_fun <- function(data, ...) {
    list(fd_mean = mean(data$FD), fa_mean = mean(data$FA))
  }

  ctx <- .build_wrap_context(modelarray, phen, scalar = "FD")

  # With context — initiation
  result_ctx <- analyseOneElement.wrap(
    i_element = 1,
    user_fun = my_fun,
    ctx = ctx,
    num.subj.lthr = 0,
    flag_initiate = TRUE
  )

  # Without context — legacy
  result_legacy <- analyseOneElement.wrap(
    i_element = 1,
    user_fun = my_fun,
    modelarray = modelarray,
    phenotypes = phen,
    scalar = "FD",
    num.subj.lthr = 0,
    flag_initiate = TRUE
  )

  expect_identical(result_ctx$column_names, result_legacy$column_names)

  # Runtime
  num_stat <- length(result_ctx$column_names)

  result_ctx_run <- analyseOneElement.wrap(
    i_element = 1,
    user_fun = my_fun,
    ctx = ctx,
    num.subj.lthr = 0,
    num.stat.output = num_stat,
    flag_initiate = FALSE
  )

  result_legacy_run <- analyseOneElement.wrap(
    i_element = 1,
    user_fun = my_fun,
    modelarray = modelarray,
    phenotypes = phen,
    scalar = "FD",
    num.subj.lthr = 0,
    num.stat.output = num_stat,
    flag_initiate = FALSE
  )

  expect_equal(result_ctx_run, result_legacy_run)
})

test_that("analyseOneElement.wrap with ctx attaches ALL scalars", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src, FA = src, FC = src),
                             scalars = list(
                               FD = matrix(1:8, nrow = 2),
                               FA = matrix(11:18, nrow = 2),
                               FC = matrix(21:28, nrow = 2)
                             ),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))

  check_fun <- function(data, ...) {
    # Verify all three scalars are present
    list(
      has_FD = as.numeric("FD" %in% colnames(data)),
      has_FA = as.numeric("FA" %in% colnames(data)),
      has_FC = as.numeric("FC" %in% colnames(data))
    )
  }

  ctx <- .build_wrap_context(modelarray, phen, scalar = "FD")

  result <- analyseOneElement.wrap(
    i_element = 1,
    user_fun = check_fun,
    ctx = ctx,
    num.subj.lthr = 0,
    num.stat.output = 4,
    flag_initiate = FALSE
  )

  # element_id = 0, then has_FD=1, has_FA=1, has_FC=1
  expect_equal(result[2], 1)
  expect_equal(result[3], 1)
  expect_equal(result[4], 1)
})

test_that("analyseOneElement.wrap with ctx returns NaN for insufficient subjects", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src),
                             scalars = list(FD = matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, byrow = TRUE)),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))
  my_fun <- function(data, ...) list(m = mean(data$FD))

  ctx <- .build_wrap_context(modelarray, phen, scalar = "FD")

  result_init <- analyseOneElement.wrap(
    i_element = 1,
    user_fun = my_fun,
    ctx = ctx,
    num.subj.lthr = 1000,
    flag_initiate = TRUE
  )
  expect_true(is.nan(result_init$column_names[1]))

  result_run <- analyseOneElement.wrap(
    i_element = 1,
    user_fun = my_fun,
    ctx = ctx,
    num.subj.lthr = 1000,
    num.stat.output = 2,
    flag_initiate = FALSE
  )
  expect_equal(result_run[1], 0)
  expect_true(all(is.nan(result_run[-1])))
})

test_that("analyseOneElement.wrap with ctx handles on_error = 'skip'", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src),
                             scalars = list(FD = matrix(1:8, nrow = 2)),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))
  ctx <- .build_wrap_context(modelarray, phen, scalar = "FD")

  expect_warning(
    result <- analyseOneElement.wrap(
      i_element = 1,
      user_fun = function(data, ...) stop("boom"),
      ctx = ctx,
      num.subj.lthr = 0,
      num.stat.output = 3,
      flag_initiate = FALSE,
      on_error = "skip"
    )
  )
  expect_equal(result[1], 0)
  expect_true(all(is.nan(result[-1])))
})

test_that("analyseOneElement.wrap with ctx handles various return types", {
  src <- c("s1", "s2", "s3", "s4")
  modelarray <- methods::new("ModelArray",
                             sources = list(FD = src),
                             scalars = list(FD = matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, byrow = TRUE)),
                             results = list(),
                             path = tempfile(fileext = ".h5")
  )
  phen <- data.frame(source_file = src, age = c(20, 30, 40, 50))
  ctx <- .build_wrap_context(modelarray, phen, scalar = "FD")

  # data.frame return
  df_fun <- function(data, ...) data.frame(m = mean(data$FD))
  result_df <- analyseOneElement.wrap(
    i_element = 1, user_fun = df_fun, ctx = ctx,
    num.subj.lthr = 0, num.stat.output = 2, flag_initiate = FALSE
  )
  expect_equal(length(result_df), 2)

  # named list return
  list_fun <- function(data, ...) list(m = mean(data$FD))
  result_list <- analyseOneElement.wrap(
    i_element = 1, user_fun = list_fun, ctx = ctx,
    num.subj.lthr = 0, num.stat.output = 2, flag_initiate = FALSE
  )
  expect_equal(length(result_list), 2)

  # named atomic return
  atomic_fun <- function(data, ...) c(m = mean(data$FD))
  result_atomic <- analyseOneElement.wrap(
    i_element = 1, user_fun = atomic_fun, ctx = ctx,
    num.subj.lthr = 0, num.stat.output = 2, flag_initiate = FALSE
  )
  expect_equal(length(result_atomic), 2)

  # All three should produce the same values
  expect_equal(result_df[2], result_list[2])
  expect_equal(result_df[2], result_atomic[2])
})


# ==================================================================
# SECTION 4: Initiator search
# ==================================================================

# ---- .find_initiator_element ----

test_that(".find_initiator_element succeeds on middle element", {
  fix <- make_test_modelarray(n_elements = 5, n_subjects = 8, scalars = c("FD"))
  src <- sources(fix$ma)[["FD"]]
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 5, length.out = length(src))
  )
  ctx <- .build_lm_context(FD ~ age, fix$ma, phen, scalar = "FD")

  init_args <- list(
    ctx = ctx,
    var.terms = c("estimate"),
    var.model = c("adj.r.squared"),
    num.subj.lthr = 0,
    on_error = "stop"
  )

  result <- .find_initiator_element(
    analyseOneElement.lm, 5, init_args, verbose = FALSE
  )

  expect_false(is.nan(result$outputs$column_names[1]))
  expect_true(result$i_element >= 1 && result$i_element <= 5)
})

test_that(".find_initiator_element fails when all elements have insufficient subjects", {
  fix <- make_test_modelarray(n_elements = 3, n_subjects = 4, scalars = c("FD"))
  src <- sources(fix$ma)[["FD"]]
  phen <- data.frame(
    source_file = src,
    age = seq(20, by = 5, length.out = length(src))
  )
  ctx <- .build_lm_context(FD ~ age, fix$ma, phen, scalar = "FD")

  init_args <- list(
    ctx = ctx,
    var.terms = c("estimate"),
    var.model = c("adj.r.squared"),
    num.subj.lthr = 1000,  # impossibly high threshold
    on_error = "stop"
  )

  expect_error(
    .find_initiator_element(analyseOneElement.lm, 3, init_args, verbose = FALSE),
    "no elements with sufficient valid"
  )
})


# ==================================================================
# SECTION 5: Top-level integration
#   Verifies ctx flows correctly through ModelArray.gam/wrap
#   using real HDF5-backed data from system.file()
# ==================================================================

# ---- ModelArray.gam integration ----

test_that("ModelArray.gam produces consistent results after context refactoring", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")

  modelarray <- ModelArray(h5_path, scalar_types = c("FD"))
  phenotypes <- read.csv(csv_path)
  element.subset <- 1:5

  result <- ModelArray.gam(
    FD ~ s(age, fx = TRUE) + sex,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = "FD",
    element.subset = element.subset,
    n_cores = 1,
    pbar = FALSE
  )

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(element.subset))
  expect_equal(result$element_id, 0:(length(element.subset) - 1))

  # Verify Intercept naming (not "(Intercept)")
  expect_true(any(grepl("^Intercept\\.", colnames(result))))
  expect_false(any(grepl("\\(Intercept\\)", colnames(result))))

  # Verify smooth term naming uses underscores [10]
  expect_true(any(grepl("^s_age\\.", colnames(result))))

  # No NAs for valid elements
  stat_cols <- setdiff(colnames(result), "element_id")
  expect_true(all(is.finite(as.matrix(result[, stat_cols]))))
})

# ---- ModelArray.wrap integration ----

test_that("ModelArray.wrap produces consistent results after context refactoring", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")

  modelarray <- ModelArray(h5_path, scalar_types = c("FD"))
  phenotypes <- read.csv(csv_path)
  element.subset <- 1:5

  my_fun <- function(data, ...) {
    list(fd_mean = mean(data$FD, na.rm = TRUE))
  }

  result <- ModelArray.wrap(
    FUN = my_fun,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = "FD",
    element.subset = element.subset,
    n_cores = 1,
    pbar = FALSE
  )

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(element.subset))
  expect_equal(result$element_id, 0:(length(element.subset) - 1))
  expect_true("fd_mean" %in% colnames(result))
  expect_true(all(is.finite(result$fd_mean)))
})
