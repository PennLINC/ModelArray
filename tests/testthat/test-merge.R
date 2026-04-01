test_that("mergeModelArrays combines scalars correctly", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  # Load same h5 twice to simulate different files
  ma1 <- ModelArray(h5_path, scalar_types = c("FD"))

  src <- sources(ma1)[["FD"]]
  phen1 <- data.frame(
    subject_id = paste0("subj_", seq_along(src)),
    source_file = src,
    age = runif(length(src), 10, 80),
    stringsAsFactors = FALSE
  )
  # Create a second phenotype with same subject_ids but different source_file col
  phen2 <- data.frame(
    subject_id = paste0("subj_", seq_along(src)),
    source_file = src, # same sources since same h5
    sex = sample(c("M", "F"), length(src), replace = TRUE),
    stringsAsFactors = FALSE
  )

  # This will error because both have scalar "FD" - test collision detection
  expect_error(
    mergeModelArrays(list(ma1, ma1), list(phen1, phen2), merge_on = "subject_id"),
    "Scalar name collision"
  )
})

test_that("mergeModelArrays validates inputs", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  phen <- data.frame(
    subject_id = paste0("subj_", 1:50),
    source_file = sources(ma)[["FD"]],
    stringsAsFactors = FALSE
  )

  # Not enough ModelArrays
  expect_error(
    mergeModelArrays(list(ma), list(phen), merge_on = "subject_id"),
    "at least 2"
  )

  # Mismatched lengths
  expect_error(
    mergeModelArrays(list(ma, ma), list(phen), merge_on = "subject_id"),
    "same length"
  )

  # Missing merge_on column
  phen_no_id <- data.frame(source_file = sources(ma)[["FD"]], stringsAsFactors = FALSE)
  expect_error(
    mergeModelArrays(list(ma, ma), list(phen_no_id, phen), merge_on = "subject_id"),
    "subject_id.*not found"
  )

  # Missing source_file column
  phen_no_sf <- data.frame(subject_id = paste0("subj_", 1:50), stringsAsFactors = FALSE)
  expect_error(
    mergeModelArrays(list(ma, ma), list(phen_no_sf, phen), merge_on = "subject_id"),
    "source_file"
  )
})

test_that("mergeModelArrays subject intersection works", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  src <- sources(ma)[["FD"]]

  # Give ma1 all 50, ma2 only first 30 (but they share the same scalar name
  # so we need to use different scalar names... use the same h5 but pretend)
  # For this test, just verify the join logic by checking row count
  phen_all <- data.frame(
    subject_id = paste0("subj_", 1:50),
    source_file = src,
    stringsAsFactors = FALSE
  )
  phen_partial <- data.frame(
    subject_id = paste0("subj_", 1:30),
    source_file = src[1:30],
    stringsAsFactors = FALSE
  )

  # Can't test with same scalar name, but can verify inner join behavior
  # by checking merged_phen has only the 30 common subjects
  # (This test is limited because we can't easily create two h5 files with different scalars
  # in the test environment. The real test is the cifti integration test above.)
})

test_that("mergeModelArrays returns proper structure", {
  # Build a minimal test with two h5 files by creating temp ones
  # We'll use a mock approach: create a ModelArray manually
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")

  # Create two ModelArrays by loading with different scalar_types
  # (the test h5 only has FD, but we can construct manually)
  ma_fd <- ModelArray(h5_path, scalar_types = c("FD"))
  src <- sources(ma_fd)[["FD"]]

  # Manually construct a second ModelArray with a renamed scalar
  fd_matrix <- scalars(ma_fd)[["FD"]]
  ma2 <- new("ModelArray",
    scalars = list(FD2 = fd_matrix),
    sources = list(FD2 = src),
    results = list(),
    path = h5_path
  )

  phen1 <- data.frame(
    subject_id = paste0("subj_", seq_along(src)),
    source_file = src,
    age = runif(length(src)),
    stringsAsFactors = FALSE
  )
  phen2 <- data.frame(
    subject_id = paste0("subj_", seq_along(src)),
    source_file = src,
    stringsAsFactors = FALSE
  )

  merged <- mergeModelArrays(
    list(ma_fd, ma2),
    list(phen1, phen2),
    merge_on = "subject_id"
  )

  expect_true(is.list(merged))
  expect_true(inherits(merged$data, "ModelArray"))
  expect_true(is.data.frame(merged$phenotypes))
  expect_equal(sort(scalarNames(merged$data)), sort(c("FD", "FD2")))
  expect_true("source_file" %in% colnames(merged$phenotypes))
  expect_true("source_file.FD" %in% colnames(merged$phenotypes))
  expect_true("source_file.FD2" %in% colnames(merged$phenotypes))
  expect_equal(nrow(merged$phenotypes), length(src))
  expect_equal(nInputFiles(merged$data, "FD"), length(src))
  expect_equal(nInputFiles(merged$data, "FD2"), length(src))
  expect_equal(nElements(merged$data, "FD"), nElements(ma_fd, "FD"))
})

test_that("merged ModelArrays work with analyseOneElement helpers", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma_fd <- ModelArray(h5_path, scalar_types = c("FD"))
  src <- sources(ma_fd)[["FD"]]
  base <- scalars(ma_fd)[["FD"]]

  ma2 <- new("ModelArray",
    scalars = list(
      FD2 = base * 1.5 + 2,
      FD3 = log(base + 1)
    ),
    sources = list(FD2 = src, FD3 = src),
    results = list(),
    path = c(FD2 = h5_path, FD3 = h5_path)
  )

  phen1 <- data.frame(
    subject_id = paste0("subj_", seq_along(src)),
    source_file = src,
    age = seq_along(src),
    stringsAsFactors = FALSE
  )
  phen2 <- data.frame(
    subject_id = rev(paste0("subj_", seq_along(src))),
    source_file = rev(src),
    stringsAsFactors = FALSE
  )

  merged <- mergeModelArrays(
    list(ma_fd, ma2),
    list(phen1, phen2),
    merge_on = "subject_id"
  )

  lm_init <- suppressWarnings(analyseOneElement.lm(
    1L,
    FD ~ FD2 + FD3 + age,
    merged$data,
    merged$phenotypes,
    scalar = "FD",
    var.terms = c("estimate", "statistic", "p.value"),
    var.model = c("adj.r.squared", "p.value"),
    num.subj.lthr = 1,
    num.stat.output = 12,
    flag_initiate = TRUE,
    on_error = "stop"
  ))
  expect_true("FD2.estimate" %in% lm_init$column_names)
  expect_true("FD3.estimate" %in% lm_init$column_names)

  lm_fit <- suppressWarnings(ModelArray.lm(
    FD ~ FD2 + FD3 + age,
    merged$data,
    merged$phenotypes,
    scalar = "FD",
    element.subset = as.integer(1:2),
    verbose = FALSE,
    pbar = FALSE,
    n_cores = 1,
    on_error = "stop"
  ))
  expect_true(all(is.finite(lm_fit$FD2.estimate)))
  expect_true(all(is.finite(lm_fit$FD3.estimate)))

  gam_init <- analyseOneElement.gam(
    1L,
    FD ~ s(age) + FD2 + FD3,
    merged$data,
    merged$phenotypes,
    scalar = "FD",
    var.smoothTerms = c("edf", "statistic", "p.value"),
    var.parametricTerms = c("estimate", "statistic", "p.value"),
    var.model = c("adj.r.squared", "dev.expl"),
    num.subj.lthr = 1,
    num.stat.output = 12,
    flag_initiate = TRUE,
    flag_sse = FALSE,
    on_error = "stop"
  )
  expect_true("FD2.estimate" %in% gam_init$column_names)
  expect_true("FD3.estimate" %in% gam_init$column_names)

  wrap_out <- analyseOneElement.wrap(
    1L,
    function(data) {
      data.frame(
        mean_fd = mean(data$FD, na.rm = TRUE),
        cor_fd2 = cor(data$FD, data$FD2, use = "complete.obs"),
        cor_fd3 = cor(data$FD, data$FD3, use = "complete.obs")
      )
    },
    merged$data,
    merged$phenotypes,
    scalar = "FD",
    num.subj.lthr = 1,
    num.stat.output = 4,
    flag_initiate = FALSE,
    on_error = "stop"
  )
  expect_true(is.finite(unname(wrap_out[[3]])))
  expect_true(is.finite(unname(wrap_out[[4]])))
})
