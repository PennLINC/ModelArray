test_that(".validate_modelarray_input rejects non-ModelArray", {
  expect_error(.validate_modelarray_input("not a modelarray"), "not ModelArray")
  expect_error(.validate_modelarray_input(data.frame()), "not ModelArray")
})

test_that(".validate_element_subset defaults and validates", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  # NULL defaults to all elements
  result <- .validate_element_subset(NULL, ma, "FD")
  expect_equal(length(result), numElementsTotal(ma, "FD"))

  # min < 1
  expect_error(.validate_element_subset(as.integer(0:5), ma, "FD"), "Minimal value")
  # max too large
  expect_error(
    .validate_element_subset(as.integer(c(1, 999999999)), ma, "FD"),
    "Maximal value"
  )
  # non-integer
  expect_error(.validate_element_subset(c(1.0, 2.0), ma, "FD"), "integers")

  # valid input
  result <- .validate_element_subset(as.integer(1:10), ma, "FD")
  expect_equal(result, as.integer(1:10))
})

test_that(".align_phenotypes validates and reorders", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  ma <- ModelArray(h5_path, scalar_types = c("FD"))

  src <- sources(ma)[["FD"]]
  phen <- data.frame(source_file = src, age = seq_along(src))

  # Identical order passes
  result <- .align_phenotypes(ma, phen, "FD")
  expect_identical(result$source_file, src)

  # Reversed order gets reordered
  phen_rev <- phen[rev(seq_len(nrow(phen))), ]
  result <- .align_phenotypes(ma, phen_rev, "FD")
  expect_identical(result$source_file, src)

  # Missing source_file column
  phen_no_sf <- data.frame(age = 1:50)
  expect_error(.align_phenotypes(ma, phen_no_sf, "FD"), "source_file")

  # Length mismatch
  phen_short <- phen[1:10, ]
  expect_error(.align_phenotypes(ma, phen_short, "FD"), "not the same")

  # Non-unique sources
  phen_dup <- phen
  phen_dup$source_file[2] <- phen_dup$source_file[1]
  expect_error(.align_phenotypes(ma, phen_dup, "FD"), "not unique")
})

test_that(".compute_subject_threshold computes max(rel, abs)", {
  phen <- data.frame(x = 1:100)
  # 100 * 0.2 = 20 > 10

  expect_equal(.compute_subject_threshold(phen, 10, 0.2), 20)
  # 100 * 0.05 = 5 < 10
  expect_equal(.compute_subject_threshold(phen, 10, 0.05), 10)
})

test_that(".correct_pvalues inserts corrected columns", {
  df <- data.frame(
    element_id = 0:4,
    Age.p.value = c(0.01, 0.05, 0.1, 0.5, 0.9),
    model.p.value = c(0.001, 0.01, 0.05, 0.1, 0.5)
  )

  # Term-level correction
  result <- .correct_pvalues(df, "Age", "fdr", c("p.value"))
  expect_true("Age.p.value.fdr" %in% colnames(result))
  expect_equal(result$Age.p.value.fdr, p.adjust(df$Age.p.value, "fdr"))

  # Model-level correction
  result2 <- .correct_pvalues(df, "model", "fdr", c("p.value"))
  expect_true("model.p.value.fdr" %in% colnames(result2))

  # "none" returns unchanged
  result3 <- .correct_pvalues(df, "Age", "none", c("p.value"))
  expect_identical(result3, df)

  # No "p.value" in var_list returns unchanged
  result4 <- .correct_pvalues(df, "Age", "fdr", c("estimate"))
  expect_identical(result4, df)
})
