test_that("ModelArray.wrap reproduces ModelArray.lm with a custom FUN", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")

  modelarray <- ModelArray(h5_path, scalar_types = c("FD"), analysis_names = c("my_analysis"))
  phenotypes <- read.csv(csv_path)
  scalar_name <- "FD"
  element.subset <- 1:10

  # Baseline using built-in lm wrapper
  res_lm <- ModelArray.lm(
    FD ~ age + sex,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 1,
    pbar = FALSE
  )

  # User function that mimics ModelArray.lm defaults
  my_lm_fun <- function(data, formula) {
    fit <- stats::lm(formula = formula, data = data)
    td <- broom::tidy(fit)
    if (nrow(td) > 0) {
      td$term[td$term == "(Intercept)"] <- "Intercept"
      td <- td[, intersect(colnames(td), c("term", "estimate", "statistic", "p.value"))]
    }
    gl <- broom::glance(fit)
    gl$term <- "model"
    gl <- gl[, intersect(colnames(gl), c("term", "adj.r.squared", "p.value"))]

    if (nrow(td) > 0) {
      td_w <- tidyr::pivot_wider(
        td,
        names_from = term,
        values_from = c(estimate, statistic, p.value),
        names_glue = "{term}.{.value}"
      )
    } else {
      td_w <- td
    }
    gl_w <- tidyr::pivot_wider(
      gl,
      names_from = term,
      values_from = c(adj.r.squared, p.value),
      names_glue = "{term}.{.value}"
    )
    out <- dplyr::bind_cols(td_w, gl_w)

    # Add FDR p-value columns to match ModelArray.lm defaults
    term_names <- gsub(
      "\\.estimate$|\\.statistic$|\\.p.value$",
      "",
      grep("\\.estimate$|\\.statistic$|\\.p.value$", names(out), value = TRUE)
    )
    term_names <- unique(term_names[term_names != "model"])
    for (nm in term_names) {
      pcol <- paste0(nm, ".p.value")
      if (pcol %in% names(out)) out[[paste0(pcol, ".fdr")]] <- stats::p.adjust(out[[pcol]], method = "fdr")
    }
    if ("model.p.value" %in% names(out)) {
      out[["model.p.value.fdr"]] <- stats::p.adjust(out[["model.p.value"]], method = "fdr")
    }

    tibble::as_tibble(out)
  }

  res_wrap <- ModelArray.wrap(
    FUN = my_lm_fun,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 1,
    pbar = FALSE,
    formula = FD ~ age + sex
  )

  # Same element_id and shape
  expect_equal(res_wrap$element_id, res_lm$element_id)
  expect_true(is.data.frame(res_wrap))

  # Compare on the intersection of columns (names may differ in ordering),
  # excluding FDR columns (which are applied across all elements in ModelArray.lm)
  common_cols <- intersect(colnames(res_lm), colnames(res_wrap))
  common_cols <- setdiff(common_cols, grep("\\.p.value.fdr$", common_cols, value = TRUE))
  # Ensure we have a reasonable overlap (at least the default stats)
  required_patterns <- c(
    "^Intercept\\.estimate$", "^age\\.estimate$", "^sex.*\\.estimate$",
    "^Intercept\\.statistic$", "^age\\.statistic$", "^sex.*\\.statistic$",
    "^Intercept\\.p\\.value$", "^age\\.p\\.value$", "^sex.*\\.p\\.value$",
    "^model\\.adj\\.r\\.squared$", "^model\\.p\\.value$"
  )
  expect_true(all(sapply(required_patterns, function(p) any(grepl(p, common_cols)))))

  expect_equal(res_wrap[, common_cols], res_lm[, common_cols], tolerance = 1e-8)

  rhdf5::h5closeAll()
})

test_that("ModelArray.wrap exposes all scalars and uses intersection threshold", {
  set.seed(7)
  num_elements <- 2L
  num_subj <- 6L
  subj <- paste0("sub-", seq_len(num_subj))
  phen <- data.frame(source_file = subj, z = rnorm(num_subj))

  A <- matrix(rnorm(num_elements * num_subj), nrow = num_elements)
  B <- matrix(rnorm(num_elements * num_subj), nrow = num_elements)
  # introduce NAs in different subjects for element 2
  A[2, 1] <- NA_real_
  B[2, 2] <- NA_real_

  ma <- methods::new(
    "ModelArray",
    scalars = list(A = A, B = B),
    sources = list(A = subj, B = subj),
    results = list(),
    path = ""
  )

  # FUN that checks presence of both scalars in data
  my_fun <- function(data) {
    stopifnot(all(c("A", "B") %in% names(data)))
    c(meanA = mean(data$A, na.rm = TRUE), meanB = mean(data$B, na.rm = TRUE))
  }

  out <- ModelArray.wrap(
    FUN = my_fun,
    data = ma,
    phenotypes = phen,
    scalar = "A",
    element.subset = as.integer(1:2),
    num.subj.lthr.abs = 6L,
    num.subj.lthr.rel = 0,
    verbose = FALSE, pbar = FALSE
  )
  # element 2 (element_id == 1) should be all NaN due to intersection mask
  row2 <- out[out$element_id == 1, , drop = FALSE]
  expect_true(all(is.nan(as.numeric(row2[,-1]))))
})
test_that("ModelArray.wrap errors on collisions and source mismatches for multi-scalars", {
  set.seed(11)
  num_elements <- 2L
  num_subj <- 5L
  subj <- paste0("sub-", seq_len(num_subj))
  phen <- data.frame(source_file = subj, z = rnorm(num_subj))

  A <- matrix(rnorm(num_elements * num_subj), nrow = num_elements)
  B <- matrix(rnorm(num_elements * num_subj), nrow = num_elements)

  ma <- methods::new(
    "ModelArray",
    scalars = list(A = A, B = B),
    sources = list(A = subj, B = subj),
    results = list(),
    path = ""
  )

  # Collision: add a scalar-name column to phenotypes
  phen_bad <- phen
  phen_bad$A <- 1
  expect_error(
    ModelArray.wrap(
      FUN = function(data) tibble::tibble(x = 1),
      data = ma,
      phenotypes = phen_bad,
      scalar = "A",
      element.subset = as.integer(1),
      verbose = FALSE, pbar = FALSE
    ),
    regexp = "Column name collision between phenotypes and scalar names: .* before wrapping\\."
  )

  # Source mismatch for a non-response scalar
  subj_bad <- subj
  subj_bad[1] <- paste0(subj_bad[1], "_bad")
  ma_bad <- methods::new(
    "ModelArray",
    scalars = list(A = A, B = B),
    sources = list(A = subj, B = subj_bad),
    results = list(),
    path = ""
  )
  expect_error(
    ModelArray.wrap(
      FUN = function(data) tibble::tibble(x = 1),
      data = ma_bad,
      phenotypes = phen,
      scalar = "A",
      element.subset = as.integer(1),
      verbose = FALSE, pbar = FALSE
    ),
    regexp = "The source files for scalar 'B' do not match phenotypes\\$source_file\\."
  )
})
test_that("ModelArray.wrap runs simple group means function and writes outputs", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")

  modelarray <- ModelArray(h5_path, scalar_types = c("FD"), analysis_names = c("my_analysis"))
  phenotypes <- read.csv(csv_path)
  scalar_name <- "FD"
  element.subset <- 1:10

  my_group_means_fun <- function(data, group_var, response_name) {
    if (!group_var %in% names(data)) stop("group_var not found in data")
    if (!response_name %in% names(data)) stop("response_name not found in data")
    df <- data
    vals <- df[[response_name]]
    keep <- is.finite(vals)
    df <- df[keep, , drop = FALSE]
    vals <- df[[response_name]]
    grp <- df[[group_var]]
    levs <- sort(unique(grp))
    out <- list()
    for (lv in levs) {
      out[[paste0("mean_", response_name, "_", group_var, "_", lv)]] <- mean(vals[grp == lv], na.rm = TRUE)
    }
    if (length(levs) >= 2) {
      out[[paste0("diff_mean_", response_name, "_", group_var, "_", levs[2], "_minus_", levs[1])]] <-
        out[[paste0("mean_", response_name, "_", group_var, "_", levs[2])]] -
        out[[paste0("mean_", response_name, "_", group_var, "_", levs[1])]]
    }
    tibble::as_tibble(out)
  }

  res_group <- ModelArray.wrap(
    FUN = my_group_means_fun,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 1,
    pbar = FALSE,
    group_var = "sex",
    response_name = "FD"
  )

  expect_true(is.data.frame(res_group))
  expect_equal(res_group$element_id, 0:(length(element.subset) - 1))

  # Expected columns exist (based on available levels in phenotypes$sex)
  levs <- sort(unique(phenotypes$sex))
  expected_cols <- paste0("mean_", scalar_name, "_sex_", levs)
  expect_true(all(expected_cols %in% colnames(res_group)))
  if (length(levs) >= 2) {
    expect_true(paste0("diff_mean_", scalar_name, "_sex_", levs[2], "_minus_", levs[1]) %in% colnames(res_group))
  }

  # Values should be finite (no NAs introduced by wrapper mechanics)
  expect_true(all(is.finite(as.matrix(res_group[, setdiff(colnames(res_group), "element_id")]))))

  rhdf5::h5closeAll()
})
