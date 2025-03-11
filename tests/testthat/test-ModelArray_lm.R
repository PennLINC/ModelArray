test_that("ModelArray.lm() works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")

  modelarray <- ModelArray(
    h5_path,
    "scalar_types" = c("FD"),
    analysis_names = c("my_analysis")
  )

  # h5_path <- paste0(system.file(package = "ModelArray"),
  #                   "inst/extdata/","n50_fixels.h5")
  # modelarray <- ModelArray(h5_path,
  # scalar_types = c("FD"))

  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
  # csv_path <- paste0(system.file(package = "ModelArray"),
  #                    "inst/extdata/","n50_cohort.csv")

  phenotypes <- read.csv(csv_path)
  scalar_name <- "FD"
  # List of columns to keep
  var.terms <- c("estimate", "p.value")
  var.terms.full <- c("estimate", "p.value", "std.error", "statistic")
  var.model <- c("r.squared", "p.value", "AIC")

  default.var.terms <- c("estimate", "statistic", "p.value")
  default.var.model <- c("adj.r.squared", "p.value")
  var.terms.full <- c("estimate", "std.error", "statistic", "p.value")
  # Full list of model statistics to collect
  var.model.full <- c(
    "r.squared", "adj.r.squared", "sigma", "statistic", "p.value", "df",
    "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs"
  )

  element.subset <- 1:10

  ### generate + load expected results #####
  idx.fixel.lm <- 1
  # Used when generating random numbers in expected results and below
  num.set.seed <- 5
  # Generate the expected results, and get `expected.results`, a list of the expected results
  expected.results <- helper_generate_expect_lm(
    fn.phenotypes = csv_path,
    fn.h5 = h5_path,
    idx.fixel.lm = idx.fixel.lm,
    num.set.seed = num.set.seed
  )

  #' @param idx.fixel starts from 1
  compare_expected_results <- function(actual, expected, idx.fixel = idx.fixel.lm) {
    ## check if each p.value correction is correct:
    # before removing other rows in "actual":
    col.names <- colnames(actual)
    # Iterate over different correction methods
    for (name.p.adjust in p.adjust.methods) {
      # If it contains this name.p.adjust
      if (stringr::str_detect(col.names, name.p.adjust) %>% any()) {
        # Get index of columns containing this name.p.adjust
        list.idx <- stringr::str_which(col.names, name.p.adjust)
        # For each corrected term/model's p.value
        for (idx in list.idx) {
          thecolname <- col.names[idx]
          # Remove the name.p.adjust to get the p.value's name
          thecolname.pvalue <- stringr::str_remove_all(
            thecolname,
            paste0(".", name.p.adjust)
          )
          expect_equal(
            actual[[thecolname.pvalue]] %>% stats::p.adjust(method = name.p.adjust),
            actual[[thecolname]]
          )
        }
      }
    }

    # Now, remove any fdr, etc corrections:
    actual <- actual %>% select(-ends_with(p.adjust.methods))
    col.names <- colnames(actual)

    # Select the corresponding row (which has the expected value):
    actual <- actual[actual$element_id == idx.fixel - 1, ]

    flag.belong <- col.names %in% colnames(expected) %>% all()
    if (flag.belong == FALSE) {
      stop("not all columns in actual data.frame are in expected data.frame")
    }

    ## test if actual = expected values:
    expect_equal(
      actual,
      expected %>% select(col.names)
    )
  }

  ### basic check #####
  mylm <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    n_cores = 1,
    pbar = FALSE
  )

  compare_expected_results(mylm, expected.results$age)
  # Check output$element_id
  expect_equal(
    mylm$element_id,
    0:(length(element.subset) - 1)
  )
  # Should be data.frame
  expect_true(
    is.data.frame(mylm)
  )
  # Check shape
  expect_equal(
    as.numeric(dim(mylm)),
    c(
      length(element.subset),
      1 +
        # 2*: intercept + age; +1 in length: fdr is default
        2 * (length(var.terms) + 1) +
        # 1*: lm model; +1 in length: fdr is default
        1 * (length(var.model) + 1)
    )
  )
  # FDR is saved by default
  expect_true(
    c("age.p.value.fdr", "Intercept.p.value.fdr", "model.p.value.fdr") %in% colnames(mylm) %>% all()
  )

  # Default full.outputs and var.*
  mylm_default <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mylm_default, expected.results$age)
  expect_equal(
    as.numeric(dim(mylm_default)),
    c(
      length(element.subset),
      1 +
        2 * (length(default.var.terms) + 1) +
        1 * (length(default.var.model) + 1)
    )
  )

  # Default: FALSE
  mylm_fullOutputs <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    full.outputs = TRUE,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mylm_fullOutputs, expected.results$age)
  expect_equal(
    as.numeric(dim(mylm_fullOutputs)),
    c(
      length(element.subset),
      1 +
        2 * (length(var.terms.full) + 1) +
        1 * (length(var.model.full) + 1)
    )
  )
  # FDR is saved by default
  expect_true(
    c("age.p.value.fdr", "Intercept.p.value.fdr", "model.p.value.fdr") %in% colnames(mylm) %>% all()
  )

  mylm_age_sex <- ModelArray.lm(
    FD ~ age + sex,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mylm_age_sex, expected.results$age_sex)
  # Check output$element_id
  expect_equal(
    mylm_age_sex$element_id,
    0:(length(element.subset) - 1)
  )
  expect_equal(
    as.numeric(dim(mylm_age_sex)),
    c(
      length(element.subset),
      1 +
        3 * (length(var.terms) + 1) +
        1 * (length(var.model) + 1)
    )
  )

  # Expect not identical between two models
  expect_false(
    all.equal(
      mylm %>% dplyr::select("age.estimate"),
      mylm_age_sex %>% dplyr::select("age.estimate")
    )
    %>% isTRUE()
  )

  ## Test whether the validity of list of var is checked:
  # Test misspelling
  expect_error(
    ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = c("estimator"),
      var.model = var.model,
      n_cores = 2,
      pbar = FALSE
    )
  )
  # Test handling of duplicates
  temp <- (ModelArray.lm(FD ~ age,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.terms = var.terms,
    var.model = c(var.model, "AIC"),
    n_cores = 2, pbar = FALSE
  ))
  expect_equal(mylm, temp)


  ### Test n_cores, pbar work #####
  # Test n_cores=2
  mylm_ncores2 <- ModelArray.lm(FD ~ age,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    n_cores = 2, pbar = FALSE
  )
  expect_equal(mylm, mylm_ncores2)
  # Test pbar=TRUE & n_cores=1
  mylm_pbarTRUE_ncores1 <- ModelArray.lm(FD ~ age,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    n_cores = 1, pbar = TRUE
  )
  expect_equal(mylm, mylm_pbarTRUE_ncores1)
  # Test pbar=TRUE & n_cores=2
  mylm_pbarTRUE_ncores2 <- ModelArray.lm(FD ~ age,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    n_cores = 2, pbar = TRUE
  )
  expect_equal(mylm, mylm_pbarTRUE_ncores2)

  ## Different output statistics #####
  # Expect warning because p.value was not included in var.terms, so not to perform its p.value corrections
  expect_warning(
    mylm_noTermsOutput <- ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = c(),
      var.model = var.model,
      n_cores = 1,
      pbar = FALSE
    )
  )
  compare_expected_results(mylm_noTermsOutput, expected.results$age)
  # Check shape
  expect_equal(
    as.numeric(dim(mylm_noTermsOutput)),
    c(length(element.subset), 1 + length(var.model) + 1)
  )
  # Check for NAs in output (would indicate bug when combining empty tibble in analyseOneElement.lm())
  expect_true((!mylm_noTermsOutput %>% is.na()) %>% all())

  # Expect warning because p.value was not included in var.model, so not to perform its p.value corrections
  expect_warning(
    mylm_noModelOutput <- ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = var.terms,
      var.model = c(),
      n_cores = 1,
      pbar = FALSE
    )
  )
  compare_expected_results(mylm_noModelOutput, expected.results$age)
  # Check shape
  expect_equal(
    as.numeric(dim(mylm_noModelOutput)),
    c(length(element.subset), 1 + 2 * (length(var.terms) + 1))
  )
  expect_true((!mylm_noModelOutput %>% is.na()) %>% all())

  # Error if both var.* are empty
  expect_error(
    mylm_noModelOutput <- ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = c(),
      var.model = c(),
      n_cores = 1,
      pbar = FALSE
    )
  )

  ## Whether to correct p.values: #####
  # Test terms corrections
  mylm_corr_pvalues_1 <- ModelArray.lm(FD ~ age,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.terms = var.terms, var.model = var.model,
    correct.p.value.terms = c("fdr", "bonferroni"),
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mylm_corr_pvalues_1, expected.results$age)
  # Check if requested p value corrections exist
  expect_true(
    c("age.p.value.bonferroni", "Intercept.p.value.bonferroni")
    %in% colnames(mylm_corr_pvalues_1)
    %>% all()
  )

  # Test model corrections
  mylm_corr_pvalues_2 <- ModelArray.lm(FD ~ age,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.terms = var.terms, var.model = var.model,
    correct.p.value.model = c("fdr", "bonferroni"),
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mylm_corr_pvalues_2, expected.results$age)
  # Check if requested p value corrections exist
  expect_true(
    c("model.p.value.bonferroni") %in% colnames(mylm_corr_pvalues_2) %>% all()
  )

  # Test terms + only bonferroni, no fdr
  mylm_corr_pvalues_3 <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    correct.p.value.terms = c("bonferroni"),
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mylm_corr_pvalues_3, expected.results$age)
  # Check correct.p.value.terms: bonferroni is in, but not fdr
  expect_true(
    c(
      "age.p.value.bonferroni",
      "Intercept.p.value.bonferroni",
      "model.p.value.fdr"
    ) %in% colnames(mylm_corr_pvalues_3) %>% all()
  )
  # All false
  expect_true(
    !(c("age.p.value.fdr", "Intercept.p.value.fdr") %in% colnames(mylm_corr_pvalues_3)) %>% all()
  )

  # Test model + only bonferroni, no fdr
  mylm_corr_pvalues_4 <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    correct.p.value.model = c("bonferroni"),
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mylm_corr_pvalues_4, expected.results$age)
  # Check correct.p.value.model: bonferroni is in, but not fdr
  expect_true(
    c(
      "age.p.value.fdr",
      "model.p.value.bonferroni"
    ) %in% colnames(mylm_corr_pvalues_4) %>% all()
  )
  # All false
  expect_true(
    !(c("model.p.value.fdr") %in% colnames(mylm_corr_pvalues_4)) %>% all()
  )

  # Test warning when p.value not provided in var.terms
  expect_warning(
    temp <- ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = c("estimate"),
      var.model = var.model,
      correct.p.value.terms = c("fdr", "bonferroni"),
      n_cores = 2,
      pbar = FALSE
    )
  )
  compare_expected_results(temp, expected.results$age)

  # Test warning when p.value not provided in var.model
  expect_warning(
    temp <- ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = var.terms,
      var.model = c("AIC"),
      correct.p.value.model = c("fdr", "bonferroni"),
      n_cores = 2,
      pbar = FALSE
    )
  )
  compare_expected_results(temp, expected.results$age)

  ## Test other variables as covariates #####
  # factorA is literally correlated with age; factorB is another random variable
  # Factor A is fully correlated with age, expecting testing results are NA
  mylm_age_factorA <- ModelArray.lm(
    FD ~ age + factorA,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms.full,
    var.model = var.model,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mylm_age_factorA, expected.results$age_factorA)
  # Test that $element_id is 0:(length(element.subset)-1)
  expect_equal(
    mylm_age_factorA$element_id,
    c(0:(length(element.subset) - 1))
  )

  # Anything of factorA should be NA
  expect_true(
    all(is.na(mylm_age_factorA$factorA.estimate))
  )
  expect_true(
    all(is.na(mylm_age_factorA$factorA.p.value))
  )
  expect_true(
    all(is.na(mylm_age_factorA$factorA.std.error))
  )
  expect_true(
    all(is.na(mylm_age_factorA$factorA.statistic))
  )

  # Factor B is not correlated with age, so not expecting NA
  mylm_age_factorB <- ModelArray.lm(
    FD ~ age + factorB,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms.full,
    var.model = var.model,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mylm_age_factorB, expected.results$age_factorB)
  # Check factorB results are not NA
  expect_false(
    all(is.na(mylm_age_factorB$factorB.estimate))
  )
  expect_false(
    all(is.na(mylm_age_factorB$factorB.p.value))
  )
  expect_false(
    all(is.na(mylm_age_factorB$factorB.std.error))
  )
  expect_false(
    all(is.na(mylm_age_factorB$factorB.statistic))
  )

  ## Test different optional arguments of lm #####
  # Test handling inputs with NA
  phenotypes_wNA <- phenotypes
  phenotypes_wNA$age[1] <- NA

  # Test default `lm()` method for handling NA
  mylm_phntyp_naActDef <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes_wNA,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mylm_phntyp_naActDef, expected.results[["age_phenotypeswNA"]])

  # Check differences in results after changing one value to NA
  expect_false(
    all.equal(
      mylm %>% dplyr::select(age.estimate),
      mylm_phntyp_naActDef %>% dplyr::select(age.estimate)
    )
    %>% isTRUE()
  )

  # Test different (not-default) `na.action` with inputs with NA
  # Test error if `na.action = "na.fail"`
  expect_error(
    ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes_wNA,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = var.terms,
      var.model = var.model,
      n_cores = 1,
      pbar = FALSE,
      na.action = "na.fail"
    )
  )
  # Test with different n_cores
  expect_error(
    ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes_wNA,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = var.terms,
      var.model = var.model,
      n_cores = 2,
      pbar = FALSE,
      na.action = "na.fail"
    )
  )
  # Test with different pbar settings
  expect_error(
    ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes_wNA,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = var.terms,
      var.model = var.model,
      n_cores = 1,
      pbar = TRUE,
      na.action = "na.fail"
    )
  )
  expect_error(
    ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes_wNA,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = var.terms,
      var.model = var.model,
      n_cores = 2,
      pbar = TRUE,
      na.action = "na.fail"
    )
  )

  # Test na.action = "na.omit"
  mylm_phenotypes_naActionOmit <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes_wNA,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    n_cores = 2,
    pbar = FALSE,
    na.action = "na.omit"
  )
  compare_expected_results(mylm_phenotypes_naActionOmit, expected.results[["age_phenotypeswNA_na.action-na.omit"]])
  # Check differences in results after changing one value to NA
  expect_false(
    all.equal(
      mylm %>% dplyr::select(age.estimate),
      mylm_phenotypes_naActionOmit %>% dplyr::select(age.estimate)
    ) %>% isTRUE()
  )
  # Check default option to handle NA is same as using `na.omit`
  expect_true(
    all.equal(
      mylm_phntyp_naActDef %>% dplyr::select(age.estimate),
      mylm_phenotypes_naActionOmit %>% dplyr::select(age.estimate)
    )
    %>% isTRUE()
  )

  # Test if "weights" have been successfully passed into lm
  mylm_weights1 <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    pbar = FALSE,
    n_cores = 2,
    weights = rep(1, nrow(phenotypes))
  )
  expect_equal(mylm, mylm_weights1)

  # Test random weights
  set.seed(num.set.seed)
  mylm_weightsRnorm <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = var.terms,
    var.model = var.model,
    n_cores = 2,
    pbar = FALSE,
    weights = abs(rnorm(nrow(phenotypes)))
  )
  compare_expected_results(mylm_weightsRnorm, expected.results[["age_weights-random"]])
  expect_false(
    all.equal(mylm, mylm_weightsRnorm) %>% isTRUE()
  )

  ### Test source file list sanity check #####
  # Test not the same length
  phenotypes_wrong1 <- phenotypes[-c(1), ]
  expect_error(
    ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes_wrong1,
      scalar = scalar_name,
      element.subset = element.subset,
      n_cores = 2,
      pbar = FALSE
    )
  )

  # Test not the same order (will be handled)
  phenotypes_swap <- phenotypes
  temp <- phenotypes_swap[2, ]
  phenotypes_swap[2, ] <- phenotypes_swap[1, ]
  phenotypes_swap[1, ] <- temp
  mylm_testSwap <- ModelArray.lm(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes_swap,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 2,
    pbar = FALSE
  )
  expect_equal(mylm_testSwap, mylm_default)

  # Test change one row (cannot be matched - error)
  phenotypes_wrong2 <- phenotypes
  phenotypes_wrong2[1, "source_file"] <- "wrong_file"
  expect_error(
    ModelArray.lm(
      FD ~ age,
      data = modelarray,
      phenotypes = phenotypes_wrong2,
      scalar = scalar_name,
      element.subset = element.subset,
      n_cores = 2,
      pbar = FALSE
    )
  )

  rhdf5::h5closeAll()
})
