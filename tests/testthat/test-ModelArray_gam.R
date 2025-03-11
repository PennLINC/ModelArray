test_that("test that ModelArray.gam() works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")

  scalar_name <- c("FD")
  modelarray <- ModelArray(
    h5_path,
    scalar_types = scalar_name,
    analysis_names = c("my_analysis")
  )
  # h5_path <- paste0(system.file(package = "ModelArray"),
  #                   "inst/extdata/","n50_fixels.h5")
  # scalar_name = c("FD")
  # modelarray <- ModelArray(h5_path, scalar_types = scalar_name, analysis_names = c("my_analysis"))


  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
  # csv_path <- paste0(system.file(package = "ModelArray"),
  #                    "inst/extdata/","n50_cohort.csv")

  phenotypes <- read.csv(csv_path)
  # Create ordered factor with "F" as reference group
  phenotypes$oSex <- ordered(phenotypes$sex, levels = c("F", "M"))
  # phenotypes$sexFactor <- factor(phenotypes$sex, levels = unique(phenotypes$sex))  # factor but not ordered

  # Multiple levels:
  # Create ordered factor with "A" as reference group
  phenotypes$oMultiLevels <- c(
    rep("A", 15),
    rep("B", 15),
    rep("C", 20)
  )
  phenotypes$oMultiLevels <- ordered(phenotypes$oMultiLevels, levels = c("A", "B", "C"))

  var.smoothTerms <- c("statistic", "p.value")
  var.parametricTerms <- c("estimate", "statistic", "p.value")
  var.model <- c("dev.expl", "adj.r.squared")
  # Element subset covers `idx.fixel.gam`
  element.subset <- 11:20

  var.smoothTerms.default <- c("statistic", "p.value")
  var.parametricTerms.default <- c("estimate", "statistic", "p.value")
  var.model.default <- c("dev.expl")

  var.smoothTerms.full <- c("edf", "ref.df", "statistic", "p.value")
  var.parametricTerms.full <- c("estimate", "std.error", "statistic", "p.value")
  # Full list of model statistics to collect
  var.model.full <- c(
    "adj.r.squared", "dev.expl", "sp.criterion", "scale",
    "df", "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs"
  )

  ### generate + load expected results #####
  idx.fixel.gam <- 11
  num.set.seed <- 5
  # Generate the expected results, and get `expected.results`, a list of the expected results
  expected.results <- helper_generate_expect_gam(
    fn.phenotypes = csv_path,
    fn.h5 = h5_path,
    idx.fixel.gam = idx.fixel.gam,
    num.set.seed = num.set.seed
  )

  # Extract the last several characters in a string
  # Ref: https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
  substrRight <- function(x, n) {
    substr(x, nchar(x) - n + 1, nchar(x))
  }

  compare_expected_results <- function(actual, expected, idx.fixel = idx.fixel.gam) {
    ## check if each p.value correction is correct:
    # before removing other rows in "actual":
    col.names <- colnames(actual)
    # Iterate over different correction methods
    for (name.p.adjust in p.adjust.methods) {
      # If any column name's last several characters matches
      # (to avoid name.p.adjust = "BY" vs interation term's name)
      if (name.p.adjust %in% substrRight(col.names, nchar(name.p.adjust))) {
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

    # Also remove any effect size calculations:
    actual <- actual %>% select(-ends_with(c(
      "delta.adj.rsq",
      "partial.rsq"
    )))

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

  ### functions for effect size #####
  # Calculating partial rsq for modelarray
  # This is not comprehensive; does not include additional methods
  # Focus is to check if partial rsq is calculated correctly
  partialRsq_gam_modelarray <- function(i_element, modelarray, phenotypes, scalar, full.formula, reduced.formula) {
    values <- scalars(modelarray)[[scalar]][i_element, ]
    dat <- phenotypes
    dat[[scalar]] <- values

    fullmodel <- mgcv::gam(full.formula, data = dat)
    redmodel <- mgcv::gam(reduced.formula, data = dat)

    # Calculating SSE: used observed y (i.e. excluding observations with NA), and fitted values from model object
    sse.full <- sum((fullmodel$y - fullmodel$fitted.values)^2)
    sse.red <- sum((redmodel$y - redmodel$fitted.values)^2)

    partialRsq <- (sse.red - sse.full) / sse.red
    partialRsq
  }

  wrap_partial_rsq <- function(element.subset, modelarray, phenotypes, scalar, full.formula, reduced.formula) {
    partial.rsq.list <- lapply(
      element.subset, # a list of i_element
      partialRsq_gam_modelarray, # the function
      modelarray, phenotypes, scalar, full.formula, reduced.formula
    )
    partial.rsq.vec <- do.call(rbind, partial.rsq.list)
    partial.rsq.vec <- as.numeric(partial.rsq.vec)
    partial.rsq.vec
  }

  # Calculating delta adjusted rsq for modelarray
  deltaAdjRsq_gam_modelarray <- function(i_element, modelarray, phenotypes, scalar, full.formula, reduced.formula) {
    values <- scalars(modelarray)[[scalar]][i_element, ]
    dat <- phenotypes
    dat[[scalar]] <- values

    fullmodel <- mgcv::gam(full.formula, data = dat)
    redmodel <- mgcv::gam(reduced.formula, data = dat)

    # Calculate adjusted R-squared values for both models
    adjrsq.full <- summary(fullmodel)$r.sq
    adjrsq.red <- summary(redmodel)$r.sq

    # Return the difference in adjusted R-squared values
    delta.adjrsq <- adjrsq.full - adjrsq.red
    delta.adjrsq
  }

  ### basic checks #####
  mygam <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.smoothTerms = var.smoothTerms,
    var.parametricTerms = var.parametricTerms,
    var.model = var.model,
    n_cores = 1, pbar = FALSE
  )

  compare_expected_results(mygam, expected.results[["s-age_sex"]])
  expect_equal(mygam$element_id, (min(element.subset) - 1):(max(element.subset) - 1)) # check output$element_id
  expect_true(is.data.frame(mygam)) # should be data.frame
  expect_equal(as.numeric(dim(mygam)), c(
    length(element.subset),
    1 +
      1 * (length(var.smoothTerms) + 1) +
      2 * (length(var.parametricTerms) + 1) +
      length(var.model)
  )) # check the shape

  mygam_default <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 1, pbar = FALSE
  ) # default stat outputs
  compare_expected_results(mygam_default, expected.results[["s-age_sex"]])
  expect_equal(as.numeric(dim(mygam_default)), c(
    length(element.subset),
    1 +
      1 * (length(var.smoothTerms.default) + 1) +
      2 * (length(var.parametricTerms.default) + 1) +
      length(var.model.default)
  ))

  mygam_fullOutputs <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    full.outputs = TRUE, # overwrites the var* arguments below
    var.smoothTerms = var.smoothTerms,
    var.parametricTerms = var.parametricTerms,
    var.model = var.model,
    n_cores = 1, pbar = FALSE
  )
  compare_expected_results(mygam_fullOutputs, expected.results[["s-age_sex"]])
  expect_equal(as.numeric(dim(mygam_fullOutputs)), c(
    length(element.subset),
    1 +
      1 * (length(var.smoothTerms.full) + 1) +
      2 * (length(var.parametricTerms.full) + 1) +
      length(var.model.full)
  ))

  ### when there is no term or stat output #####
  # if there is no explicit parametric term (although there will be term "Intercept") or parametric stat:
  mygam_noExplicitParamTerm <- ModelArray.gam(FD ~ s(age),
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_noExplicitParamTerm, expected.results[["s-age"]])
  expect_warning(
    mygam_noParamStat <- ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      # there will be warning: p.value was not included in var.parametricTerms,
      # so not to perform its p.value corrections
      var.parametricTerms = c(),
      n_cores = 2,
      pbar = FALSE
    )
  )

  # if there is no smooth term or stat:
  mygam_noSmoothTerm <- ModelArray.gam(
    FD ~ age,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 2,
    pbar = FALSE
  )
  # should without error or warning, as it's just a message
  # expect_output() to test output "there is no smooth term in
  # the requested formula" does not work
  compare_expected_results(mygam_noSmoothTerm, expected.results[["age"]])
  expect_true(c("age.statistic", "age.estimate") %in% colnames(mygam_noSmoothTerm) %>% all())

  expect_warning(
    mygam_noSmoothStat <- ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.smoothTerms = c(),
      n_cores = 2, pbar = FALSE
    )
  )
  compare_expected_results(mygam_noSmoothStat, expected.results[["s-age_sex"]])
  expect_false("age.statistic" %in% colnames(mygam_noSmoothStat))

  # if there is only one var.* not empty:
  # (test with request of reduced model - where there will calls as:
  # c(), c(), c("adj.r.squared"))
  w <- capture_warnings(
    mygam_onlyModelStat <- ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.smoothTerms = c(),
      var.parametricTerms = c(),
      var.model = c("dev.expl"),
      changed.rsq.term.index = c(1),
      n_cores = 2,
      pbar = FALSE
    )
  )
  # warning: p.value was not included in var.smoothTerms OR var.parametricTerms,
  # so not to perform its p.value corrections
  compare_expected_results(mygam_onlyModelStat, expected.results[["s-age_sex"]])
  # all=FALSE: only one element instead of all elements in actual value need to match)
  expect_match(w, "p.value was not included in var.smoothTerms", all = FALSE)
  expect_match(w, "p.value was not included in var.parametricTerms", all = FALSE)
  w <- capture_warnings(
    mygam_onlyParametricTermsStat <- ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      # warning: p.value was not included in var.xxx
      var.smoothTerms = c(),
      var.parametricTerms = c("estimate"),
      var.model = c(),
      changed.rsq.term.index = c(1),
      n_cores = 2,
      pbar = FALSE
    )
  )
  compare_expected_results(mygam_onlyParametricTermsStat, expected.results[["s-age_sex"]])
  # all=FALSE: only one element instead of all elements in actual value need to match)
  expect_match(w, "p.value was not included in var.smoothTerms", all = FALSE)
  expect_match(w, "p.value was not included in var.parametricTerms", all = FALSE)
  # warning: p.value was not included in var.xxx
  w <- capture_warnings(
    mygam_onlySmoothTermsStat <- ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.smoothTerms = c("statistic"),
      var.parametricTerms = c(),
      var.model = c(),
      changed.rsq.term.index = c(1),
      n_cores = 2,
      pbar = FALSE
    )
  )
  compare_expected_results(mygam_onlySmoothTermsStat, expected.results[["s-age_sex"]])
  # all=FALSE: only one element instead of all elements in actual value need to match)
  expect_match(w, "p.value was not included in var.smoothTerms", all = FALSE)
  expect_match(w, "p.value was not included in var.parametricTerms", all = FALSE)
  # there shouldn't be any NA in these outputs: #
  # (otherwise, there is a bug when combining empty tibble in analyseOneElement.gam())
  # if there is any NA in the output data.frame, the result is FALSE
  expect_true((!mygam_onlyParametricTermsStat %>% is.na()) %>% all())
  expect_true((!mygam_onlyModelStat %>% is.na()) %>% all())
  expect_true((!mygam_onlySmoothTermsStat %>% is.na()) %>% all())

  # no any model stat:
  expect_error(ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.smoothTerms = c(), var.parametricTerms = c(), var.model = c(),
    n_cores = 2, pbar = FALSE
  ))


  ## multiple smooth terms:
  # s(age) + s(factorA)   # cannot use s(sex) as sex are characters (M or F), not values!
  mygam_sAge_sFactorA <- ModelArray.gam(FD ~ s(age) + s(factorA),
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_sAge_sFactorA, expected.results[["s-age_s-factorA"]])

  # # s(age + factorA)   # not good example
  # expect_warning(
  #   ModelArray.gam(
  #     FD ~ s(age + factorA),
  #     data = modelarray,
  #     phenotypes = phenotypes,
  #     scalar = scalar_name,
  #     element.subset = element.subset,
  #     n_cores = 2, pbar = FALSE
  #   )
  # )  # note: the column name will be s-age instead of s-age-sex


  ### Test whether the validity of list of var is checked: #####
  expect_error(ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.smoothTerms = c("wrong_name"),
    n_cores = 2, pbar = FALSE
  ))
  # duplicated
  temp <- ModelArray.gam(
    FD ~ s(age) + sex,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.parametricTerms = c(var.parametricTerms, "p.value"),
    n_cores = 2,
    pbar = FALSE
  )
  expect_equal(temp, mygam_default)



  ### different arguments in GAM #####
  ## different settings in formula:
  # s(k=?):
  mygam_k4 <- ModelArray.gam(FD ~ s(age, k = 4) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_k4, expected.results[["s-age-k-4_sex"]])
  expect_false(dplyr::all_equal(
    mygam_default %>% dplyr::select("s_age.statistic"),
    mygam_k4 %>% dplyr::select("s_age.statistic")
  )
  %>% isTRUE()) # should be different when k is different

  # s(fx=TRUE) vs default (fx=FALSE): 1) different stat;
  mygam_fxTRUE <- ModelArray.gam(FD ~ s(age, fx = TRUE) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_fxTRUE, expected.results[["s-age-fx-T_sex"]])
  expect_false(dplyr::all_equal(
    mygam_default %>% dplyr::select("s_age.statistic"),
    mygam_fxTRUE %>% dplyr::select("s_age.statistic")
  )
  %>% isTRUE())

  # TODO: if force saving edf when fx=FALSE: 2) check if it's saved or not

  # s(bs=?)
  mygam_bsCR <- ModelArray.gam(FD ~ s(age, bs = "cr") + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_bsCR, expected.results[["s-age-bs-cr_sex"]])
  expect_false(dplyr::all_equal(
    mygam_default %>% dplyr::select("s_age.statistic"),
    mygam_bsCR %>% dplyr::select("s_age.statistic")
  )
  %>% isTRUE())

  ## different settings in mgcv::gam()'s additional arguments:
  # test if the arguments have been passed into analyseOneElement.gam()
  # handling inputs with NA:
  phenotypes_wNA <- phenotypes
  phenotypes_wNA$age[1] <- NA

  # using default methods in `gam()` to handle NA:
  mygam_agewNA <- ModelArray.gam(
    FD ~ s(age) + sex,
    data = modelarray,
    phenotypes = phenotypes_wNA,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mygam_agewNA, expected.results[["s-age-wNA"]])
  expect_false(dplyr::all_equal(
    mygam_default %>% dplyr::select("s_age.statistic"),
    mygam_agewNA %>% dplyr::select("s_age.statistic")
  )
  %>% isTRUE())
  # using default methods in `gam()` to handle NA + in a formula with interaction terms:
  # WARNING: should not request effect size for `age` without removing NA observations from it!
  # s(x, by=orderedFactor):
  formula <- FD ~ oSex + s(age, k = 4, fx = TRUE) + s(age, by = oSex, fx = TRUE) + factorB # ordered factor
  mygam_agewNA_sby <- ModelArray.gam(
    formula = formula,
    data = modelarray,
    phenotypes = phenotypes_wNA,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mygam_agewNA_sby, expected.results[["oSex_s-age-wNA-k-4-fx-T_s-age-byoSex-fx-T_factorB"]])
  # ti(x) + ti(z) + ti(x,z), where x has NA:
  formula <- FD ~ ti(age, fx = TRUE) + ti(factorB, fx = TRUE) + ti(age, factorB, fx = TRUE) + factorA
  mygam_agewNA_tiInteract <- ModelArray.gam(
    formula = formula,
    data = modelarray,
    phenotypes = phenotypes_wNA,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  )

  compare_expected_results(
    mygam_agewNA_tiInteract,
    expected.results[["ti-age-wNA-fx-T_ti-factorB-fx-T_ti-age-factorB-fx-T_factorA"]]
  )
  # s(x) + s(z) + ti(x,z), where x has NA:
  formula <- FD ~ s(age, fx = TRUE) + s(factorB, fx = TRUE) + ti(age, factorB, fx = TRUE) + factorA
  mygam_agewNA_s_tiInteract <- ModelArray.gam(
    formula = formula,
    data = modelarray,
    phenotypes = phenotypes_wNA,
    scalar = scalar_name,
    element.subset = element.subset,
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(
    mygam_agewNA_s_tiInteract,
    expected.results[["s-age-wNA-fx-T_s-factorB-fx-T_ti-age-factorB-fx-T_factorA"]]
  )

  # using "na.omit" to handle NA:
  mygam_agewNA_na.omit <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = element.subset,
    na.action = "na.omit",
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_agewNA_na.omit, expected.results[["s-age-wNA_na.action-na.omit"]])
  # default option to handle NA should be the same as using `na.omit`:
  expect_true(dplyr::all_equal(
    mygam_agewNA %>% dplyr::select("s_age.statistic"),
    mygam_agewNA_na.omit %>% dplyr::select("s_age.statistic")
  )
  %>% isTRUE())
  # using "na.fail" to handle NA:
  expect_error(ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = element.subset,
    na.action = "na.fail",
    n_cores = 2, pbar = FALSE
  ))

  # method:
  mygam_methodREML <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    method = "REML", # default = "GCV.Cp"
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_methodREML, expected.results[["s-age_sex_method-REML"]])
  expect_false(dplyr::all_equal(
    mygam_default %>% dplyr::select("s_age.statistic"),
    mygam_methodREML %>% dplyr::select("s_age.statistic")
  )
  %>% isTRUE())
  ### Test n_cores, pbar work: ######
  # n_cores = 2:
  mygam_pbarFalse_ncores2 <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  ) # default stat outputs
  expect_equal(mygam_default, mygam_pbarFalse_ncores2)

  # pbar:
  mygam_pbarTrue_ncores2 <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = TRUE
  ) # default stat outputs
  expect_equal(mygam_default, mygam_pbarTrue_ncores2)

  mygam_pbarTrue_ncores1 <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    n_cores = 1, pbar = TRUE
  ) # default stat outputs
  expect_equal(mygam_default, mygam_pbarTrue_ncores1)


  ### Test: p.value correction: #####
  mygam_parametric_pCorrect <- ModelArray.gam(
    FD ~ s(age) + sex,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    correct.p.value.parametricTerms = c("fdr", "bonferroni"),
    n_cores = 2,
    pbar = FALSE
  )
  # this includes check *.p.value.fdr and *.p.value.bonferroni values are correct
  compare_expected_results(
    mygam_parametric_pCorrect,
    expected.results[["s-age_sex"]]
  )
  # check if requested p value corrections exist:
  expect_true(
    c("sexM.p.value.bonferroni", "Intercept.p.value.bonferroni")
    %in% colnames(mygam_parametric_pCorrect) %>% all()
  )

  mygam_smooth_pCorrect <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    correct.p.value.smoothTerms = c("fdr", "bonferroni"),
    n_cores = 2, pbar = FALSE
  )
  # this includes check *.p.value.fdr and *.p.value.bonferroni values are correct
  compare_expected_results(
    mygam_smooth_pCorrect,
    expected.results[["s-age_sex"]]
  )

  # check if requested p value corrections exist:
  expect_true(
    c("s_age.p.value.bonferroni")
    %in% colnames(mygam_smooth_pCorrect) %>% all()
  )

  expect_error(
    ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      correct.p.value.parametricTerms = c("wrong_correction"),
      n_cores = 2,
      pbar = FALSE
    )
  )
  expect_warning(
    temp <- ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.smoothTerms = c("statistic"), # no p.value
      correct.p.value.smoothTerms = c("fdr"),
      n_cores = 2,
      pbar = FALSE
    )
  )
  compare_expected_results(temp, expected.results[["s-age_sex"]])

  ### Test: changed R-squared (delta.adj.rsq, partial.rsq) #####
  #### one term of interest: reduced model will be FD ~ 1 #####
  # also, to test whether the delta.adj.rsq and partial.rsq are calculated correctly
  # also, to test whether without "," in s(), the column name could be correctly "s_age.delta.adj.rsq"
  mygam_changedrsq_oneSmoothTerm <- ModelArray.gam(FD ~ s(age),
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.model = c("dev.expl", "adj.r.squared"),
    changed.rsq.term.index = c(1),
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_changedrsq_oneSmoothTerm, expected.results[["s-age"]])

  mygam_intercept <- ModelArray.gam(FD ~ 1,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.model = c("adj.r.squared"),
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_intercept, expected.results[["intercept"]])
  # delta.adj.rsq:
  expect_equal(
    mygam_changedrsq_oneSmoothTerm$s_age.delta.adj.rsq,
    mygam_changedrsq_oneSmoothTerm$model.adj.r.squared - mygam_intercept$model.adj.r.squared
  )
  # partial.rsq:
  expected <- wrap_partial_rsq(
    element.subset, modelarray, phenotypes, scalar_name,
    FD ~ s(age), FD ~ 1
  )
  expect_equal(mygam_changedrsq_oneSmoothTerm$s_age.partial.rsq, expected)
  # expect that delta.adj.rsq and partial.rsq are highly correlated
  expect_gt(
    cor(mygam_changedrsq_oneSmoothTerm$s_age.partial.rsq, mygam_changedrsq_oneSmoothTerm$s_age.delta.adj.rsq),
    0.95
  )

  ##### the statistics should be consistent, when with or without requesting changed.rsq: #####
  mygam_rsq_smooth_sex <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    full.outputs = TRUE, correct.p.value.smoothTerms = "fdr", correct.p.value.parametricTerms = "fdr",
    changed.rsq.term.index = c(1),
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_rsq_smooth_sex, expected.results[["s-age_sex"]])
  # compared to statistics in full outputs:
  colnames_intersect <- intersect(
    mygam_rsq_smooth_sex %>% colnames(),
    mygam_fullOutputs %>% colnames()
  )
  expect_equal(
    mygam_rsq_smooth_sex %>% select(colnames_intersect),
    mygam_fullOutputs %>% select(colnames_intersect)
  )
  # compared to fdr:
  expect_equal(
    mygam_rsq_smooth_sex$sex.p.value.fdr,
    mygam_parametric_pCorrect$sex.p.value.fdr
  )
  expect_equal(
    mygam_rsq_smooth_sex$s_age.p.value.fdr,
    mygam_smooth_pCorrect$s_age.p.value.fdr
  )


  #### more than one term of interest; also, parametric term or smooth term: #####
  mygam_changedRsq_twoSmoothTerm <- ModelArray.gam(FD ~ factorB + s(age) + s(factorA),
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    # var.model = c("dev.expl", "adj.r.squared"),
    full.outputs = TRUE,
    changed.rsq.term.index = c(1, 2, 3),
    correct.p.value.smoothTerms = "fdr", correct.p.value.parametricTerms = "fdr",
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_changedRsq_twoSmoothTerm, expected.results[["factorB_s-age_s-factorA"]])
  mygam_rsq_red1 <- ModelArray.gam(
    FD ~ factorB + s(factorA),
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.model = c("adj.r.squared"),
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(
    mygam_rsq_red1,
    expected.results[["factorB_s-factorA"]]
  )
  mygam_rsq_red2 <- ModelArray.gam(
    FD ~ factorB + s(age),
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.model = c("adj.r.squared"),
    n_cores = 2,
    pbar = FALSE
  )
  compare_expected_results(mygam_rsq_red2, expected.results[["factorB_s-age"]])
  mygam_chRsq_2SmTrm_red3 <- ModelArray.gam(FD ~ s(age) + s(factorA),
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    var.model = c("adj.r.squared"),
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_chRsq_2SmTrm_red3, expected.results[["s-age_s-factorA"]])
  # delta.adj.rsq:
  expect_equal(
    mygam_changedRsq_twoSmoothTerm$s_age.delta.adj.rsq,
    mygam_changedRsq_twoSmoothTerm$model.adj.r.squared - mygam_rsq_red1$model.adj.r.squared
  )
  # partial.rsq:
  # test red1:
  expected <- wrap_partial_rsq(
    element.subset, modelarray, phenotypes, scalar_name,
    FD ~ factorB + s(age) + s(factorA), FD ~ factorB + s(factorA)
  )
  expect_equal(mygam_changedRsq_twoSmoothTerm$s_age.partial.rsq, expected)
  # test red2:
  expected <- wrap_partial_rsq(
    element.subset, modelarray, phenotypes, scalar_name,
    FD ~ factorB + s(age) + s(factorA), FD ~ factorB + s(age)
  )
  expect_equal(mygam_changedRsq_twoSmoothTerm$s_factorA.partial.rsq, expected)
  # test red3:
  expected <- wrap_partial_rsq(
    element.subset, modelarray, phenotypes, scalar_name,
    FD ~ factorB + s(age) + s(factorA), FD ~ s(age) + s(factorA)
  )
  expect_equal(mygam_changedRsq_twoSmoothTerm$factorB.partial.rsq, expected)


  ##### the statistics should be consistent, when with or without requesting changed.rsq: #####
  mygam_smooth_base <- ModelArray.gam(FD ~ factorB + s(age) + s(factorA),
    data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
    # var.model = c("dev.expl", "adj.r.squared"),
    full.outputs = TRUE,
    correct.p.value.smoothTerms = "fdr", correct.p.value.parametricTerms = "fdr",
    n_cores = 2, pbar = FALSE
  )
  compare_expected_results(mygam_smooth_base, expected.results[["factorB_s-age_s-factorA"]])
  # compared to statistics in full outputs:
  colnames_intersect <- intersect(
    mygam_changedRsq_twoSmoothTerm %>% colnames(),
    mygam_smooth_base %>% colnames()
  )
  expect_equal(
    mygam_changedRsq_twoSmoothTerm %>% select(colnames_intersect),
    mygam_smooth_base %>% select(colnames_intersect)
  )
  # compared to fdr:
  expect_equal(
    mygam_changedRsq_twoSmoothTerm$sex.p.value.fdr,
    mygam_smooth_base$sex.p.value.fdr
  )
  expect_equal(
    mygam_changedRsq_twoSmoothTerm$s_age.p.value.fdr,
    mygam_smooth_base$s_age.p.value.fdr
  )

  #### test out the functions for generating gam functions: #####
  ## Formula #1:
  myFormula_1 <- generator_gamFormula_factorXsmooth(
    response.var = "FD", factor.var = "oSex", smooth.var = "age",
    phenotypes = phenotypes
  )
  myFormula_1$formula # requires visually check

  myFormula_2 <- generator_gamFormula_factorXsmooth(
    response.var = "FD", factor.var = "sex", smooth.var = "age",
    phenotypes = phenotypes, reference.group = "F"
  )
  myFormula_2$formula # requires visually check
  expect_true("osex" %in% colnames(myFormula_2$phenotypes))
  phenotypes_updated <- myFormula_2$phenotypes
  osex.class <- class(phenotypes_updated[["osex"]])
  expect_true((length(osex.class) == 2) & (osex.class[1] == "ordered") & (osex.class[2] == "factor"))

  expect_error(generator_gamFormula_factorXsmooth(
    response.var = "FD", factor.var = "sex", smooth.var = "age",
    phenotypes = phenotypes
  )) # not ordered factor, and did not provide reference.group

  # change k and fx:
  myFormula_3 <- generator_gamFormula_factorXsmooth(
    response.var = "FD", factor.var = "oSex", smooth.var = "age",
    phenotypes = phenotypes, fx = FALSE, k = 4
  )
  myFormula_3$formula # requires visually check


  ## Formula #2:
  myFormula_4 <- generator_gamFormula_continuousInteraction(
    response.var = "FD",
    cont1.var = "age",
    cont2.var = "factorA"
  )
  myFormula_4 # requires visually check

  # change k and fx:
  myFormula_5 <- generator_gamFormula_continuousInteraction(
    response.var = "FD", cont1.var = "age", cont2.var = "factorA",
    fx = FALSE, k = 3
  )
  myFormula_5 # requires visually check


  ### source file list sanity check #####
  # not the same length:
  phenotypes_wrong1 <- phenotypes[-c(1), ]
  expect_error(ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes_wrong1, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  ))

  # not the same order --> will be handled!
  phenotypes_swap <- phenotypes
  temp <- phenotypes_swap[2, ] # swap row1 and row2
  phenotypes_swap[2, ] <- phenotypes_swap[1, ]
  phenotypes_swap[1, ] <- temp
  mygam_testSwap <- ModelArray.gam(FD ~ s(age) + sex,
    data = modelarray, phenotypes = phenotypes_swap, scalar = scalar_name, element.subset = element.subset,
    n_cores = 2, pbar = FALSE
  )
  expect_equal(mygam_testSwap, mygam_default)

  # change one row --> cannot be matched --> error:
  phenotypes_wrong2 <- phenotypes
  phenotypes_wrong2[1, "source_file"] <- "wrong_file"
  expect_error(
    ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes_wrong2,
      scalar = scalar_name,
      element.subset = element.subset,
      n_cores = 2,
      pbar = FALSE
    )
  )




  ### debugging:
  #  Error in term[i] <- attr(terms(reformulate(term[i])), "term.labels") :
  #   replacement has length zero
  # may because of invalid arguments in smooth term (e.g. d in s(age, d = 1))
})

test_that("ModelArray.gam works with effect size calculations", {
  # Test partial R-squared calculation
  temp <- ModelArray.gam(
    FD ~ s(age) + sex,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = c("estimate"),
    var.model = var.model,
    effect.size.terms = c("partial.rsq", "delta.adj.rsq"),
    n_cores = 2,
    pbar = FALSE
  )

  # Check if effect size columns exist
  expect_true("sex.partial.rsq" %in% colnames(temp))
  expect_true("sex.delta.adj.rsq" %in% colnames(temp))

  # Calculate expected values for comparison
  expected.partial.rsq <- partialRsq_gam_modelarray(
    idx.fixel.gam,
    modelarray,
    phenotypes,
    scalar_name,
    FD ~ s(age) + sex,
    FD ~ s(age)
  )

  expected.delta.adjrsq <- deltaAdjRsq_gam_modelarray(
    idx.fixel.gam,
    modelarray,
    phenotypes,
    scalar_name,
    FD ~ s(age) + sex,
    FD ~ s(age)
  )

  # Compare calculated values with expected values
  expect_equal(
    temp$sex.partial.rsq[temp$element_id == idx.fixel.gam - 1],
    expected.partial.rsq
  )

  expect_equal(
    temp$sex.delta.adj.rsq[temp$element_id == idx.fixel.gam - 1],
    expected.delta.adjrsq
  )
})

test_that("ModelArray.gam works with p-value corrections", {
  # Test multiple p-value correction methods
  temp <- ModelArray.gam(
    FD ~ s(age) + sex,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = c("p.value"),
    var.model = var.model,
    correct.p.value.terms = c("fdr", "bonferroni"),
    n_cores = 2,
    pbar = FALSE
  )

  # Check if corrected p-value columns exist
  expect_true("sex.p.value.fdr" %in% colnames(temp))
  expect_true("sex.p.value.bonferroni" %in% colnames(temp))

  # Verify that p-value corrections are calculated correctly
  compare_expected_results(temp, expected.gam, idx.fixel.gam)
})

test_that("ModelArray.gam handles errors appropriately", {
  # Test error for incorrect p-value correction method
  expect_error(
    ModelArray.gam(
      FD ~ s(age) + sex,
      data = modelarray,
      phenotypes = phenotypes,
      scalar = scalar_name,
      element.subset = element.subset,
      var.terms = c("p.value"),
      var.model = var.model,
      correct.p.value.terms = c("invalid_method"),
      n_cores = 2,
      pbar = FALSE
    ),
    "Invalid p-value correction method"
  )

  # Test error for missing required arguments
  expect_error(
    ModelArray.gam(
      FD ~ s(age) + sex,
      phenotypes = phenotypes,
      scalar = scalar_name
    ),
    "argument .* is missing"
  )
})

test_that("ModelArray.gam handles NA values correctly", {
  # Create test data with NA values
  phenotypes_with_na <- phenotypes
  phenotypes_with_na$age[1] <- NA

  # Test handling of NA values in covariates
  temp <- ModelArray.gam(
    FD ~ s(age) + sex,
    data = modelarray,
    phenotypes = phenotypes_with_na,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = c("estimate", "p.value"),
    var.model = var.model,
    n_cores = 2,
    pbar = FALSE
  )

  # Verify that results exclude NA observations
  expect_equal(nrow(temp), length(element.subset))
  expect_true(all(!is.na(temp$sex.estimate)))
  expect_true(all(!is.na(temp$sex.p.value)))
})

test_that("ModelArray.gam works with additional covariates", {
  # Test model with multiple covariates
  temp <- ModelArray.gam(
    FD ~ s(age) + sex + factorA + factorB,
    data = modelarray,
    phenotypes = phenotypes,
    scalar = scalar_name,
    element.subset = element.subset,
    var.terms = c("estimate", "p.value"),
    var.model = var.model,
    n_cores = 2,
    pbar = FALSE
  )

  # Check if all covariate columns exist
  expect_true(all(c(
    "sex.estimate", "sex.p.value",
    "factorA.estimate", "factorA.p.value",
    "factorB.estimate", "factorB.p.value"
  ) %in% colnames(temp)))

  # Verify results match expected values
  compare_expected_results(temp, expected.gam.with.factors, idx.fixel.gam)
})
