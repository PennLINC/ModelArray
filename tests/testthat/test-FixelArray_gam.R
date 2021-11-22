test_that("test that FixelArray.gam() works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "FixelArray")   # TODO: ask Tinashe
 
  fa <- FixelArray(h5_path,
                    scalar_types = c("FD"),
                    analysis_names = c("my_analysis"))
  
  # h5_path <- paste0(system.file(package = "FixelArray"),
  #                   "inst/extdata/","n50_fixels.h5")
  # fa <- FixelArray(h5_path,
  # scalar_types = c("FD"))
  
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "FixelArray")   # TODO: ask Tinashe
  # csv_path <- paste0(system.file(package = "FixelArray"),
  #                    "inst/extdata/","n50_cohort.csv")
  
  phenotypes <- read.csv(csv_path)
  scalar_name <- "FD"
  var.smoothTerms = c("statistic","p.value")
  var.parametricTerms = c("estimate", "statistic", "p.value")
  var.model = c("dev.expl")
  fixel.subset = 1:10

  mygam <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                          var.smoothTerms = var.smoothTerms,
                          var.parametricTerms = var.parametricTerms,
                          var.model = var.model,
                          n_cores = 1, pbar = FALSE)

  expect_equal(mygam$fixel_id, 0:(max(fixel.subset)-1))  # check output$fixel_id
  expect_true(is.data.frame(mygam))   # should be data.frame
  expect_equal(as.numeric(dim(mygam)), c(length(fixel.subset), 1+1*length(var.smoothTerms) + 2*length(var.parametricTerms) + length(var.model)))  # check the shape

  mygam_default <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                  n_cores = 2, pbar = FALSE)   # default stat outputs
  expect_equal(as.numeric(dim(mygam_default)), c(length(fixel.subset),10))
  expect_equal(mygam, mygam_default)
  
  mygam_fullOutputs <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                      full.outputs = TRUE,    # overwrites the var* arguments below
                                      var.smoothTerms = var.smoothTerms,
                                      var.parametricTerms = var.parametricTerms,
                                      var.model = var.model,
                                      n_cores = 1, pbar = FALSE)
  expect_equal(as.numeric(dim(mygam_fullOutputs)), c(length(fixel.subset),24))
  
  # if there is no explicit parametric term (although there will be term "Intercept") or parametric stat:
  mygam_noExplicitParamTerm <- FixelArray.gam(FD ~ s(age), data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                  n_cores = 2, pbar = FALSE)  
  mygam_noParamStat <- FixelArray.gam(FD ~ s(age)+sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                  var.parametricTerms = c(),
                                  n_cores = 2, pbar = FALSE)  
  # if there is no smooth term or stat:
  mygam_noSmoothTerm <- FixelArray.gam(FD ~ age, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                  n_cores = 2, pbar = FALSE) 
  mygam_noSmoothStat <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                       var.smoothTerms = c(),
                                       n_cores = 2, pbar = FALSE) 
  
  # if there is no model stat:
  mygam_noModelStat <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                      var.model = c(),
                                      n_cores = 2, pbar = FALSE) 
  
  # no any model stat:
  expect_error(FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                      var.smoothTerms = c(), var.parametricTerms = c(), var.model = c(),
                                      n_cores = 2, pbar = FALSE))
})

