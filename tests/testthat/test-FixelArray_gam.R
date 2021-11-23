test_that("test that FixelArray.gam() works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "FixelArray")   # TODO: ask Tinashe
 
  fa <- FixelArray(h5_path,
                    scalar_types = c("FD"),
                    analysis_names = c("my_analysis"))
  
  # h5_path <- paste0(system.file(package = "FixelArray"),
  #                   "inst/extdata/","n50_fixels.h5")
  # scalar_name = c("FD")
  # fa <- FixelArray(h5_path, scalar_types = scalar_name, analysis_names = c("my_analysis"))
  

  csv_path <- system.file("extdata", "n50_cohort.csv", package = "FixelArray")   # TODO: ask Tinashe
  # csv_path <- paste0(system.file(package = "FixelArray"),
  #                    "inst/extdata/","n50_cohort.csv")
  
  phenotypes <- read.csv(csv_path)
  var.smoothTerms = c("statistic","p.value")
  var.parametricTerms = c("estimate", "statistic", "p.value")
  var.model = c("dev.expl", "AIC")
  fixel.subset = 1:10

  ### basic checks #####
  mygam <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                          var.smoothTerms = var.smoothTerms,
                          var.parametricTerms = var.parametricTerms,
                          var.model = var.model,
                          n_cores = 1, pbar = FALSE)

  expect_equal(mygam$fixel_id, 0:(max(fixel.subset)-1))  # check output$fixel_id
  expect_true(is.data.frame(mygam))   # should be data.frame
  expect_equal(as.numeric(dim(mygam)), c(length(fixel.subset), 1+1*length(var.smoothTerms) + 2*length(var.parametricTerms) + length(var.model)))  # check the shape

  mygam_default <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                  n_cores = 1, pbar = FALSE)   # default stat outputs
  expect_equal(as.numeric(dim(mygam_default)), c(length(fixel.subset),10))

  mygam_fullOutputs <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                      full.outputs = TRUE,    # overwrites the var* arguments below
                                      var.smoothTerms = var.smoothTerms,
                                      var.parametricTerms = var.parametricTerms,
                                      var.model = var.model,
                                      n_cores = 1, pbar = FALSE)
  expect_equal(as.numeric(dim(mygam_fullOutputs)), c(length(fixel.subset),24))
  
  ### when there is no term or stat output #####
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
  
  
  ## multiple smooth terms:
  # s(age) + s(factorA)   # cannot use s(sex) as sex are characters (M or F), not values!
  mygam_sAge_sFactorA <- FixelArray.gam(FD ~ s(age) + s(factorA), data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                n_cores = 2, pbar = FALSE)
  
  # s(age + factorA)   # should be rare case
  expect_warning(FixelArray.gam(FD ~ s(age + factorA), data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                 n_cores = 2, pbar = FALSE))  # TODO: fix this problem: the column name will be s-age instead of s-age-sex
  
  ### Test whether the validity of list of var is checked: #####
  expect_error(FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                              var.smoothTerms = c("wrong_name"),
                              n_cores = 2, pbar = FALSE))
  # duplicated
  temp <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                         var.parametricTerms = c(var.parametricTerms, "p.value"),
                         n_cores = 2, pbar = FALSE)
  expect_equal(temp, mygam_default)
  
  
  ### different arguments in GAM #####
  # s(k=?):
  mygam_k4 <- FixelArray.gam(FD ~ s(age, k=4) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                             n_cores = 2, pbar = FALSE)
  expect_false(all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_k4 %>% dplyr::select("s_age.statistic"))
               %>% isTRUE())    # should be different when k is different
  
  # s(fx=TRUE) vs default (fx=FALSE): 1) different stat; 2) edf is saved or not
  # TODO here.......
  
  # s(bs=?)
  mygam_bsCR <- FixelArray.gam(FD ~ s(age, bs="cr") + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                               n_cores = 2, pbar = FALSE)
  expect_false(all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_bsCR %>% dplyr::select("s_age.statistic"))
               %>% isTRUE()) 
  
  ### Test n_cores, pbar work: ######
  # n_cores = 2: 
  mygam_pbarFalse_ncores2 <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                           n_cores = 2, pbar = FALSE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarFalse_ncores2)
  
  # pbar:
  mygam_pbarTrue_ncores2 <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                  n_cores = 2, pbar = TRUE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarTrue_ncores2)
  
  mygam_pbarTrue_ncores1 <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                                           n_cores = 1, pbar = TRUE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarTrue_ncores1)
  
  
  ### Test: p.value correction: #####
  mygam_fdr <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                 correct.p.value.parametricTerms = c("fdr"),
                 n_cores = 2, pbar = FALSE) 
  expect_equal(mygam_fdr$sexM.p.value.fdr,
               mygam_fdr$sexM.p.value %>% p.adjust("fdr"))
  
  mygam_bonferroni <- FixelArray.gam(FD ~ s(age) + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = fixel.subset,
                              correct.p.value.smoothTerms = c("bonferroni"),
                              n_cores = 2, pbar = FALSE) 
  expect_equal(mygam_bonferroni$s_age.p.value.bonferroni,
               mygam_bonferroni$s_age.p.value %>% p.adjust("bonferroni"))
  
  
  ### Test: eff.size
  # one term of interest:
  
  # more than one term of interest:
  
  # invalid request - see my checker in FixelArray.gam()
})

