test_that("FixelArray's lm works as expected", {
  # h5_path <- system.file("extdata", "n50_fixels.h5", package = "FixelArray")   # TODO: ask Tinashe
  
  # fa <- FixelArray(h5_path, 
  #                  "scalar_types" = c("FD"), 
  #                  analysis_names = c("my_analysis"))
  
  h5_path <- paste0(system.file(package = "FixelArray"),
                    "inst/extdata/","n50_fixels.h5")
  fa <- FixelArray(h5_path,
                   scalar_types = c("FD"))
  
  # csv_path <- system.file("extdata", "n50_cohort.csv", package = "FixelArray")   # TODO: ask Tinashe
  csv_path <- paste0(system.file(package = "FixelArray"),
                     "inst/extdata/","n50_cohort.csv")
  
  phenotypes <- read.csv(csv_path)
  scalar_name <- "FD"
  var.terms <- c("estimate", "p.value")   # list of columns to keep  | , "std.error","statistic"
  var.terms.full <- c("estimate", "p.value", "std.error","statistic")
  var.model <- c("r.squared", "p.value", "AIC")
  mylm <- FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = 1:100, 
                        var.terms = var.terms,
                        var.model = var.model,
                        n_cores = 1, pbar=FALSE)
  
  expect_equal(mylm$fixel_id, 0:99)   # check output$fixel_id 
  expect_true(is.data.frame(mylm))  # should be data.frame
  expect_equal(as.numeric(dim(mylm)), c(100,1+2*length(var.terms)+length(var.model))) # check shape
  
  mylm_age_sex <- FixelArray.lm(FD ~ age + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = 1:100, 
                                var.terms = var.terms,
                                var.model = var.model,
                                n_cores = 1, pbar=FALSE)
  expect_equal(mylm_age_sex$fixel_id, 0:99)   # check output$fixel_id 
  expect_equal(as.numeric(dim(mylm_age_sex)), c(100,1+3*length(var.terms)+length(var.model))) 
  
  expect_false(all.equal(mylm %>% dplyr::select("age.estimate"), 
                         mylm_age_sex %>% dplyr::select("age.estimate"))
               %>% isTRUE())  # expect not identical between two models

  
  ## Test n_cores, pbar work:
  # n_cores=2:
  mylm_ncores2 <- FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = 1:100, 
                                var.terms = var.terms,
                                var.model = var.model, 
                                n_cores = 2, pbar=FALSE)
  expect_equal(mylm, mylm_ncores2)
  # pbar=TRUE & n_cores=1
  mylm_pbarTRUE_ncores1 <- FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = 1:100, 
                                         var.terms = var.terms,
                                         var.model = var.model,
                                         n_cores = 1, pbar=TRUE)
  expect_equal(mylm, mylm_pbarTRUE_ncores1)
  # pbar=TRUE & n_cores=2
  mylm_pbarTRUE_ncores2 <- FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = 1:100, 
                                         var.terms = var.terms,
                                         var.model = var.model,
                                         n_cores = 2, pbar=TRUE)
  expect_equal(mylm, mylm_pbarTRUE_ncores2)
  
  
  ## How about other variables as covariate? factorA is literally correlated with age; factorB is another random variable
  # factor A is fully correlated with age, expecting testing results are NA:
  mylm_age_factorA <- FixelArray.lm(FD ~ age + factorA, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = 1:100, 
                                    var.terms = var.terms.full,
                                    var.model = var.model, 
                                    n_cores = 2, pbar=FALSE)
  # temp <- mylm_age_factorA %>% dplyr::filter(term=="factorA")%>% dplyr::select(-term)  # only extracting column "factorA"
  expect_equal(mylm_age_factorA$fixel_id, c(0:99))    # test that $fixel_id is 0:99
  
  expect_true(all(is.na(mylm_age_factorA$factorA.estimate)))  # anything of factorA should be NA
  expect_true(all(is.na(mylm_age_factorA$factorA.p.value)))
  expect_true(all(is.na(mylm_age_factorA$factorA.std.error)))
  expect_true(all(is.na(mylm_age_factorA$factorA.statistic)))

  
  # factor B is not correlated with age, so not expecting NA:
  mylm_age_factorB <- FixelArray.lm(FD ~ age + factorB, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = 1:100, 
                                    var.terms = var.terms.full,
                                    var.model = var.model, 
                                    n_cores = 2, pbar=FALSE)
  # temp <- mylm_age_factorB %>% dplyr::filter(term=="factorB")%>% dplyr::select(-term) %>% dplyr::select(-fixel_id)  
  expect_false(all(is.na(mylm_age_factorB$factorB.estimate)))
  expect_false(all(is.na(mylm_age_factorB$factorB.p.value)))
  expect_false(all(is.na(mylm_age_factorB$factorB.std.error)))
  expect_false(all(is.na(mylm_age_factorB$factorB.statistic)))
  
  ## Different optional arguments of lm: in order to test that the additional arguments have really been passed into the lm:
  # test "na.action" with inputs with NA
  phenotypes_wNA <- phenotypes
  phenotypes_wNA$age[1] <- NA
  expect_error(FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes_wNA, scalar = scalar_name, fixel.subset = 1:100, 
                             var.terms = var.terms, var.model = var.model, n_cores = 1, pbar=FALSE,
                             na.action="na.fail"))  # expect error of "missing values in object". If na.action was not passed into lm, there will not be error
    # with different n_cores: (n_cores = either 1 or 2: additional arguments of lm have been passed into)
  expect_error(FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes_wNA, scalar = scalar_name, fixel.subset = 1:100, 
                               var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=FALSE,
                               na.action="na.fail"))  # NOTE: after updating FixelArray.lm with one row, specifying column names to keep, expect error (instead of warning for each core with error, expect_warning) 
    # with different pbar: (TRUE or FALSE: additional arguments of lm have been passed into)
  expect_error(FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes_wNA, scalar = scalar_name, fixel.subset = 1:100, 
                             var.terms = var.terms, var.model = var.model, n_cores = 1, pbar=TRUE,
                             na.action="na.fail"))  
  expect_error(FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes_wNA, scalar = scalar_name, fixel.subset = 1:100, 
                               var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=TRUE,
                               na.action="na.fail"))   # after updating FixelArray.lm with one row, specifying column names to keep, expect error (instead of warning for each core with error, expect_warning) 
  
  #mylm_phenotypes_naActionDefault <- FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes_wNA, scalar = scalar_name, fixel.subset = 1:100, n_cores = 2, pbar=FALSE)
  mylm_phenotypes_naActionOmit <- FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes_wNA, scalar = scalar_name, fixel.subset = 1:100, 
                                                var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=FALSE, 
                                                na.action="na.omit")
  # there should be differences in results, after changing one value to NA:
  expect_false(all.equal(mylm                         %>% dplyr::select(age.estimate),
                         mylm_phenotypes_naActionOmit %>% dplyr::select(age.estimate)) 
               %>% isTRUE())
  

  # NOTE: we can add more tests regarding other lm's arguments
  
  ## TODO: check if "weights" have been successfully passed into lm:
  # FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes, scalar = scalar_name, fixel.subset = 1:100, n_cores = 1, pbar=FALSE,  
  #            weights = rep(1,length(phenotypes$subject_id)) )
  
  
  ## TODO: check var.terms = c(); var.model =c() --> any error
  
  h5closeAll()
  
})


