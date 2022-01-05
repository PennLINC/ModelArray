test_that("ModelArray.lm() works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
 
  modelarray <- ModelArray(h5_path,
                    "scalar_types" = c("FD"),
                    analysis_names = c("my_analysis"))
  
  # h5_path <- paste0(system.file(package = "ModelArray"),
  #                   "inst/extdata/","n50_fixels.h5")
  # modelarray <- ModelArray(h5_path,
  # scalar_types = c("FD"))
  
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
  # csv_path <- paste0(system.file(package = "ModelArray"),
  #                    "inst/extdata/","n50_cohort.csv")
  
  phenotypes <- read.csv(csv_path)
  scalar_name <- "FD"
  var.terms <- c("estimate", "p.value")   # list of columns to keep  | , "std.error","statistic"
  var.terms.full <- c("estimate", "p.value", "std.error","statistic")
  var.model <- c("r.squared", "p.value", "AIC")
  
  
  ### basic check #####
  mylm <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                        var.terms = var.terms,
                        var.model = var.model,
                        n_cores = 1, pbar=FALSE)
  
  expect_equal(mylm$element_id, 0:99)   # check output$element_id 
  expect_true(is.data.frame(mylm))  # should be data.frame
  expect_equal(as.numeric(dim(mylm)), c(100,1+2*length(var.terms)+length(var.model))) # check shape
  
  mylm_default <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                n_cores = 2, pbar=FALSE)   # default full.outputs and var.*
  expect_equal(as.numeric(dim(mylm_default)), c(100,1+2*3+2)) # check shape  
  
  mylm_fullOutputs <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                    full.outputs = TRUE,   # default: FALSE
                                n_cores = 2, pbar=FALSE)   
  expect_equal(as.numeric(dim(mylm_fullOutputs)), c(100,21))
  
  
  mylm_age_sex <- ModelArray.lm(FD ~ age + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                var.terms = var.terms,
                                var.model = var.model,
                                n_cores = 2, pbar=FALSE)
  expect_equal(mylm_age_sex$element_id, 0:99)   # check output$element_id 
  expect_equal(as.numeric(dim(mylm_age_sex)), c(100,1+3*length(var.terms)+length(var.model))) 
  
  expect_false(all.equal(mylm %>% dplyr::select("age.estimate"), 
                         mylm_age_sex %>% dplyr::select("age.estimate"))
               %>% isTRUE())  # expect not identical between two models

  
  ## Test whether the validity of list of var is checked:
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                             var.terms = c("estimator"),    # misspelling
                             var.model = var.model, 
                             n_cores = 2, pbar=FALSE))
  temp <- (ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                             var.terms = var.terms,
                             var.model = c(var.model, "AIC"),   # duplicated, should be handled
                             n_cores = 2, pbar=FALSE))
  expect_equal(mylm, temp)
  
  
  ### Test n_cores, pbar work #####
  # n_cores=2:
  mylm_ncores2 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                var.terms = var.terms,
                                var.model = var.model, 
                                n_cores = 2, pbar=FALSE)
  expect_equal(mylm, mylm_ncores2)
  # pbar=TRUE & n_cores=1
  mylm_pbarTRUE_ncores1 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                         var.terms = var.terms,
                                         var.model = var.model,
                                         n_cores = 1, pbar=TRUE)
  expect_equal(mylm, mylm_pbarTRUE_ncores1)
  # pbar=TRUE & n_cores=2
  mylm_pbarTRUE_ncores2 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                         var.terms = var.terms,
                                         var.model = var.model,
                                         n_cores = 2, pbar=TRUE)
  expect_equal(mylm, mylm_pbarTRUE_ncores2)
  
  ## Different output statistics #####
  mylm_noTermsOutput <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                              var.terms = c(),
                              var.model = var.model,
                              n_cores = 1, pbar=FALSE)
  expect_equal(as.numeric(dim(mylm_noTermsOutput)), c(100,1+length(var.model))) # check shape
  
  mylm_noModelOutput <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                      var.terms = var.terms,
                                      var.model = c(),
                                      n_cores = 1, pbar=FALSE)
  expect_equal(as.numeric(dim(mylm_noModelOutput)), c(100,1+2*length(var.terms))) # check shape
  
  ## Whether to correct p.values:   #####
  # terms:
  mylm_corr_pvalues_1 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                       var.terms = var.terms, var.model = var.model,
                                       correct.p.value.terms = c("fdr","bonferroni"),
                                       n_cores = 2, pbar=FALSE)
  
  expect_equal(mylm_corr_pvalues_1$age.p.value.fdr,
               mylm_corr_pvalues_1$age.p.value %>% stats::p.adjust("fdr"))
  expect_equal(mylm_corr_pvalues_1$age.p.value.bonferroni,
               mylm_corr_pvalues_1$age.p.value %>% stats::p.adjust("bonferroni"))
  
  # model:
  mylm_corr_pvalues_2 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                       var.terms = var.terms, var.model = var.model,
                                       correct.p.value.model = c("fdr","bonferroni"),
                                       n_cores = 2, pbar=FALSE)
  
  expect_equal(mylm_corr_pvalues_2$model.p.value.fdr,
               mylm_corr_pvalues_2$model.p.value %>% stats::p.adjust("fdr"))
  expect_equal(mylm_corr_pvalues_2$model.p.value.bonferroni,
               mylm_corr_pvalues_2$model.p.value %>% stats::p.adjust("bonferroni"))
  
  
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                             var.terms = var.terms, var.model = var.model,
                             correct.p.value.terms = c("fdr_wrong","bonferroni"),   # wrong name
                             n_cores = 2, pbar=FALSE))
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                             var.terms = var.terms, var.model = var.model,
                             correct.p.value.model = c("fdr_wrong","bonferroni"),   # wrong name
                             n_cores = 2, pbar=FALSE))
  
  expect_warning( ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                var.terms = c("estimate"), var.model = var.model,  # did not provide p.value
                                correct.p.value.terms = c("fdr","bonferroni"),
                                n_cores = 2, pbar=FALSE))
  expect_warning( ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                var.terms = var.terms, var.model = c("AIC"),  # did not provide p.value
                                correct.p.value.model = c("fdr","bonferroni"),
                                n_cores = 2, pbar=FALSE))
  
  ## How about other variables as covariate? factorA is literally correlated with age; factorB is another random variable
  # factor A is fully correlated with age, expecting testing results are NA:
  mylm_age_factorA <- ModelArray.lm(FD ~ age + factorA, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                    var.terms = var.terms.full,
                                    var.model = var.model, 
                                    n_cores = 2, pbar=FALSE)
  # temp <- mylm_age_factorA %>% dplyr::filter(term=="factorA")%>% dplyr::select(-term)  # only extracting column "factorA"
  expect_equal(mylm_age_factorA$element_id, c(0:99))    # test that $element_id is 0:99
  
  expect_true(all(is.na(mylm_age_factorA$factorA.estimate)))  # anything of factorA should be NA
  expect_true(all(is.na(mylm_age_factorA$factorA.p.value)))
  expect_true(all(is.na(mylm_age_factorA$factorA.std.error)))
  expect_true(all(is.na(mylm_age_factorA$factorA.statistic)))

  
  # factor B is not correlated with age, so not expecting NA:
  mylm_age_factorB <- ModelArray.lm(FD ~ age + factorB, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                    var.terms = var.terms.full,
                                    var.model = var.model, 
                                    n_cores = 2, pbar=FALSE)
  # temp <- mylm_age_factorB %>% dplyr::filter(term=="factorB")%>% dplyr::select(-term) %>% dplyr::select(-element_id)  
  expect_false(all(is.na(mylm_age_factorB$factorB.estimate)))
  expect_false(all(is.na(mylm_age_factorB$factorB.p.value)))
  expect_false(all(is.na(mylm_age_factorB$factorB.std.error)))
  expect_false(all(is.na(mylm_age_factorB$factorB.statistic)))
  
  ## Different optional arguments of lm: in order to test that the additional arguments have really been passed into the lm: #####
  # test "na.action" with inputs with NA
  phenotypes_wNA <- phenotypes
  phenotypes_wNA$age[1] <- NA
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = 1:100, 
                             var.terms = var.terms, var.model = var.model, n_cores = 1, pbar=FALSE,
                             na.action="na.fail"))  # expect error of "missing values in object". If na.action was not passed into lm, there will not be error
    # with different n_cores: (n_cores = either 1 or 2: additional arguments of lm have been passed into)
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = 1:100, 
                               var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=FALSE,
                               na.action="na.fail"))  # NOTE: after updating ModelArray.lm with one row, specifying column names to keep, expect error (instead of warning for each core with error, expect_warning) 
    # with different pbar: (TRUE or FALSE: additional arguments of lm have been passed into)
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = 1:100, 
                             var.terms = var.terms, var.model = var.model, n_cores = 1, pbar=TRUE,
                             na.action="na.fail"))  
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = 1:100, 
                               var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=TRUE,
                               na.action="na.fail"))   # after updating ModelArray.lm with one row, specifying column names to keep, expect error (instead of warning for each core with error, expect_warning) 
  
  #mylm_phenotypes_naActionDefault <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = 1:100, n_cores = 2, pbar=FALSE)
  mylm_phenotypes_naActionOmit <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = 1:100, 
                                                var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=FALSE, 
                                                na.action="na.omit")
  # there should be differences in results, after changing one value to NA:
  expect_false(all.equal(mylm                         %>% dplyr::select(age.estimate),
                         mylm_phenotypes_naActionOmit %>% dplyr::select(age.estimate)) 
               %>% isTRUE())
  

  # check if "weights" have been successfully passed into lm:
  
  
  mylm_weights1 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                        var.terms = var.terms, var.model = var.model, 
                        pbar=FALSE, n_cores = 2, 
                        weights = rep(1,nrow(phenotypes)) )   # length(phenotypes$subject_id)   # weights = rep(1,nrow(phenotypes))
  expect_equal(mylm, mylm_weights1)
  
  
  set.seed(5)
  mylm_weightsRnorm <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = 1:100, 
                                     var.terms = var.terms,
                                     var.model = var.model,
                                     n_cores = 2, pbar=FALSE, 
                                     weights = abs(rnorm(nrow(phenotypes))) )  
  expect_false(all.equal(mylm, mylm_weightsRnorm) %>% isTRUE() )
  
  
  # NOTE: we can add more tests regarding other lm's arguments
  
  ### source file list sanity check #####
  # not the same length:
  phenotypes_wrong1 <- phenotypes[-c(1),]
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wrong1, scalar = scalar_name, element.subset = 1:100, 
                             n_cores = 2, pbar=FALSE))
  
  # not the same order --> will be handled!
  phenotypes_swap <- phenotypes
  temp <- phenotypes_swap[2,]   # swap row1 and row2
  phenotypes_swap[2,] <- phenotypes_swap[1,]
  phenotypes_swap[1,] <- temp
  mylm_testSwap <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_swap, scalar = scalar_name, element.subset = 1:100, 
                             n_cores = 2, pbar=FALSE)
  expect_equal(mylm_testSwap, mylm_default)

  # change one row --> cannot be matched --> error:
  phenotypes_wrong2 <- phenotypes
  phenotypes_wrong2[1,"source_file"] <- "wrong_file"
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wrong2, scalar = scalar_name, element.subset = 1:100, 
                             n_cores = 2, pbar=FALSE))
  
  rhdf5::h5closeAll()
  
})


