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
  
  default.var.terms = c("estimate", "statistic", "p.value")
  default.var.model = c("adj.r.squared", "p.value")
  var.terms.full = c("estimate","std.error","statistic","p.value")
  var.model.full = c("r.squared", "adj.r.squared", "sigma", "statistic", "p.value", "df", "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs")
  
  element.subset = 1:10
  
  ### generate + load expected results #####
  idx.fixel.lm <- 1
  num.set.seed <- 5  # this will be used when generating random numbers (in expected results and below)
  # generate the expected results, and get `expected.results`, a list of the expected results
  expected.results <- helper_generate_expect_lm(fn.phenotypes = csv_path, 
                                                fn.h5 = h5_path,
                                                idx.fixel.lm = idx.fixel.lm,
                                                num.set.seed = num.set.seed)
  
  #' @param idx.fixel starts from 1
  #' 
  compare_expected_results <- function(actual, expected, idx.fixel = idx.fixel.lm) {
    ## check if each p.value correction is correct:
    # before removing other rows in "actual":
    col.names <- colnames(actual)
    for (name.p.adjust in p.adjust.methods) {   # iterate over different correction methods
      #message(name.p.adjust)
      if (stringr::str_detect(col.names, name.p.adjust) %>% any()) {   # if it contains this name.p.adjust
        #message("detected!")
        list.idx <- stringr::str_which(col.names, name.p.adjust)  # the index of the columns that contains this name.p.adjust
        for (idx in list.idx) {   # for each corrected term/model's p.value
          #message(toString(idx))
          thecolname <- col.names[idx]
          thecolname.pvalue <- stringr::str_remove_all(thecolname, 
                                                        paste0(".",name.p.adjust))  # remove the name.p.adjust to get the p.value's name
          expect_equal(actual[[thecolname.pvalue]] %>% stats::p.adjust(method = name.p.adjust),
                       actual[[thecolname]])
          
        }
      }
    }
    
    # now, remove any fdr, etc corrections:
    actual <- actual %>% select(-ends_with(p.adjust.methods))
    col.names <- colnames(actual)
    
    # select the corresponding row (which has the expected value):
    actual <- actual[actual$element_id ==idx.fixel - 1, ]
    
    flag.belong <- col.names %in% colnames(expected) %>% all()
    if (flag.belong == FALSE) {
      stop("not all columns in actual data.frame are in expected data.frame")
    }
    
    ## test if actual = expected values:
    expect_equal(actual,
                 expected %>% select(col.names))
  }
  
  ### basic check #####
  mylm <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                        var.terms = var.terms,
                        var.model = var.model,
                        n_cores = 1, pbar=FALSE)
  
  compare_expected_results(mylm, expected.results$age)
  expect_equal(mylm$element_id, 0:(length(element.subset)-1))   # check output$element_id 
  expect_true(is.data.frame(mylm))  # should be data.frame
  expect_equal(as.numeric(dim(mylm)), c(length(element.subset),   # check shape
                                        1+
                                          2*(length(var.terms)+1)+   # 2*: intercept + age; +1 in length: fdr is default
                                          1*(length(var.model)+1)))  # 1*: lm model; +1 in length: fdr is default
  expect_true(c("age.p.value.fdr", "Intercept.p.value.fdr", "model.p.value.fdr") %in% colnames(mylm) %>% all() )  # fdr is saved by default
  
  mylm_default <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                n_cores = 2, pbar=FALSE)   # default full.outputs and var.*
  compare_expected_results(mylm_default, expected.results$age)
  expect_equal(as.numeric(dim(mylm_default)), c(length(element.subset),     # check shape  
                                                1+
                                                  2*(length(default.var.terms)+1)+
                                                  1*(length(default.var.model)+1))) 
  
  mylm_fullOutputs <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                    full.outputs = TRUE,   # default: FALSE
                                n_cores = 2, pbar=FALSE)   
  compare_expected_results(mylm_fullOutputs, expected.results$age)
  expect_equal(as.numeric(dim(mylm_fullOutputs)), c(length(element.subset),
                                                    1+
                                                      2*(length(var.terms.full) + 1) + 
                                                      1*(length(var.model.full) + 1)))
  expect_true(c("age.p.value.fdr", "Intercept.p.value.fdr", "model.p.value.fdr") %in% colnames(mylm) %>% all())  # fdr is saved by default
  
  mylm_age_sex <- ModelArray.lm(FD ~ age + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                var.terms = var.terms,
                                var.model = var.model,
                                n_cores = 2, pbar=FALSE)
  compare_expected_results(mylm_age_sex, expected.results$age_sex)
  expect_equal(mylm_age_sex$element_id, 0:(length(element.subset)-1))   # check output$element_id 
  expect_equal(as.numeric(dim(mylm_age_sex)), c(length(element.subset),
                                                1+
                                                  3*(length(var.terms)+1)+
                                                  1*(length(var.model)+1))) 
  
  expect_false(all.equal(mylm %>% dplyr::select("age.estimate"), 
                         mylm_age_sex %>% dplyr::select("age.estimate"))
               %>% isTRUE())  # expect not identical between two models

  
  ## Test whether the validity of list of var is checked:
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                             var.terms = c("estimator"),    # misspelling
                             var.model = var.model, 
                             n_cores = 2, pbar=FALSE))
  temp <- (ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                             var.terms = var.terms,
                             var.model = c(var.model, "AIC"),   # duplicated, should be handled
                             n_cores = 2, pbar=FALSE))
  expect_equal(mylm, temp)
  
  
  ### Test n_cores, pbar work #####
  # n_cores=2:
  mylm_ncores2 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                var.terms = var.terms,
                                var.model = var.model, 
                                n_cores = 2, pbar=FALSE)
  expect_equal(mylm, mylm_ncores2)
  # pbar=TRUE & n_cores=1
  mylm_pbarTRUE_ncores1 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                         var.terms = var.terms,
                                         var.model = var.model,
                                         n_cores = 1, pbar=TRUE)
  expect_equal(mylm, mylm_pbarTRUE_ncores1)
  # pbar=TRUE & n_cores=2
  mylm_pbarTRUE_ncores2 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                         var.terms = var.terms,
                                         var.model = var.model,
                                         n_cores = 2, pbar=TRUE)
  expect_equal(mylm, mylm_pbarTRUE_ncores2)
  
  ## Different output statistics #####
  expect_warning(mylm_noTermsOutput <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                      var.terms = c(), var.model = var.model,
                                      n_cores = 1, pbar=FALSE))  # expect warning because there will be a warning from ModelArray.lm(): p.value was not included in var.terms, so not to perform its p.value corrections
  compare_expected_results(mylm_noTermsOutput, expected.results$age)
  expect_equal(as.numeric(dim(mylm_noTermsOutput)), c(length(element.subset),
                                                      1+
                                                        length(var.model)+1)) # check shape
    # there shouldn't be any NA in this output: (otherwise, there is a bug when combining empty tibble in analyseOneElement.lm())
  expect_true( (!mylm_noTermsOutput %>% is.na()) %>% all() )   # if there is any NA in the output data.frame, the result is FALSE
  
  expect_warning(mylm_noModelOutput <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                      var.terms = var.terms, var.model = c(),
                                      n_cores = 1, pbar=FALSE))   # expect warning because there will be a warning from ModelArray.lm(): p.value was not included in var.model, so not to perform its p.value corrections 
  compare_expected_results(mylm_noModelOutput, expected.results$age)
  expect_equal(as.numeric(dim(mylm_noModelOutput)), c(length(element.subset),
                                                      1+
                                                        2*(length(var.terms)+1))) # check shape
  expect_true( (!mylm_noModelOutput %>% is.na()) %>% all() )   
  
  # there will be error if both var.* are empty:
  expect_error(mylm_noModelOutput <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                                   var.terms = c(), var.model = c(),
                                                   n_cores = 1, pbar=FALSE))
  
  ## Whether to correct p.values:   #####
  # terms:
  mylm_corr_pvalues_1 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                       var.terms = var.terms, var.model = var.model,
                                       correct.p.value.terms = c("fdr","bonferroni"),
                                       n_cores = 2, pbar=FALSE)
  compare_expected_results(mylm_corr_pvalues_1, expected.results$age)   # this includes check *.p.value.fdr and *.p.value.bonferroni values are correct
  # check if requested p value corrections exist:
  expect_true(c("age.p.value.bonferroni","Intercept.p.value.bonferroni")
              %in% colnames(mylm_corr_pvalues_1)
              %>% all())
  
  # model:
  mylm_corr_pvalues_2 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                       var.terms = var.terms, var.model = var.model,
                                       correct.p.value.model = c("fdr","bonferroni"),
                                       n_cores = 2, pbar=FALSE)
  compare_expected_results(mylm_corr_pvalues_2, expected.results$age)
  # check if requested p value corrections exist:
  expect_true(c("model.p.value.bonferroni")
              %in% colnames(mylm_corr_pvalues_2)
              %>% all())
  
  # terms + only bonferroni, no fdr:
  mylm_corr_pvalues_3 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                       var.terms = var.terms, var.model = var.model,
                                       correct.p.value.terms = c("bonferroni"),
                                       n_cores = 2, pbar=FALSE)
  compare_expected_results(mylm_corr_pvalues_3, expected.results$age)
  expect_true(c("age.p.value.bonferroni", "Intercept.p.value.bonferroni",
                "model.p.value.fdr") %in% colnames(mylm_corr_pvalues_3) %>% all())   # correct.p.value.terms: bonferroni is in, but not fdr
  expect_true( !( c("age.p.value.fdr", "Intercept.p.value.fdr") %in% colnames(mylm_corr_pvalues_3) ) %>% all() )  # all false
  
  # model + only bonferroni, no fdr:
  mylm_corr_pvalues_4 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                       var.terms = var.terms, var.model = var.model,
                                       correct.p.value.model = c("bonferroni"),
                                       n_cores = 2, pbar=FALSE)
  compare_expected_results(mylm_corr_pvalues_4, expected.results$age)
  expect_true(c("age.p.value.fdr",
                "model.p.value.bonferroni") %in% colnames(mylm_corr_pvalues_4) %>% all())   # correct.p.value.model: bonferroni is in, but not fdr
  expect_true( !(c("model.p.value.fdr") %in% colnames(mylm_corr_pvalues_4) ) %>% all() )   # all false
  
  
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                             var.terms = var.terms, var.model = var.model,
                             correct.p.value.terms = c("fdr_wrong","bonferroni"),   # wrong name
                             n_cores = 2, pbar=FALSE))
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                             var.terms = var.terms, var.model = var.model,
                             correct.p.value.model = c("fdr_wrong","bonferroni"),   # wrong name
                             n_cores = 2, pbar=FALSE))
  
  expect_warning( temp <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                var.terms = c("estimate"), var.model = var.model,  # did not provide p.value
                                correct.p.value.terms = c("fdr","bonferroni"),
                                n_cores = 2, pbar=FALSE))
  compare_expected_results(temp, expected.results$age)
  
  expect_warning( temp <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                var.terms = var.terms, var.model = c("AIC"),  # did not provide p.value
                                correct.p.value.model = c("fdr","bonferroni"),
                                n_cores = 2, pbar=FALSE))
  compare_expected_results(temp, expected.results$age)
  
  ## How about other variables as covariate? #####
  # factorA is literally correlated with age; factorB is another random variable
  # factor A is fully correlated with age, expecting testing results are NA:
  mylm_age_factorA <- ModelArray.lm(FD ~ age + factorA, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                    var.terms = var.terms.full,
                                    var.model = var.model, 
                                    n_cores = 2, pbar=FALSE)
  compare_expected_results(mylm_age_factorA, expected.results$age_factorA)
  # temp <- mylm_age_factorA %>% dplyr::filter(term=="factorA")%>% dplyr::select(-term)  # only extracting column "factorA"
  expect_equal(mylm_age_factorA$element_id, c(0:(length(element.subset)-1)))    # test that $element_id is 0:(length(element.subset)-1)
  
  expect_true(all(is.na(mylm_age_factorA$factorA.estimate)))  # anything of factorA should be NA
  expect_true(all(is.na(mylm_age_factorA$factorA.p.value)))
  expect_true(all(is.na(mylm_age_factorA$factorA.std.error)))
  expect_true(all(is.na(mylm_age_factorA$factorA.statistic)))

  
  # factor B is not correlated with age, so not expecting NA:
  mylm_age_factorB <- ModelArray.lm(FD ~ age + factorB, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                    var.terms = var.terms.full,
                                    var.model = var.model, 
                                    n_cores = 2, pbar=FALSE)
  compare_expected_results(mylm_age_factorB, expected.results$age_factorB)
  # temp <- mylm_age_factorB %>% dplyr::filter(term=="factorB")%>% dplyr::select(-term) %>% dplyr::select(-element_id)  
  expect_false(all(is.na(mylm_age_factorB$factorB.estimate)))
  expect_false(all(is.na(mylm_age_factorB$factorB.p.value)))
  expect_false(all(is.na(mylm_age_factorB$factorB.std.error)))
  expect_false(all(is.na(mylm_age_factorB$factorB.statistic)))
  
  ## Different optional arguments of lm: in order to test that the additional arguments have really been passed into the lm: #####
  # hanlding inputs with NA:
  phenotypes_wNA <- phenotypes
  phenotypes_wNA$age[1] <- NA
  
  # test default `lm()` method for handling NA:
  mylm_phenotypes_naActionDefault <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, 
                                                   scalar = scalar_name, element.subset = element.subset, n_cores = 2, pbar=FALSE)
  compare_expected_results(mylm_phenotypes_naActionDefault, expected.results[["age_phenotypeswNA"]])
  
  # there should be differences in results, after changing one value to NA:
  expect_false(all.equal(mylm                         %>% dplyr::select(age.estimate),
                         mylm_phenotypes_naActionDefault %>% dplyr::select(age.estimate)) 
               %>% isTRUE())

  # test different (not-default) `na.action` with inputs with NA:
  # error if `na.action = "na.fail"`:
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = element.subset, 
                             var.terms = var.terms, var.model = var.model, n_cores = 1, pbar=FALSE,
                             na.action="na.fail"))  # expect error of "missing values in object". If na.action was not passed into lm, there will not be error
    # with different n_cores: (n_cores = either 1 or 2: additional arguments of lm have been passed into)
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = element.subset, 
                               var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=FALSE,
                               na.action="na.fail"))  # NOTE: after updating ModelArray.lm with one row, specifying column names to keep, expect error (instead of warning for each core with error, expect_warning) 
    # with different pbar: (TRUE or FALSE: additional arguments of lm have been passed into)
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = element.subset, 
                             var.terms = var.terms, var.model = var.model, n_cores = 1, pbar=TRUE,
                             na.action="na.fail"))  
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = element.subset, 
                               var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=TRUE,
                               na.action="na.fail"))   # after updating ModelArray.lm with one row, specifying column names to keep, expect error (instead of warning for each core with error, expect_warning) 
  
  # okay if `na.action = "na.omit"`:
  mylm_phenotypes_naActionOmit <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wNA, scalar = scalar_name, element.subset = element.subset, 
                                                var.terms = var.terms, var.model = var.model, n_cores = 2, pbar=FALSE, 
                                                na.action="na.omit")
  compare_expected_results(mylm_phenotypes_naActionOmit, expected.results[["age_phenotypeswNA_na.action-na.omit"]])
  # there should be differences in results, after changing one value to NA:
  expect_false(all.equal(mylm                         %>% dplyr::select(age.estimate),
                         mylm_phenotypes_naActionOmit %>% dplyr::select(age.estimate)) 
               %>% isTRUE())
  # default option to handle NA should be the same as using `na.omit`:
  expect_true(all.equal(mylm_phenotypes_naActionDefault %>% dplyr::select(age.estimate),
                        mylm_phenotypes_naActionOmit %>% dplyr::select(age.estimate)) 
               %>% isTRUE())

  
  # check if "weights" have been successfully passed into lm:
  mylm_weights1 <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                        var.terms = var.terms, var.model = var.model, 
                        pbar=FALSE, n_cores = 2, 
                        weights = rep(1,nrow(phenotypes)) )   # length(phenotypes$subject_id)   # weights = rep(1,nrow(phenotypes))
  expect_equal(mylm, mylm_weights1)
  
  
  set.seed(num.set.seed)
  mylm_weightsRnorm <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset, 
                                     var.terms = var.terms,
                                     var.model = var.model,
                                     n_cores = 2, pbar=FALSE, 
                                     weights = abs(rnorm(nrow(phenotypes))) )  
  compare_expected_results(mylm_weightsRnorm, expected.results[["age_weights-random"]])
  expect_false(all.equal(mylm, mylm_weightsRnorm) %>% isTRUE() )
  
  
  # NOTE: we can add more tests regarding other lm's arguments
  
  ### source file list sanity check #####
  # not the same length:
  phenotypes_wrong1 <- phenotypes[-c(1),]
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wrong1, scalar = scalar_name, element.subset = element.subset, 
                             n_cores = 2, pbar=FALSE))
  
  # not the same order --> will be handled!
  phenotypes_swap <- phenotypes
  temp <- phenotypes_swap[2,]   # swap row1 and row2
  phenotypes_swap[2,] <- phenotypes_swap[1,]
  phenotypes_swap[1,] <- temp
  mylm_testSwap <- ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_swap, scalar = scalar_name, element.subset = element.subset, 
                             n_cores = 2, pbar=FALSE)
  expect_equal(mylm_testSwap, mylm_default)

  # change one row --> cannot be matched --> error:
  phenotypes_wrong2 <- phenotypes
  phenotypes_wrong2[1,"source_file"] <- "wrong_file"
  expect_error(ModelArray.lm(FD ~ age, data = modelarray, phenotypes = phenotypes_wrong2, scalar = scalar_name, element.subset = element.subset, 
                             n_cores = 2, pbar=FALSE))
  
  rhdf5::h5closeAll()
  
})


