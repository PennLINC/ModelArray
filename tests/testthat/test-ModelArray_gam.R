test_that("test that ModelArray.gam() works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
 
  scalar_name <- c("FD")
  modelarray <- ModelArray(h5_path,
                    scalar_types = scalar_name,
                    analysis_names = c("my_analysis"))
  # h5_path <- paste0(system.file(package = "ModelArray"),
  #                   "inst/extdata/","n50_fixels.h5")
  # scalar_name = c("FD")
  # modelarray <- ModelArray(h5_path, scalar_types = scalar_name, analysis_names = c("my_analysis"))
  

  csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
  # csv_path <- paste0(system.file(package = "ModelArray"),
  #                    "inst/extdata/","n50_cohort.csv")
  
  phenotypes <- read.csv(csv_path)
  phenotypes$oSex <- ordered(phenotypes$sex, levels = c("F", "M"))  # ordered factor, "F" as reference group
  phenotypes$sexFactor <- factor(phenotypes$sex, levels = unique(phenotypes$sex))   # factor but not ordered
  
  var.smoothTerms = c("statistic","p.value")
  var.parametricTerms = c("estimate", "statistic", "p.value")
  var.model = c("dev.expl", "adj.r.squared")
  element.subset = 1:10
  
  var.smoothTerms.default <- c("statistic","p.value")
  var.parametricTerms.default <- c("estimate", "statistic", "p.value")
  var.model.default <- c("dev.expl")
  
  var.smoothTerms.full <- c("edf","ref.df","statistic","p.value")
  var.parametricTerms.full <- c("estimate", "std.error","statistic","p.value")
  var.model.full <- c("adj.r.squared","dev.expl", "sp.criterion", "scale",
                      "df", "logLik","AIC", "BIC", "deviance", "df.residual", "nobs")
  
  # calculating partial rsq for modelarray (this is not comprehensive; does not include additional methods; focus is to check if partial rsq is calculated correctly)
  partialRsq_gam_modelarray <- function(i_element, modelarray, phenotypes, scalar, full.formula, reduced.formula) {
    values <- scalars(modelarray)[[scalar]][i_element,]
    dat <- phenotypes
    dat[[scalar]] <- values

    fullmodel <- mgcv::gam(full.formula, data=dat)
    redmodel <- mgcv::gam(reduced.formula, data=dat)

    # calculating SSE: used observed y (i.e. excluding observations with NA), and fitted values, directly from model object
    
    sse.full <- sum( (fullmodel$y - fullmodel$fitted.values)^2 )
    #message(sse.full)
    sse.red <- sum( (redmodel$y - redmodel$fitted.values)^2 )
    #message(sse.red)
    #message("")
    
    partialRsq <- (sse.red - sse.full) / sse.red
    
    return(partialRsq)
    # toReturn <- list(partialRsq = partialRsq,
    #                  sse.full = sse.full,
    #                  sse.red = sse.red)
    # return(toReturn)
  }
  
  wrapper_partialRsq_gam_modelarray <- function(element.subset, modelarray, phenotypes, scalar, full.formula, reduced.formula) {
    partial.rsq.list <- lapply(element.subset,   # a list of i_element
                                  partialRsq_gam_modelarray,  # the function
                                  modelarray, phenotypes, scalar, full.formula, reduced.formula)
    partial.rsq.vec <- do.call(rbind, partial.rsq.list)
    partial.rsq.vec <- as.numeric(partial.rsq.vec)
    return(partial.rsq.vec)

  }
  
  ### basic checks #####
  mygam <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                          var.smoothTerms = var.smoothTerms,
                          var.parametricTerms = var.parametricTerms,
                          var.model = var.model,
                          n_cores = 1, pbar = FALSE)

  expect_equal(mygam$element_id, 0:(max(element.subset)-1))  # check output$element_id
  expect_true(is.data.frame(mygam))   # should be data.frame
  expect_equal(as.numeric(dim(mygam)), c(length(element.subset), 
                                         1+
                                           1*(length(var.smoothTerms)+1) + 
                                           2*(length(var.parametricTerms)+1) + 
                                           length(var.model)))  # check the shape

  mygam_default <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                  n_cores = 1, pbar = FALSE)   # default stat outputs
  expect_equal(as.numeric(dim(mygam_default)), c(length(element.subset),
                                                 1+
                                                   1*(length(var.smoothTerms.default)+1) + 
                                                   2*(length(var.parametricTerms.default)+1) + 
                                                   length(var.model.default)))

  mygam_fullOutputs <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                      full.outputs = TRUE,    # overwrites the var* arguments below
                                      var.smoothTerms = var.smoothTerms,
                                      var.parametricTerms = var.parametricTerms,
                                      var.model = var.model,
                                      n_cores = 1, pbar = FALSE)
  expect_equal(as.numeric(dim(mygam_fullOutputs)), c(length(element.subset),
                                                     1+
                                                       1*(length(var.smoothTerms.full)+1) + 
                                                       2*(length(var.parametricTerms.full)+1)+
                                                       length(var.model.full)))
  
  ### when there is no term or stat output #####
  # if there is no explicit parametric term (although there will be term "Intercept") or parametric stat:
  mygam_noExplicitParamTerm <- ModelArray.gam(FD ~ s(age), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                  n_cores = 2, pbar = FALSE)  
  expect_warning(mygam_noParamStat <- ModelArray.gam(FD ~ s(age)+sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                  var.parametricTerms = c(),     # there will be warning: p.value was not included in var.parametricTerms, so not to perform its p.value corrections
                                  n_cores = 2, pbar = FALSE))
  # if there is no smooth term or stat:
  mygam_noSmoothTerm <- ModelArray.gam(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                  n_cores = 2, pbar = FALSE)   # should without error # expect_output() to test output "there is no smooth term in the requested formula" does not work
  expect_true(c("age.statistic", "age.estimate") %in% colnames(mygam_noSmoothTerm) %>% all())
  
  expect_warning(mygam_noSmoothStat <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                       var.smoothTerms = c(),   # warning: p.value was not included in var.smoothTerms, so not to perform its p.value corrections
                                       n_cores = 2, pbar = FALSE))
  expect_false("age.statistic" %in% colnames(mygam_noSmoothStat))
  
  # if there is only one var.* not empty: (test with request of reduced model - where there will calls as: c(), c(), c("adj.r.squared"))
  w <- capture_warnings(mygam_onlyModelStat <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                      var.smoothTerms=c(), var.parametricTerms=c(), var.model = c("dev.expl"),
                                      changed.rsq.term.index=c(1),
                                      n_cores = 2, pbar = FALSE) )    # warning: p.value was not included in var.smoothTerms OR var.parametricTerms, so not to perform its p.value corrections
    expect_match(w, "p.value was not included in var.smoothTerms", all=FALSE)   # all=FALSE: only one element instead of all elements in actual value need to match)
    expect_match(w, "p.value was not included in var.parametricTerms", all=FALSE)
  w <- capture_warnings(mygam_onlyParametricTermsStat <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                        var.smoothTerms=c(), var.parametricTerms=c("estimate"), var.model = c(),   # warning: p.value was not included in var.xxx
                                        changed.rsq.term.index=c(1),
                                        n_cores = 2, pbar = FALSE))
    expect_match(w, "p.value was not included in var.smoothTerms", all=FALSE)   # all=FALSE: only one element instead of all elements in actual value need to match)
    expect_match(w, "p.value was not included in var.parametricTerms", all=FALSE)
  w <- capture_warnings(mygam_onlySmoothTermsStat <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                        var.smoothTerms=c("statistic"), var.parametricTerms=c(), var.model = c(),
                                        changed.rsq.term.index=c(1),
                                        n_cores = 2, pbar = FALSE) )   # warning: p.value was not included in var.xxx
    expect_match(w, "p.value was not included in var.smoothTerms", all=FALSE)   # all=FALSE: only one element instead of all elements in actual value need to match)
    expect_match(w, "p.value was not included in var.parametricTerms", all=FALSE)
      # there shouldn't be any NA in these outputs: (otherwise, there is a bug when combining empty tibble in analyseOneElement.gam())
  expect_true( (!mygam_onlyModelStat %>% is.na()) %>% all())   # if there is any NA in the output data.frame, the result is FALSE
  expect_true( (!mygam_onlyParametricTermsStat %>% is.na()) %>% all()) 
  expect_true( (!mygam_onlySmoothTermsStat %>% is.na()) %>% all())
  
  # no any model stat:
  expect_error(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                      var.smoothTerms = c(), var.parametricTerms = c(), var.model = c(),
                                      n_cores = 2, pbar = FALSE))
  
  
  ## multiple smooth terms:
  # s(age) + s(factorA)   # cannot use s(sex) as sex are characters (M or F), not values!
  mygam_sAge_sFactorA <- ModelArray.gam(FD ~ s(age) + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                n_cores = 2, pbar = FALSE)
  
  # # s(age + factorA)   # not good example
  # expect_warning(ModelArray.gam(FD ~ s(age + factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
  #                n_cores = 2, pbar = FALSE))  # note: the column name will be s-age instead of s-age-sex
  
  
  
  ### Test whether the validity of list of var is checked: #####
  expect_error(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                              var.smoothTerms = c("wrong_name"),
                              n_cores = 2, pbar = FALSE))
  # duplicated
  temp <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                         var.parametricTerms = c(var.parametricTerms, "p.value"),
                         n_cores = 2, pbar = FALSE)
  expect_equal(temp, mygam_default)
  
  
  
  ### different arguments in GAM #####
  ## different settings in formula:
  # s(k=?):
  mygam_k4 <- ModelArray.gam(FD ~ s(age, k=4) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                             n_cores = 2, pbar = FALSE)
  expect_false(dplyr::all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_k4 %>% dplyr::select("s_age.statistic"))
               %>% isTRUE())    # should be different when k is different
  
  # s(fx=TRUE) vs default (fx=FALSE): 1) different stat; 
  mygam_fxTRUE <- ModelArray.gam(FD ~ s(age, fx=TRUE) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                 n_cores = 2, pbar = FALSE)
  expect_false(dplyr::all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_fxTRUE %>% dplyr::select("s_age.statistic"))
               %>% isTRUE())
  
  # TODO: if force saving edf when fx=FALSE: 2) check if it's saved or not

  # s(bs=?)
  mygam_bsCR <- ModelArray.gam(FD ~ s(age, bs="cr") + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                               n_cores = 2, pbar = FALSE)
  expect_false(dplyr::all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_bsCR %>% dplyr::select("s_age.statistic"))
               %>% isTRUE()) 
  
  ## different settings in mgcv::gam()'s additional arguments: test if the arguments have been passed into analyseOneElement.gam()
  # method: 
  mygam_methodREML <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                     method = "REML",  # default = "GCV.Cp"
                                     n_cores = 2, pbar = FALSE)
  expect_false(dplyr::all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_methodREML %>% dplyr::select("s_age.statistic"))
               %>% isTRUE())
  ### Test n_cores, pbar work: ######
  # n_cores = 2: 
  mygam_pbarFalse_ncores2 <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                           n_cores = 2, pbar = FALSE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarFalse_ncores2)
  
  # pbar:
  mygam_pbarTrue_ncores2 <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                  n_cores = 2, pbar = TRUE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarTrue_ncores2)
  
  mygam_pbarTrue_ncores1 <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                           n_cores = 1, pbar = TRUE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarTrue_ncores1)
  
  
  ### Test: p.value correction: #####
  mygam_parametric_pCorrect <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                 correct.p.value.parametricTerms = c("fdr","bonferroni"),
                 n_cores = 2, pbar = FALSE) 
  expect_equal(mygam_parametric_pCorrect$sexM.p.value.fdr,
               mygam_parametric_pCorrect$sexM.p.value %>% stats::p.adjust("fdr"))
  expect_equal(mygam_parametric_pCorrect$sexM.p.value.bonferroni,
               mygam_parametric_pCorrect$sexM.p.value %>% stats::p.adjust("bonferroni"))
  
  mygam_smooth_pCorrect <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                              correct.p.value.smoothTerms = c("fdr","bonferroni"),
                              n_cores = 2, pbar = FALSE) 
  expect_equal(mygam_smooth_pCorrect$s_age.p.value.fdr,
               mygam_smooth_pCorrect$s_age.p.value %>% stats::p.adjust("fdr"))
  expect_equal(mygam_smooth_pCorrect$s_age.p.value.bonferroni,
               mygam_smooth_pCorrect$s_age.p.value %>% stats::p.adjust("bonferroni"))

  expect_error(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                              correct.p.value.parametricTerms = c("wrong_correct"),
                              n_cores = 2, pbar = FALSE))
  expect_warning(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                var.smoothTerms = c("statistic"),  # no p.value
                                correct.p.value.smoothTerms = c("fdr"),
                                n_cores = 2, pbar = FALSE))
  
  ### Test: changed R-squared (delta.adj.rsq, partial.rsq) #####
  #### one term of interest: reduced model will be FD ~ 1 #####
  # also, to test whether the delta.adj.rsq and partial.rsq are calculated correctly
  # also, to test whether without "," in s(), the column name could be correctly "s_age.delta.adj.rsq"
  mygam_changedrsq_oneSmoothTerm <- ModelArray.gam(FD ~ s(age), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                               var.model = c("dev.expl", "adj.r.squared"),
                                               changed.rsq.term.index = c(1),
                                               n_cores = 2, pbar = FALSE)
  mygam_intercept <- ModelArray.gam(FD ~ 1, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                    var.model = c("adj.r.squared"),
                                    n_cores = 2, pbar = FALSE)
  # delta.adj.rsq:
  expect_equal(mygam_changedrsq_oneSmoothTerm$s_age.delta.adj.rsq,
               mygam_changedrsq_oneSmoothTerm$model.adj.r.squared - mygam_intercept$model.adj.r.squared)
  # partial.rsq:
  expected <- wrapper_partialRsq_gam_modelarray(element.subset, modelarray, phenotypes, scalar_name,
                                              FD ~ s(age), FD ~ 1)
  expect_equal(mygam_changedrsq_oneSmoothTerm$s_age.partial.rsq, expected)
  # expect that delta.adj.rsq and partial.rsq are highly correlated
  expect_gt(cor(mygam_changedrsq_oneSmoothTerm$s_age.partial.rsq, mygam_changedrsq_oneSmoothTerm$s_age.delta.adj.rsq),
            0.95)
  
  ##### the statistics should be consistent, when with or without requesting changed.rsq: #####
  mygam_changedrsq_oneSmoothTerm_sex <- ModelArray.gam(FD ~ s(age)+sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                                      full.outputs = TRUE, correct.p.value.smoothTerms = "fdr", correct.p.value.parametricTerms = "fdr",
                                                      changed.rsq.term.index = c(1),
                                                      n_cores = 2, pbar = FALSE)
    # compared to statistics in full outputs:
  colnames_intersect <- intersect(mygam_changedrsq_oneSmoothTerm_sex %>% colnames(),
                                  mygam_fullOutputs %>% colnames())
  expect_equal(mygam_changedrsq_oneSmoothTerm_sex %>% select(colnames_intersect), 
               mygam_fullOutputs %>% select(colnames_intersect))
    # compared to fdr:
  expect_equal(mygam_changedrsq_oneSmoothTerm_sex$sex.p.value.fdr,
               mygam_parametric_pCorrect$sex.p.value.fdr)
  expect_equal(mygam_changedrsq_oneSmoothTerm_sex$s_age.p.value.fdr,
               mygam_smooth_pCorrect$s_age.p.value.fdr)
  
  
  #### more than one term of interest; also, parametric term or smooth term: #####
  mygam_changedRsq_twoSmoothTerm <- ModelArray.gam(FD ~ factorB + s(age) + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                               #var.model = c("dev.expl", "adj.r.squared"),
                                               full.outputs = TRUE,
                                               changed.rsq.term.index = c(1,2,3),
                                               correct.p.value.smoothTerms = "fdr", correct.p.value.parametricTerms = "fdr",
                                               n_cores = 2, pbar = FALSE)
  mygam_changedRsq_twoSmoothTerm_red1 <- ModelArray.gam(FD ~ factorB + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                                     var.model = c( "adj.r.squared"),
                                                     n_cores = 2, pbar = FALSE)
  mygam_changedRsq_twoSmoothTerm_red2 <- ModelArray.gam(FD ~ factorB + s(age), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                                     var.model = c( "adj.r.squared"),
                                                     n_cores = 2, pbar = FALSE)
  mygam_changedRsq_twoSmoothTerm_red3 <- ModelArray.gam(FD ~ s(age) + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                                     var.model = c( "adj.r.squared"),
                                                     n_cores = 2, pbar = FALSE)
  # delta.adj.rsq:
  expect_equal(mygam_changedRsq_twoSmoothTerm$s_age.delta.adj.rsq,
               mygam_changedRsq_twoSmoothTerm$model.adj.r.squared - mygam_changedRsq_twoSmoothTerm_red1$model.adj.r.squared)
  expect_equal(mygam_changedRsq_twoSmoothTerm$s_factorA.delta.adj.rsq,
               mygam_changedRsq_twoSmoothTerm$model.adj.r.squared - mygam_changedRsq_twoSmoothTerm_red2$model.adj.r.squared)
  expect_equal(mygam_changedRsq_twoSmoothTerm$factorB.delta.adj.rsq,
               mygam_changedRsq_twoSmoothTerm$model.adj.r.squared - mygam_changedRsq_twoSmoothTerm_red3$model.adj.r.squared)
  # partial.rsq:
    # test red1:
  expected <- wrapper_partialRsq_gam_modelarray(element.subset, modelarray, phenotypes, scalar_name,
                                                FD ~ factorB + s(age) + s(factorA), FD ~ factorB + s(factorA))
  expect_equal(mygam_changedRsq_twoSmoothTerm$s_age.partial.rsq, expected)
    # test red2:
  expected <- wrapper_partialRsq_gam_modelarray(element.subset, modelarray, phenotypes, scalar_name,
                                                FD ~ factorB + s(age) + s(factorA), FD ~ factorB + s(age))
  expect_equal(mygam_changedRsq_twoSmoothTerm$s_factorA.partial.rsq, expected)
    # test red3:
  expected <- wrapper_partialRsq_gam_modelarray(element.subset, modelarray, phenotypes, scalar_name,
                                                FD ~ factorB + s(age) + s(factorA), FD ~ s(age) + s(factorA))
  expect_equal(mygam_changedRsq_twoSmoothTerm$factorB.partial.rsq, expected)

  
  ##### the statistics should be consistent, when with or without requesting changed.rsq: #####
  mygam_twoSmoothTerm_withoutChangedRsq <- ModelArray.gam(FD ~ factorB + s(age) + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                                          #var.model = c("dev.expl", "adj.r.squared"),
                                                          full.outputs = TRUE,
                                                          correct.p.value.smoothTerms = 'fdr', correct.p.value.parametricTerms = "fdr",
                                                          n_cores = 2, pbar = FALSE)
  # compared to statistics in full outputs:
  colnames_intersect <- intersect(mygam_changedRsq_twoSmoothTerm %>% colnames(),
                                  mygam_twoSmoothTerm_withoutChangedRsq %>% colnames())
  expect_equal(mygam_changedRsq_twoSmoothTerm %>% select(colnames_intersect), 
               mygam_twoSmoothTerm_withoutChangedRsq %>% select(colnames_intersect))
  # compared to fdr:
  expect_equal(mygam_changedRsq_twoSmoothTerm$sex.p.value.fdr,
               mygam_twoSmoothTerm_withoutChangedRsq$sex.p.value.fdr)
  expect_equal(mygam_changedRsq_twoSmoothTerm$s_age.p.value.fdr,
               mygam_twoSmoothTerm_withoutChangedRsq$s_age.p.value.fdr)
  
  #### test that s(age, k=4) with "," in the term label --> see if the column name (s_age.delta.adj.rsq) is as expected #####
  mygam_changedRsq_withComma <- ModelArray.gam(FD ~ s(age, k=4), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                 var.model = c("dev.expl", "adj.r.squared"),
                 changed.rsq.term.index = c(1),
                 n_cores = 2, pbar = FALSE)
  expect_true(c("s_age.delta.adj.rsq", 
                "s_age.partial.rsq") %in% colnames(mygam_changedRsq_withComma) %>% all())
  
  #### invalid request - see my checker in ModelArray.gam() #####
  expect_error(ModelArray.gam(FD ~ s(age, k=4), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                              changed.rsq.term.index = c(0),   # index is < 1
                              n_cores = 2, pbar = FALSE))
  expect_error(ModelArray.gam(FD ~ s(age, k=4), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                              changed.rsq.term.index = c(2),   # index is more than # of terms
                              n_cores = 2, pbar = FALSE))
  expect_error(ModelArray.gam(FD ~ 1, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                              changed.rsq.term.index = c(1),   # invalid formula for changed.rsq
                              n_cores = 2, pbar = FALSE))
  
  #### output of arguments in smooth terms: #####
  # test out te(xxx) instead of s(xxx)
  mygam_changedRsq_te <- ModelArray.gam(FD ~ te(age), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                     changed.rsq.term.index = c(1),  
                                     n_cores = 2, pbar = FALSE)
  expect_true(c("te_age.delta.adj.rsq",
                "te_age.partial.rsq") %in% colnames(mygam_changedRsq_te) %>% all())
  
  # invalid formula: invalid parameters in s such as d = ?
  expect_error(ModelArray.gam(FD ~ s(age, d=1), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                 changed.rsq.term.index = c(1),
                 n_cores = 2, pbar = FALSE))
  
  #### MANUALLY CHECK: whether the printed output is as expected (not able to automatically test via expect_output, probably due to crayon package) #####
  formula <- FD ~ s(age, factorA, fx = FALSE, bs = c("tp", "cr"))   # s(), two terms in s()
  ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                 changed.rsq.term.index = c(1),
                 n_cores = 2, pbar = FALSE)
  # to expect: "s(age,factorA):   k = -1 (default);   fx = FALSE (default);   bs = tp, cr"
  
  formula <- FD ~ s(age, factorA, k=4, bs = c("tp", "tp"))   # s(), two terms in s()
  ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                         changed.rsq.term.index = c(1),
                          n_cores = 2, pbar = FALSE)
  # to expect: "s(age,factorA):   k = 4;   fx = FALSE (default);   bs = tp, tp (default)"

  formula <- FD ~ ti(age, fx = FALSE, bs = c("cr"))   # ti()
  ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                 changed.rsq.term.index = c(1),
                 n_cores = 2, pbar = FALSE)
  # to expect: "ti(age):   fx = FALSE (default);   bs = cr (default)"
  
  formula <- FD ~ ti(age, factorA, fx = TRUE, bs = c("cr", "tp"))   # ti(), two terms in ti()
  ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                 changed.rsq.term.index = c(1),
                 n_cores = 2, pbar = FALSE)
  # to expect: "ti(age,factorA):   fx = TRUE, TRUE;   bs = cr, tp"
  
  # fx should only has one value; otherwise there will be a warning from mgcv::gam()
  
  
  ### invalid formula: invalid interaction term:
  
  
  # formula <- FD ~ s(age*factorA)
  # expect_error(ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
  #                             n_cores = 2, pbar = FALSE))
  # 
  # formula <- FD ~ s(age*factorA) + s(age)  # will generate duplicated smooth terms of s(age) (as it cannot digest * or +)
  # expect_error(ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
  #                             n_cores = 2, pbar = FALSE))
  # 
  # formula <- FD ~ s(age + factorA) + s(age)  # will generate duplicated smooth terms of s(age) (as it cannot digest * or +)
  # 
  # formula <- FD ~ s(age + factorA)    
  
  # formula <- FD ~ s()
  # gam.formula.breakdown <- interpret.gam(formula)
  # ofInterest <- gam.formula.breakdown$smooth.spec[[1]]
  
  
  # s(age * factorA) and s(age + factorA): 
  # mgcv::gam throws a warning; column name is NOT correct for .statistics etc; but correct for delta.adj.rsq; also need to test out writing to .h5
  # may try out interpret.gam(formula) to fix .statistics
  # ModelArray.gam(FD ~ s(age * factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
  #                changed.rsq.term.index = c(1),
  #                n_cores = 2, pbar = FALSE)
  # ModelArray.gam(FD ~ s(age + factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
  #                changed.rsq.term.index = c(1),
  #                n_cores = 2, pbar = FALSE)
  # 
  
  ### check for formula with interaction term #####
  ## s(age, by=oSex):
  formula <- FD ~ oSex + s(age,k=4, fx=TRUE) + s(age, by=oSex, fx=TRUE) + factorB  # ordered factor
  mygam_sby <- ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                        changed.rsq.term.index = c(1,2,3), var.model = c("dev.expl","adj.r.squared"),
                        n_cores = 2, pbar = FALSE)
  # column names as expected:
  expect_true("oSex.L.estimate" %in% colnames(mygam_sby))   # parametric term | L: linear parameter (vs Q: quadratic; C: cubic)
  expect_true("s_age.statistic" %in% colnames(mygam_sby))   # (regular) smooth term
  expect_true("s_age_BYoSexM.p.value" %in% colnames(mygam_sby))  # interaction term | ordered factor, displayed group other than reference group
  expect_true("s_age_BYoSex.delta.adj.rsq" %in% colnames(mygam_sby))  # interaction term's changed.rsq, as it's the term itself, there is no label for group name (such as "M")
  expect_true("s_age_BYoSex.partial.rsq" %in% colnames(mygam_sby)) 
  
  red.formula <- FD ~ oSex + s(age,k=4, fx=TRUE) + factorB
  mygam_sby.red1 <- ModelArray.gam(formula = red.formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                  var.model = c("dev.expl","adj.r.squared"),
                                  n_cores = 2, pbar = FALSE)
  # check if delta.adj.rsq is as expected for interaction term (manually calculate the diff of adj.r.sq):
  expect_equal(mygam_sby$model.adj.r.squared - mygam_sby.red1$model.adj.r.squared,   
               mygam_sby$s_age_BYoSex.delta.adj.rsq)  
  
  red.formula <- FD ~ oSex + s(age, by=oSex, fx=TRUE) + factorB
  mygam_sby.red2 <- ModelArray.gam(formula = red.formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                  var.model = c("dev.expl","adj.r.squared"),
                                  n_cores = 2, pbar = FALSE)
  # check if delta.adj.rsq is as expected for regular smooth term (manually calculate the diff of adj.r.sq):
  # HOWEVER IT MAY DEVIATED FROM ITS TRUE DEFINITION WHEN FORMULA CONTAINS INTERACTION VARIABLES....
  expect_equal(mygam_sby$model.adj.r.squared - mygam_sby.red2$model.adj.r.squared,   
               mygam_sby$s_age.delta.adj.rsq)  
  
  red.formula <- FD ~ s(age,k=4, fx=TRUE) + s(age, by=oSex, fx=TRUE) + factorB
  mygam_sby.red3 <- ModelArray.gam(formula = red.formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                   var.model = c("dev.expl","adj.r.squared"),
                                   n_cores = 2, pbar = FALSE)
  # check if delta.adj.rsq is as expected for parametric term (ordered factor) (manually calculate the diff of adj.r.sq):
  # HOWEVER IT MAY DEVIATED FROM ITS TRUE DEFINITION WHEN FORMULA CONTAINS INTERACTION VARIABLES....
  expect_equal(mygam_sby$model.adj.r.squared - mygam_sby.red3$model.adj.r.squared,   
               mygam_sby$oSex.delta.adj.rsq)  
  
  
  ## ti(x,z):
  formula <- FD ~ ti(age, fx=TRUE) + ti(factorB, fx=TRUE) + ti(age, factorB, fx=TRUE) + factorA
  mygam_tiInteract <- ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                              changed.rsq.term.index = c(3), var.model = c("dev.expl","adj.r.squared"),
                              n_cores = 2, pbar = FALSE)
  expect_true("ti_age.statistic" %in% colnames(mygam_tiInteract))
  expect_true("ti_age_factorB.p.value" %in% colnames(mygam_tiInteract))
  expect_true("ti_age_factorB.delta.adj.rsq" %in% colnames(mygam_tiInteract))
  expect_true("ti_age_factorB.partial.rsq" %in% colnames(mygam_tiInteract))

  red.formula <- FD ~ ti(age, fx=TRUE) + ti(factorB, fx=TRUE)+ factorA
  mygam_tiInteract_red1 <- ModelArray.gam(formula = red.formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                                          var.model = c("adj.r.squared"),
                                          n_cores = 2, pbar = FALSE)
  expect_equal(mygam_tiInteract$model.adj.r.squared -mygam_tiInteract_red1$model.adj.r.squared,
               mygam_tiInteract$ti_age_factorB.delta.adj.rsq)  

  
  ## factorized, but not ordered: - NOT RECOMMEND
  formula <- FD ~ sexFactor + s(age) + s(age, by=sexFactor, fx=TRUE)
  mygam_sby_unordered <- ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
                              changed.rsq.term.index = c(1,2,3), var.model = c("dev.expl","adj.r.squared"),
                              n_cores = 2, pbar = FALSE)
  
  expect_true("s_age.delta.adj.rsq" %in% colnames(mygam_sby_unordered))  # check if correct colname - without other specification in s()
  expect_true("s_age.partial.rsq" %in% colnames(mygam_sby_unordered))
  expect_true("s_age_BYsexFactorF.statistic" %in% colnames(mygam_sby_unordered))  # as unordered, there are terms ending with "F" and "M" afer var name "sex_factor"
  expect_true("sexFactorM.estimate" %in% colnames(mygam_sby_unordered)) 
  expect_true("s_age_BYsexFactor.delta.adj.rsq" %in% colnames(mygam_sby_unordered)) 
  expect_true("s_age_BYsexFactor.partial.rsq" %in% colnames(mygam_sby_unordered)) 
  
  ### test out the functions for generating gam functions: #####
  ## Formula #1:
  myFormula_1 <- generator_gamFormula_factorXsmooth(response.var = "FD", factor.var = "oSex", smooth.var = "age",
                                                    phenotypes = phenotypes)
  myFormula_1$formula # requires visually check
  
  myFormula_2 <- generator_gamFormula_factorXsmooth(response.var = "FD", factor.var = "sex", smooth.var = "age",
                                                    phenotypes = phenotypes, reference.group = "F")
  myFormula_2$formula # requires visually check
  expect_true("osex" %in% colnames(myFormula_2$phenotypes))
  phenotypes_updated <- myFormula_2$phenotypes
  osex.class <- class(phenotypes_updated[["osex"]])
  expect_true(  (length(osex.class) == 2) & (osex.class[1] == "ordered") & (osex.class[2] == "factor")  )
  
  expect_error(generator_gamFormula_factorXsmooth(response.var = "FD", factor.var = "sex", smooth.var = "age",
                                                     phenotypes = phenotypes))   # not ordered factor, and did not provide reference.group
  
  # change k and fx:
  myFormula_3 <- generator_gamFormula_factorXsmooth(response.var = "FD", factor.var = "oSex", smooth.var = "age",
                                                    phenotypes = phenotypes, fx=FALSE, k=4)
  myFormula_3$formula   # requires visually check
   
  
  ## Formula #2:
  myFormula_4 <- generator_gamFormula_continuousInteraction(response.var = "FD", cont1.var = "age", cont2.var = "factorA")
  myFormula_4  # requires visually check
  
  # change k and fx:
  myFormula_5 <- generator_gamFormula_continuousInteraction(response.var = "FD", cont1.var = "age", cont2.var = "factorA",
                                                            fx=FALSE, k=3)
  myFormula_5  # requires visually check
  
  
  ### source file list sanity check #####
  # not the same length:
  phenotypes_wrong1 <- phenotypes[-c(1),]
  expect_error(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes_wrong1, scalar = scalar_name, element.subset = element.subset,
                              n_cores = 2, pbar = FALSE))

# not the same order --> will be handled!
  phenotypes_swap <- phenotypes
  temp <- phenotypes_swap[2,]   # swap row1 and row2
  phenotypes_swap[2,] <- phenotypes_swap[1,]
  phenotypes_swap[1,] <- temp
  mygam_testSwap <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes_swap, scalar = scalar_name, element.subset = element.subset,
                              n_cores = 2, pbar = FALSE)
  expect_equal(mygam_testSwap, mygam_default)

  # change one row --> cannot be matched --> error:
  phenotypes_wrong2 <- phenotypes
  phenotypes_wrong2[1,"source_file"] <- "wrong_file"
  expect_error(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes_wrong2, scalar = scalar_name, element.subset = element.subset,
                              n_cores = 2, pbar = FALSE))



  
  ### debugging:
  #  Error in term[i] <- attr(terms(reformulate(term[i])), "term.labels") : 
  #   replacement has length zero 
  # may because of invalid arguments in smooth term (e.g. d in s(age, d = 1))
  
})

