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
  grid.subset = 1:10

  ### basic checks #####
  mygam <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                          var.smoothTerms = var.smoothTerms,
                          var.parametricTerms = var.parametricTerms,
                          var.model = var.model,
                          n_cores = 1, pbar = FALSE)

  expect_equal(mygam$grid_id, 0:(max(grid.subset)-1))  # check output$grid_id
  expect_true(is.data.frame(mygam))   # should be data.frame
  expect_equal(as.numeric(dim(mygam)), c(length(grid.subset), 1+1*length(var.smoothTerms) + 2*length(var.parametricTerms) + length(var.model)))  # check the shape

  mygam_default <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                  n_cores = 1, pbar = FALSE)   # default stat outputs
  expect_equal(as.numeric(dim(mygam_default)), c(length(grid.subset),10))

  mygam_fullOutputs <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                      full.outputs = TRUE,    # overwrites the var* arguments below
                                      var.smoothTerms = var.smoothTerms,
                                      var.parametricTerms = var.parametricTerms,
                                      var.model = var.model,
                                      n_cores = 1, pbar = FALSE)
  expect_equal(as.numeric(dim(mygam_fullOutputs)), c(length(grid.subset),24))
  
  ### when there is no term or stat output #####
  # if there is no explicit parametric term (although there will be term "Intercept") or parametric stat:
  mygam_noExplicitParamTerm <- ModelArray.gam(FD ~ s(age), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                  n_cores = 2, pbar = FALSE)  
  mygam_noParamStat <- ModelArray.gam(FD ~ s(age)+sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                  var.parametricTerms = c(),
                                  n_cores = 2, pbar = FALSE)  
  # if there is no smooth term or stat:
  mygam_noSmoothTerm <- ModelArray.gam(FD ~ age, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                  n_cores = 2, pbar = FALSE)   # should without error # expect_output() to test output "there is no smooth term in the requested formula" does not work
  expect_true(c("age.statistic", "age.estimate") %in% colnames(mygam_noSmoothTerm) %>% all())
  
  mygam_noSmoothStat <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                       var.smoothTerms = c(),
                                       n_cores = 2, pbar = FALSE) 
  expect_false("age.statistic" %in% colnames(mygam_noSmoothStat))
  
  # if there is no model stat:
  mygam_noModelStat <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                      var.model = c(),
                                      n_cores = 2, pbar = FALSE) 
  
  # no any model stat:
  expect_error(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                      var.smoothTerms = c(), var.parametricTerms = c(), var.model = c(),
                                      n_cores = 2, pbar = FALSE))
  
  
  ## multiple smooth terms:
  # s(age) + s(factorA)   # cannot use s(sex) as sex are characters (M or F), not values!
  mygam_sAge_sFactorA <- ModelArray.gam(FD ~ s(age) + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                n_cores = 2, pbar = FALSE)
  
  # # s(age + factorA)   # not good example
  # expect_warning(ModelArray.gam(FD ~ s(age + factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
  #                n_cores = 2, pbar = FALSE))  # note: the column name will be s-age instead of s-age-sex
  
  
  
  ### Test whether the validity of list of var is checked: #####
  expect_error(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                              var.smoothTerms = c("wrong_name"),
                              n_cores = 2, pbar = FALSE))
  # duplicated
  temp <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                         var.parametricTerms = c(var.parametricTerms, "p.value"),
                         n_cores = 2, pbar = FALSE)
  expect_equal(temp, mygam_default)
  
  
  
  ### different arguments in GAM #####
  ## different settings in formula:
  # s(k=?):
  mygam_k4 <- ModelArray.gam(FD ~ s(age, k=4) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                             n_cores = 2, pbar = FALSE)
  expect_false(all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_k4 %>% dplyr::select("s_age.statistic"))
               %>% isTRUE())    # should be different when k is different
  
  # s(fx=TRUE) vs default (fx=FALSE): 1) different stat; 
  mygam_fxTRUE <- ModelArray.gam(FD ~ s(age, fx=TRUE) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                 n_cores = 2, pbar = FALSE)
  expect_false(all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_fxTRUE %>% dplyr::select("s_age.statistic"))
               %>% isTRUE())
  
  # TODO: if force saving edf when fx=FALSE: 2) check if it's saved or not

  # s(bs=?)
  mygam_bsCR <- ModelArray.gam(FD ~ s(age, bs="cr") + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                               n_cores = 2, pbar = FALSE)
  expect_false(all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_bsCR %>% dplyr::select("s_age.statistic"))
               %>% isTRUE()) 
  
  ## different settings in mgcv::gam()'s additional arguments: test if the arguments have been passed into analyseOneGrid.gam()
  # method: 
  mygam_methodREML <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                     method = "REML",  # default = "GCV.Cp"
                                     n_cores = 2, pbar = FALSE)
  expect_false(all_equal(mygam_default %>% dplyr::select("s_age.statistic"),
                         mygam_methodREML %>% dplyr::select("s_age.statistic"))
               %>% isTRUE())
  ### Test n_cores, pbar work: ######
  # n_cores = 2: 
  mygam_pbarFalse_ncores2 <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                           n_cores = 2, pbar = FALSE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarFalse_ncores2)
  
  # pbar:
  mygam_pbarTrue_ncores2 <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                  n_cores = 2, pbar = TRUE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarTrue_ncores2)
  
  mygam_pbarTrue_ncores1 <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                           n_cores = 1, pbar = TRUE)   # default stat outputs
  expect_equal(mygam_default, mygam_pbarTrue_ncores1)
  
  
  ### Test: p.value correction: #####
  mygam_parametric_pCorrect <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                 correct.p.value.parametricTerms = c("fdr","bonferroni"),
                 n_cores = 2, pbar = FALSE) 
  expect_equal(mygam_parametric_pCorrect$sexM.p.value.fdr,
               mygam_parametric_pCorrect$sexM.p.value %>% stats::p.adjust("fdr"))
  expect_equal(mygam_parametric_pCorrect$sexM.p.value.bonferroni,
               mygam_parametric_pCorrect$sexM.p.value %>% stats::p.adjust("bonferroni"))
  
  mygam_smooth_pCorrect <- ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                              correct.p.value.smoothTerms = c("fdr","bonferroni"),
                              n_cores = 2, pbar = FALSE) 
  expect_equal(mygam_smooth_pCorrect$s_age.p.value.fdr,
               mygam_smooth_pCorrect$s_age.p.value %>% stats::p.adjust("fdr"))
  expect_equal(mygam_smooth_pCorrect$s_age.p.value.bonferroni,
               mygam_smooth_pCorrect$s_age.p.value %>% stats::p.adjust("bonferroni"))

  expect_error(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                              correct.p.value.parametricTerms = c("wrong_correct"),
                              n_cores = 2, pbar = FALSE))
  expect_warning(ModelArray.gam(FD ~ s(age) + sex, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                var.smoothTerms = c("statistic"),  # no p.value
                                correct.p.value.smoothTerms = c("fdr"),
                                n_cores = 2, pbar = FALSE))
  
  ### Test: eff.size #####
  # one term of interest: reduced model will be FD ~ 1
  # also, to test whether the eff.size is calculated correctly
  # also, to test whether without "," in s(), the column name could be correctly "s_age.eff.size"
  mygam_effsize_oneSmoothTerm <- ModelArray.gam(FD ~ s(age), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                               var.model = c("dev.expl", "adj.r.squared"),
                                               eff.size.term.index = c(1),
                                               n_cores = 2, pbar = FALSE)
  mygam_intercept <- ModelArray.gam(FD ~ 1, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                    var.model = c("adj.r.squared"),
                                    n_cores = 2, pbar = FALSE)
  expect_equal(mygam_effsize_oneSmoothTerm$s_age.eff.size,
               mygam_effsize_oneSmoothTerm$model.adj.r.squared - mygam_intercept$model.adj.r.squared)
  
  # more than one term of interest; also, parametric term or smooth term:
  mygam_effsize_twoSmoothTerm <- ModelArray.gam(FD ~ factorB + s(age) + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                               var.model = c("dev.expl", "adj.r.squared"),
                                               eff.size.term.index = c(1,2,3),
                                               n_cores = 2, pbar = FALSE)
  mygam_effsize_twoSmoothTerm_red1 <- ModelArray.gam(FD ~ factorB + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                                     var.model = c( "adj.r.squared"),
                                                     n_cores = 2, pbar = FALSE)
  mygam_effsize_twoSmoothTerm_red2 <- ModelArray.gam(FD ~ factorB + s(age), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                                     var.model = c( "adj.r.squared"),
                                                     n_cores = 2, pbar = FALSE)
  mygam_effsize_twoSmoothTerm_red3 <- ModelArray.gam(FD ~ s(age) + s(factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                                     var.model = c( "adj.r.squared"),
                                                     n_cores = 2, pbar = FALSE)
  expect_equal(mygam_effsize_twoSmoothTerm$s_age.eff.size,
               mygam_effsize_twoSmoothTerm$model.adj.r.squared - mygam_effsize_twoSmoothTerm_red1$model.adj.r.squared)
  expect_equal(mygam_effsize_twoSmoothTerm$s_factorA.eff.size,
               mygam_effsize_twoSmoothTerm$model.adj.r.squared - mygam_effsize_twoSmoothTerm_red2$model.adj.r.squared)
  expect_equal(mygam_effsize_twoSmoothTerm$factorB.eff.size,
               mygam_effsize_twoSmoothTerm$model.adj.r.squared - mygam_effsize_twoSmoothTerm_red3$model.adj.r.squared)

  # test that s(age, k=4) with "," in the term label --> see if the column name (s_age.eff.size) is as expected
  mygam_effsize_withComma <- ModelArray.gam(FD ~ s(age, k=4), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                 var.model = c("dev.expl", "adj.r.squared"),
                 eff.size.term.index = c(1),
                 n_cores = 2, pbar = FALSE)
  expect_true("s_age.eff.size" %in% colnames(mygam_effsize_withComma))
  
  # invalid request - see my checker in ModelArray.gam()
  expect_error(ModelArray.gam(FD ~ s(age, k=4), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                              eff.size.term.index = c(0),   # index is < 1
                              n_cores = 2, pbar = FALSE))
  expect_error(ModelArray.gam(FD ~ s(age, k=4), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                              eff.size.term.index = c(2),   # index is more than # of terms
                              n_cores = 2, pbar = FALSE))
  expect_error(ModelArray.gam(FD ~ 1, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                              eff.size.term.index = c(1),   # invalid formula for effect size
                              n_cores = 2, pbar = FALSE))
  
  ### output of arguments in smooth terms:
  # test out te(xxx) instead of s(xxx)
  mygam_effsize_te <- ModelArray.gam(FD ~ te(age), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                     eff.size.term.index = c(1),  
                                     n_cores = 2, pbar = FALSE)
  expect_true("te_age.eff.size" %in% colnames(mygam_effsize_te))
  
  # invalid formula: invalid parameters in s such as d = ?
  expect_error(ModelArray.gam(FD ~ s(age, d=1), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                 eff.size.term.index = c(1),
                 n_cores = 2, pbar = FALSE))
  
  ## whether the printed output is as expected (not able to automatically test via expect_output, probably due to crayon package)
  formula <- FD ~ s(age, factorA, fx = FALSE, bs = c("tp", "cr"))   # s(), two terms in s()
  ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                 eff.size.term.index = c(1),
                 n_cores = 2, pbar = FALSE)
  # to expect: "s(age,factorA):   k = -1 (default);   fx = FALSE (default);   bs = tp, cr"
  
  formula <- FD ~ s(age, factorA, k=4, bs = c("tp", "tp"))   # s(), two terms in s()
  ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                         eff.size.term.index = c(1),
                          n_cores = 2, pbar = FALSE)
  # to expect: "s(age,factorA):   k = 4;   fx = FALSE (default);   bs = tp, tp (default)"

  formula <- FD ~ ti(age, fx = FALSE, bs = c("cr"))   # ti()
  ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                 eff.size.term.index = c(1),
                 n_cores = 2, pbar = FALSE)
  # to expect: "ti(age):   fx = FALSE (default);   bs = cr (default)"
  
  formula <- FD ~ ti(age, factorA, fx = TRUE, bs = c("cr", "tp"))   # ti(), two terms in ti()
  ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                 eff.size.term.index = c(1),
                 n_cores = 2, pbar = FALSE)
  # to expect: "ti(age,factorA):   fx = TRUE, TRUE;   bs = cr, tp"
  
  # fx should only has one value; otherwise there will be a warning from mgcv::gam()
  
  
  ### invalid formula: invalid interaction term:
  
  
  # formula <- FD ~ s(age*factorA)
  # expect_error(ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
  #                             n_cores = 2, pbar = FALSE))
  # 
  # formula <- FD ~ s(age*factorA) + s(age)  # will generate duplicated smooth terms of s(age) (as it cannot digest * or +)
  # expect_error(ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
  #                             n_cores = 2, pbar = FALSE))
  # 
  # formula <- FD ~ s(age + factorA) + s(age)  # will generate duplicated smooth terms of s(age) (as it cannot digest * or +)
  # 
  # formula <- FD ~ s(age + factorA)    
  
  # formula <- FD ~ s()
  # gam.formula.breakdown <- interpret.gam(formula)
  # ofInterest <- gam.formula.breakdown$smooth.spec[[1]]
  
  
  # s(age * factorA) and s(age + factorA): 
  # mgcv::gam throws a warning; column name is NOT correct for .statistics etc; but correct for eff.size; also need to test out writing to .h5
  # may try out interpret.gam(formula) to fix .statistics
  # ModelArray.gam(FD ~ s(age * factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
  #                eff.size.term.index = c(1),
  #                n_cores = 2, pbar = FALSE)
  # ModelArray.gam(FD ~ s(age + factorA), data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
  #                eff.size.term.index = c(1),
  #                n_cores = 2, pbar = FALSE)
  # 
  
  ### check for formula with interaction term #####
  ## s(age, by=oSex):
  formula <- FD ~ oSex + s(age,k=4, fx=TRUE) + s(age, by=oSex, fx=TRUE) + factorB  # ordered factor
  mygam_sby <- ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                        eff.size.term.index = c(1,2,3), var.model = c("dev.expl","adj.r.squared"),
                        n_cores = 2, pbar = FALSE)
  # column names as expected:
  expect_true("oSex.L.estimate" %in% colnames(mygam_sby))   # parametric term | L: linear parameter (vs Q: quadratic; C: cubic)
  expect_true("s_age.statistic" %in% colnames(mygam_sby))   # (regular) smooth term
  expect_true("s_age_BYoSexM.p.value" %in% colnames(mygam_sby))  # interaction term | ordered factor, displayed group other than reference group
  expect_true("s_age_BYoSex.eff.size" %in% colnames(mygam_sby))  # interaction term's effect size, as it's the term itself, there is no label for group name (such as "M")
  
  red.formula <- FD ~ oSex + s(age,k=4, fx=TRUE) + factorB
  mygam_sby.red1 <- ModelArray.gam(formula = red.formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                  var.model = c("dev.expl","adj.r.squared"),
                                  n_cores = 2, pbar = FALSE)
  # check if effect size is as expected for interaction term (manually calculate the diff of adj.r.sq):
  expect_equal(mygam_sby$model.adj.r.squared - mygam_sby.red1$model.adj.r.squared,   
               mygam_sby$s_age_BYoSex.eff.size)  
  
  red.formula <- FD ~ oSex + s(age, by=oSex, fx=TRUE) + factorB
  mygam_sby.red2 <- ModelArray.gam(formula = red.formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                  var.model = c("dev.expl","adj.r.squared"),
                                  n_cores = 2, pbar = FALSE)
  # check if effect size is as expected for regular smooth term (manually calculate the diff of adj.r.sq):
  # HOWEVER IT MAY DEVIATED FROM ITS TRUE DEFINITION WHEN FORMULA CONTAINS INTERACTION VARIABLES....
  expect_equal(mygam_sby$model.adj.r.squared - mygam_sby.red2$model.adj.r.squared,   
               mygam_sby$s_age.eff.size)  
  
  red.formula <- FD ~ s(age,k=4, fx=TRUE) + s(age, by=oSex, fx=TRUE) + factorB
  mygam_sby.red3 <- ModelArray.gam(formula = red.formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                   var.model = c("dev.expl","adj.r.squared"),
                                   n_cores = 2, pbar = FALSE)
  # check if effect size is as expected for parametric term (ordered factor) (manually calculate the diff of adj.r.sq):
  # HOWEVER IT MAY DEVIATED FROM ITS TRUE DEFINITION WHEN FORMULA CONTAINS INTERACTION VARIABLES....
  expect_equal(mygam_sby$model.adj.r.squared - mygam_sby.red3$model.adj.r.squared,   
               mygam_sby$oSex.eff.size)  
  
  
  ## ti(x,z):
  formula <- FD ~ ti(age, fx=TRUE) + ti(factorB, fx=TRUE) + ti(age, factorB, fx=TRUE) + factorA
  mygam_tiInteract <- ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                              eff.size.term.index = c(3), var.model = c("dev.expl","adj.r.squared"),
                              n_cores = 2, pbar = FALSE)
  expect_true("ti_age.statistic" %in% colnames(mygam_tiInteract))
  expect_true("ti_age_factorB.p.value" %in% colnames(mygam_tiInteract))
  expect_true("ti_age_factorB.eff.size" %in% colnames(mygam_tiInteract))

  red.formula <- FD ~ ti(age, fx=TRUE) + ti(factorB, fx=TRUE)+ factorA
  mygam_tiInteract_red1 <- ModelArray.gam(formula = red.formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                                          var.model = c("adj.r.squared"),
                                          n_cores = 2, pbar = FALSE)
  expect_equal(mygam_tiInteract$model.adj.r.squared -mygam_tiInteract_red1$model.adj.r.squared,
               mygam_tiInteract$ti_age_factorB.eff.size)  

  
  ## factorized, but not ordered: - NOT RECOMMEND
  formula <- FD ~ sexFactor + s(age) + s(age, by=sexFactor, fx=TRUE)
  mygam_sby_unordered <- ModelArray.gam(formula = formula, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, grid.subset = grid.subset,
                              eff.size.term.index = c(1,2,3), var.model = c("dev.expl","adj.r.squared"),
                              n_cores = 2, pbar = FALSE)
  
  expect_true("s_age.eff.size" %in% colnames(mygam_sby_unordered))  # check if correct colname - without other specification in s()
  expect_true("s_age_BYsexFactorF.statistic" %in% colnames(mygam_sby_unordered))  # as unordered, there are terms ending with "F" and "M" afer var name "sex_factor"
  expect_true("sexFactorM.estimate" %in% colnames(mygam_sby_unordered)) 
  expect_true("s_age_BYsexFactor.eff.size" %in% colnames(mygam_sby_unordered)) 

  
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
  
  ### debugging:
  #  Error in term[i] <- attr(terms(reformulate(term[i])), "term.labels") : 
  #   replacement has length zero 
  # may because of invalid arguments in smooth term (e.g. d in s(age, d = 1))
  
})

