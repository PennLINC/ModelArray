test_that("FixelArray's lm works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "FixelArray")
  fa <- FixelArray(h5_path, 
                   "scalar_types" = c("FD"), 
                   analysis_names = c("my_analysis"))
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "FixelArray")
  phenotypes <- read.csv(csv_path)
  scalar_name <- "FD"
  mylm <- FixelArray.lm(FD ~ age, data = fa, phenotypes = phenotypes, scalar = scalar_name, idx = 1:100, n_cores = 2, pbar=FALSE)
  
  expect_equal(mylm$fixel_id, rep(0:99, each=2))   # check output$fixel_id 
  expect_true(is.data.frame(mylm))  # should be data.frame
  expect_true(tibble::is_tibble(mylm))   # should be tibble
  expect_equal(as.numeric(dim(mylm)), c(100*2,6)) # check shape
  
  mylm_age_sex <- FixelArray.lm(FD ~ age + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, idx = 1:100, n_cores = 2, pbar=FALSE)
  expect_equal(mylm_age_sex$fixel_id, rep(0:99, each=3))   # check output$fixel_id 
  expect_equal(as.numeric(dim(mylm_age_sex)), c(100*3,6)) 
  
  expect_false(all.equal(mylm %>% dplyr::filter(term=="age") %>% dplyr::select("estimate"), 
                         mylm_age_sex %>% dplyr::filter(term=="age") %>% dplyr::select("estimate"))
               %>% isTRUE())   # expect not identical between two models
  
  
  # TODO: how about other variables? factorA is literally correlated with age; factorB is another random variable
  # factor A is fully correlated with age, expecting testing results are NA:
  mylm_age_factorA <- FixelArray.lm(FD ~ age + factorA, data = fa, phenotypes = phenotypes, scalar = scalar_name, idx = 1:100, n_cores = 2, pbar=FALSE)
  temp <- mylm_age_factorA %>% dplyr::filter(term=="factorA")%>% dplyr::select(-term)  # only extracting column "factorA"
  expect_equal(temp$fixel_id, c(0:99))    # test that $fixel_id is 0:99
  temp <- mylm_age_factorA %>% dplyr::filter(term=="factorA")%>% dplyr::select(-term) %>% dplyr::select(-fixel_id)   # also removing $fixel_id, now should be all NA
  expect_true(all(is.na(temp)))
  
  # factor B is not correlated with age, so not expecting NA:
  mylm_age_factorB <- FixelArray.lm(FD ~ age + factorB, data = fa, phenotypes = phenotypes, scalar = scalar_name, idx = 1:100, n_cores = 2, pbar=FALSE)
  temp <- mylm_age_factorB %>% dplyr::filter(term=="factorB")%>% dplyr::select(-term) %>% dplyr::select(-fixel_id)  
  expect_false(all(is.na(temp)))
  
  
  # TODO: how about other arguments of lm?
})


# TODO:  check if additional arguments are passed in