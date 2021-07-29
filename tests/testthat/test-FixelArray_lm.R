test_that("FixelArray's lm works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "FixelArray")
  fa <- FixelArray(h5_path, 
                   "scalar_types" = c("FD"), 
                   analysis_names = c("my_analysis"))
  csv_path <- system.file("extdata", "n50_cohort.csv", package = "FixelArray")
  phenotypes <- read.csv(csv_path)
  scalar_name <- "FD"
  mylm <- FixelArray.lm(FD ~ age + sex, data = fa, phenotypes = phenotypes, scalar = scalar_name, idx = 1:100, n_cores = 2, pbar=FALSE)
  
  expect_equal(mylm$fixel_id, rep(0:99, each=3))
  
  
  # TODO:  check if additional arguments are passed in 
})