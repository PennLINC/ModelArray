test_that("FixelArray interface works as expected", {
  # TODO implement tests for accessors
  
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "FixelArray")
  fa <- FixelArray(h5_path, 
                   "scalar_types" = c("FD"), 
                   analysis_names = c("my_analysis"))
  
  expect_s4_class(fa, class = "FixelArray")
  
})
