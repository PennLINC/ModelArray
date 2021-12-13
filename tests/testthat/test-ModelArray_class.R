test_that("ModelArray interface works as expected", {
  # TODO implement tests for accessors
  
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  modelarray <- ModelArray(h5_path, 
                   "scalar_types" = c("FD"), 
                   analysis_names = c("my_analysis"))
  
  expect_s4_class(modelarray, class = "ModelArray")
  
})
