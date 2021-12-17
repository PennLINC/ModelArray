test_that("ModelArray interface works as expected", {
  
  h5_path_n25x2 <- system.file("extdata", "n25x2_fixels.h5", package = "ModelArray")
  num.fixels <- 182581
  
  ### test loading #####
  # only loading one group of scalar data: FD, with 25subjects
  modelarray_n25 <- ModelArray(h5_path_n25x2, 
                               scalar_types = c("FD"))
  
  # loading two groups of scalar data: FD and FD_fake, each with 25 subjects
  modelarray_n25x2 <- ModelArray(h5_path_n25x2, 
                               scalar_types = c("FD", "FD_fake"))
  
  ### basic checks of ModelArray class and methods #####
  expect_s4_class(modelarray_n25, class = "ModelArray")
  expect_s4_class(modelarray_n25x2, class = "ModelArray")
  
  # check show():
  modelarray_n25
  modelarray_n25x2
  show(modelarray_n25)
  show(modelarray_n25x2)
  
  ## check each accessor:
  # scalars():
  expect_equal(scalars(modelarray_n25) %>% length(), 1)  # one list of FD
  expect_equal(scalars(modelarray_n25x2) %>% length(), 2)  # FD and FD_fake
  expect_equal(scalars(modelarray_n25)[["FD"]] %>% dim(), c(num.fixels, 25)) 
  expect_equal(scalars(modelarray_n25x2)[["FD_fake"]] %>% dim(), c(num.fixels, 25)) 
  # @scalars:
  expect_equal(modelarray_n25@scalars %>% length(), 1)  # one list of FD
  expect_equal(modelarray_n25x2@scalars %>% length(), 2)  # FD and FD_fake
  
  # results():
  expect_equal(results(modelarray_n25), list())    # empty
  expect_equal(results(modelarray_n25x2), list()) 
  # @results:
  expect_equal(modelarray_n25@results, list()) 
  expect_equal(modelarray_n25x2@results, list()) 
  
  # subjects(): 
  expect_subjlist_FD <- paste0("sub",as.character(c(1:25)))
  expect_subjlist_FDfake <- paste0("sub",as.character(c(26:50)))
    
  expect_equal(subjects(modelarray_n25)$FD, expect_subjlist_FD)
  expect_equal(subjects(modelarray_n25x2)$FD, expect_subjlist_FD)
  expect_equal(subjects(modelarray_n25x2)$FD_fake, expect_subjlist_FDfake)
  # @subjects:
  expect_equal(modelarray_n25@subjects$FD, expect_subjlist_FD)
  
  ## other slot:
  # @path:  # there is no accessor for it
  expect_equal(modelarray_n25@path, h5_path_n25x2)

  
  ### test out writing results and reloading; also test out results() accessor #####
  tryCatch(
    {
      ## try:
      # copy the .h5 file:
      h5_path_n25x2_output <- gsub(".h5","_output.h5",h5_path_n25x2)
      if (h5_path_n25x2 == h5_path_n25x2_output) {
        stop("h5_path_n25x2_output has the same filename as h5_path_n25x2!")
      }
      file.copy(from=h5_path_n25x2, to=h5_path_n25x2_output, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
      modelarray_output <- ModelArray(h5_path_n25x2_output, scalar_types = c("FD","FD_fake"))
      
      csv_path <- system.file("extdata", "n25x2_cohort.csv", package = "ModelArray")
      phenotypes <- read.csv(csv_path)
      phenotypes_FD <- phenotypes[phenotypes$scalar_name == "FD",]
      phenotypes_FDfake <- phenotypes[phenotypes$scalar_name == "FD_fake",]
        
      colname_subject_id <- "subject_id"
      element.subset <- 1:10
      # fit lm:
      mylm <- ModelArray.lm(FD ~ age, data = modelarray_output, phenotypes = phenotypes_FD, scalar = "FD", element.subset = element.subset, 
                            colname.subjid = colname_subject_id,
                            n_cores = 2, pbar=FALSE)
      mylm_fulloutput <- ModelArray.lm(FD ~ age, data = modelarray_output, phenotypes = phenotypes_FD, scalar = "FD", element.subset = element.subset, 
                                       colname.subjid = colname_subject_id,
                                       full.outputs = TRUE,
                                       n_cores = 2, pbar=FALSE)
      mylm_fulloutput_fakeFD <- ModelArray.lm(FD_fake ~ age, data = modelarray_output, phenotypes = phenotypes_FDfake, scalar = "FD_fake", element.subset = element.subset, 
                                              colname.subjid = colname_subject_id,
                                              n_cores = 2, pbar=FALSE)
      ## test: write results:  # multiple results:
      writeResults(h5_path_n25x2_output, df.output = mylm, analysis_name="result_lm", overwrite=TRUE)
      writeResults(h5_path_n25x2_output, df.output = mylm, analysis_name="result_lm_FDfake", overwrite=TRUE)
      
      # re-load:
      modelarray_new <- ModelArray(h5_path_n25x2_output, scalar_types = c("FD","FD_fake"),
                                   analysis_names = c("result_lm", "result_lm_FDfake"))
      
      ## test accessor of results():
      # results():
      expect_equal(results(modelarray_new)[["result_lm"]][["results_matrix"]] %>% dim(),
                   dim(mylm))
      # @results [slot]; also another analysis group:
      expect_equal(modelarray_new@results$result_lm_FDfake[["results_matrix"]] %>% dim(),
                   dim(mylm_fulloutput_fakeFD))
      
      
      ## test if writeResults's overwrite works:    
      writeResults(h5_path_n25x2_output, df.output = mylm_fulloutput, analysis_name="result_lm", overwrite=TRUE)
      modelarray_new <- ModelArray(h5_path_n25x2_output, scalar_types = c("FD","FD_fake"),
                                   analysis_names = c("result_lm", "result_lm_FDfake"))
      expect_equal(results(modelarray_new)[["result_lm"]][["results_matrix"]] %>% dim(),
                   dim(mylm_fulloutput))
      
    },
    finally = {  # regardless try is successful or not:
      # delete the saved .h5 file with results
      file.remove(h5_path_n25x2_output)
      
    }
  )
  
  
  
  # TODO: results with a column of strings - need to be converted into numeric | currently there is no need for .lm and .gam
  

})
