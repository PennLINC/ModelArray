test_that("ModelArray interface works as expected", {
  h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")
  num.fixels <- 182581
  num.subj <- 50

  ### get expected values ######
  tryCatch(
    {
      ## try:
      scalar.value.full <- helper_generate_expect_accessors(h5_path)
      # ^ has been realized as matrix, taking ~70MB of memory if using n50_fixels.h5 file
    },
    finally = { # regardless try is successful or not:
      h5closeAll()
    }
  )

  # make sure the expected results have got successfully:
  expect_equal(
    scalar.value.full %>% dim(),
    c(num.fixels, num.subj)
  )

  ### test loading #####
  # loading data:
  modelarray <- ModelArray(h5_path,
    scalar_types = c("FD")
  )

  ### basic checks of ModelArray class and methods #####
  expect_s4_class(modelarray, class = "ModelArray")

  # check show():
  modelarray
  show(modelarray)

  ## check each accessor:
  # scalars():
  expect_equal(scalars(modelarray) %>% length(), 1) # one list of FD
  expect_equal(scalars(modelarray)[["FD"]] %>% dim(), c(num.fixels, num.subj))
  # expect values:
  actual <- scalars(modelarray)[["FD"]] %>% as.matrix() # realize as in-memory matrix, instead of DelayedArray matrix
  dimnames(actual) <- NULL # reset the dimnames (to make it the same as `scalar.value.full`), as the expect is value only, and the dimnames (actually only colnames) will be checked later in sources()
  expect_equal(
    actual, # only checking the values, not the column names (will be checked later in sources())
    scalar.value.full
  )
  # @scalars:
  expect_equal(modelarray@scalars %>% length(), 1) # one list of FD


  # results():
  expect_equal(results(modelarray), list()) # empty

  # @results:
  expect_equal(modelarray@results, list())


  # sources():
  expect_sources_FD <- paste0("FD/sub", as.character(c(1:num.subj)), "_fd.mif")
  expect_equal(
    sources(modelarray)$FD,
    expect_sources_FD
  )

  # @sources:
  expect_equal(modelarray@sources$FD, expect_sources_FD)


  ## other slot:
  # @path:  # there is no accessor for it
  expect_equal(modelarray@path, h5_path)


  ### test out writing results and reloading; also test out results() accessor #####
  tryCatch(
    {
      ## try:
      # copy the .h5 file:
      h5_path_output <- gsub(".h5", "_output.h5", h5_path)
      if (h5_path == h5_path_output) {
        stop("h5_path_output has the same filename as h5_path!")
      }
      file.copy(from = h5_path, to = h5_path_output, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
      modelarray_output <- ModelArray(h5_path_output, scalar_types = c("FD"))

      csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
      phenotypes <- read.csv(csv_path)

      element.subset <- 1:10
      # fit lm:
      mylm <- ModelArray.lm(FD ~ age,
        data = modelarray_output, phenotypes = phenotypes, scalar = "FD", element.subset = element.subset,
        n_cores = 2, pbar = FALSE
      )
      mylm_fulloutput <- ModelArray.lm(FD ~ age,
        data = modelarray_output, phenotypes = phenotypes, scalar = "FD", element.subset = element.subset,
        full.outputs = TRUE,
        n_cores = 2, pbar = FALSE
      )

      ## test: write results:  # multiple results:
      writeResults(h5_path_output, df.output = mylm, analysis_name = "result_lm", overwrite = TRUE)
      writeResults(h5_path_output, df.output = mylm_fulloutput, analysis_name = "result_lm_fulloutput", overwrite = TRUE)

      # re-load:
      modelarray_new <- ModelArray(h5_path_output,
        scalar_types = c("FD"),
        analysis_names = c("result_lm", "result_lm_fulloutput")
      )

      ## test accessor of results():
      # results():
      expect_equal(
        results(modelarray_new)[["result_lm"]][["results_matrix"]] %>% dim(),
        dim(mylm)
      )
      # expected values:
      expect_equal(
        results(modelarray_new)[["result_lm"]][["results_matrix"]] %>% as.data.frame(),
        mylm
      )

      # @results [slot]; also another analysis group:
      expect_equal(
        modelarray_new@results$result_lm_fulloutput[["results_matrix"]] %>% dim(),
        dim(mylm_fulloutput)
      )
      # expected values:
      expect_equal(
        results(modelarray_new)[["result_lm_fulloutput"]][["results_matrix"]] %>% as.data.frame(),
        mylm_fulloutput
      )


      ## test if writeResults's overwrite works:
      writeResults(h5_path_output, df.output = mylm_fulloutput, analysis_name = "result_lm", overwrite = TRUE)
      modelarray_new <- ModelArray(h5_path_output,
        scalar_types = c("FD"),
        analysis_names = c("result_lm", "result_lm_fulloutput")
      )
      expect_equal(
        results(modelarray_new)[["result_lm"]][["results_matrix"]] %>% dim(),
        dim(mylm_fulloutput)
      )
      # expected values:
      expect_equal(
        results(modelarray_new)[["result_lm"]][["results_matrix"]] %>% as.data.frame(),
        mylm_fulloutput
      )
    },
    finally = { # regardless try is successful or not:
      # delete the saved .h5 file with results
      file.remove(h5_path_output)
    }
  )



  # TODO: results with a column of strings - need to be converted into numeric | currently there is no need for .lm and .gam



  ### test out other functions/ utils #####
  # numElementsTotal():
  expect_equal(
    numElementsTotal(modelarray, "FD"),
    num.fixels
  )
})
