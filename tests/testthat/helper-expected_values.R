# This is to calculate expected statistical results but without using ModelArray (e.g. ModelArray(), ModelArray.lm() or ModelArray.gam()).

#' The function for calculating the expected values for ModelArray.lm()
#' 
#' Details:
#' this won't include p-value corrections - as corrections depend on other calculated fixels; to check p-value corrections, directly calculate in the testthat tests.
calcu_stat_lm <- function(formula, data, i_element, ...){
  
  arguments_lm <- list(...)
  arguments_lm$formula <- formula
  arguments_lm$data <- data
  
  onemodel <- do.call(stats::lm, arguments_lm)   # explicitly passing arguments into lm, to avoid error of argument "weights"
  
  onemodel.tidy <- onemodel %>% broom::tidy()
  onemodel.glance <- onemodel %>% broom::glance()
  
  # adjust:
  onemodel.tidy$term[onemodel.tidy$term == "(Intercept)"] <- "Intercept"  # change the term name from "(Intercept)" to "Intercept"
  onemodel.glance <- onemodel.glance %>% mutate(term="model")   # add a column 
  
  # list of column names to keep:
  var.terms <- colnames(onemodel.tidy)
  var.terms <- var.terms[var.terms != "term"];   # remove "term" which is not a valid stat output
  
  var.model <- colnames(onemodel.glance)
  var.model <- var.model[var.model != "term"];   # remove "term" which is not a valid stat output
  
  # turn into one row:
  onemodel.tidy.onerow <- onemodel.tidy %>% tidyr::pivot_wider(names_from = term,
                                                               values_from = all_of(var.terms),
                                                               names_glue = "{term}.{.value}")
  onemodel.glance.onerow <- onemodel.glance %>%  tidyr::pivot_wider(names_from = term, 
                                                                    values_from = all_of(var.model),
                                                                    names_glue = "{term}.{.value}")
  
  onemodel.onerow <- dplyr::bind_cols(onemodel.tidy.onerow, onemodel.glance.onerow)
  
  # add a column of element ids:
  colnames.temp <- colnames(onemodel.onerow)
  onemodel.onerow <- onemodel.onerow %>% tibble::add_column(element_id = i_element-1, .before = colnames.temp[1])   # add as the first column
  
  # change from tibble to a data.frame:
  onemodel.onerow <- onemodel.onerow %>% as.data.frame()  
  
  onemodel.onerow
}


#' The function for calculating the expected values for ModelArray.gam()
#' 
calcu_stat_gam <- function(formula, data, i_element = idx.fixel.gam, ...) {
  onemodel <- mgcv::gam(formula = formula, data = data, ...)
  
  # in this expected value function, I use "summary.gam()"
  
  onemodel.tidy.smoothTerms <- onemodel %>% broom::tidy(parametric = FALSE) # needs to use broom::tidy instead of broom::tidy.gam() - the latter one is not an exported object from 'namespace:broom'
  onemodel.tidy.parametricTerms <- onemodel %>% broom::tidy(parametric = TRUE)
  onemodel.glance <- onemodel %>% broom::glance()   # needs to use broom::glance instead of broom::glance.gam()!
  onemodel.summary <- onemodel %>% summary.gam()
  # add additional model's stat to onemodel.glance():
  onemodel.glance[["adj.r.squared"]] <- onemodel.summary$r.sq
  onemodel.glance[["dev.expl"]] <- onemodel.summary$dev.expl
  
  sp.criterion.attr.name <- onemodel.summary$sp.criterion %>% attr(which = "name")   # get the attr name  # maually checked: default (GCV.Cp) and "REML", and the *.attr.name string is correct
  onemodel.glance[["sp.criterion"]] <- onemodel.summary$sp.criterion[[ sp.criterion.attr.name ]]   # use the attr name to extract the value
  onemodel.glance[["scale"]] <- onemodel.summary$scale   # scale estimate
  
  #num.smoothTerms <- onemodel.summary$m   # The number of smooth terms in the model.
  
  # adjust:
  if (nrow(onemodel.tidy.smoothTerms) > 0) {   # if there is any smooth term   # NOTE: I used different method for detecting if there is smooth term
    onemodel.tidy.smoothTerms$term[onemodel.tidy.smoothTerms$term == "(Intercept)"] <- "Intercept"  # change the term name from "(Intercept)" to "Intercept"  
  }
  if (nrow(onemodel.tidy.parametricTerms) > 0) {  # if there is any parametric term
    onemodel.tidy.parametricTerms$term[onemodel.tidy.parametricTerms$term == "(Intercept)"] <- "Intercept"  # change the term name from "(Intercept)" to "Intercept"
  }
  
  
  # change from s(x) to s_x: (could be s, te, etc); from s(x):oFactor to s_x_BYoFactor; from ti(x,z) to ti_x_z
  if (nrow(onemodel.tidy.smoothTerms) > 0) {   # if there is any smooth term   # NOTE: I used different method for detecting if there is smooth term
    for (i_row in 1:nrow(onemodel.tidy.smoothTerms)) {  
      # step 1: change from s(x) to s_x
      term_name <- onemodel.tidy.smoothTerms$term[i_row]
      str_list <- strsplit(term_name, split="[()]")[[1]]
      
      str <- str_list[2]   # extract string between ()
      smooth_name <- str_list[1]   # "s" or some other smooth method type such as "te"
      str_valid <- paste0(smooth_name, "_",str)
      
      if (length(str_list)>2) {   # there is string after variable name
        str_valid <- paste0(str_valid, "_",
                            paste(str_list[3:length(str_list)], collapse=""))   # combine rest of strings
      }   
      
      # detect ":", and change to "BY"   # there is "_" replacing for ")" in "s()" already
      str_valid <- gsub(":", "BY", str_valid, fixed=TRUE)
      
      # detect ",", and change to "_"
      str_valid <- gsub(",", "_", str_valid, fixed=TRUE)
      
      onemodel.tidy.smoothTerms$term[i_row] <- str_valid
    }
  }
  
  onemodel.glance <- onemodel.glance %>% mutate(term="model")   # add a column 
  
  
  # check if the onemodel.* does not have real statistics (but only a column of 'term')
  temp_colnames <- onemodel.tidy.smoothTerms %>% colnames()
  temp <- union(temp_colnames, "term")    # union of colnames and "term"; if colnames only has "term" or lengt of 0 (tibble()), union = "term", all(union)=TRUE; otherwise, if there is colnames other than "term", all(union) = c(TRUE, FALSE, ...)
  if (all(temp == "term")) onemodel.tidy.smoothTerms <- tibble()   # just an empty tibble (so below, all(dim(onemodel.tidy.smoothTerms)) = FALSE)
  
  temp_colnames <- onemodel.tidy.parametricTerms %>% colnames()
  temp <- union(temp_colnames, "term")     
  if (all(temp == "term")) onemodel.tidy.parametricTerms <- tibble()   # just an empty tibble
  
  temp_colnames <- onemodel.glance %>% colnames()
  temp <- union(temp_colnames, "term")    
  if (all(temp == "term")) onemodel.glance <- tibble()   # just an empty tibble
  
  
  ## flatten:
  # list of column names to keep: for each var.*, remove "term" which is not a valid stat output
  var.smoothTerms <- colnames(onemodel.tidy.smoothTerms)
  var.smoothTerms <- var.smoothTerms[var.smoothTerms != "term"]  
  
  var.parametricTerms <- colnames(onemodel.tidy.parametricTerms)
  var.parametricTerms <- var.parametricTerms[var.parametricTerms != "term"]
  
  var.model <- colnames(onemodel.glance)
  var.model <- var.model[var.model != "term"]
  
  # flatten:
  if (nrow(onemodel.tidy.smoothTerms) > 0) {
    onemodel.tidy.smoothTerms.onerow <- onemodel.tidy.smoothTerms %>% tidyr::pivot_wider(names_from = term,
                                                                                         values_from = all_of(var.smoothTerms),
                                                                                         names_glue = "{term}.{.value}")
  } else {
    onemodel.tidy.smoothTerms.onerow <- onemodel.tidy.smoothTerms
  }
  
  if (nrow(onemodel.tidy.parametricTerms)>0) {
    onemodel.tidy.parametricTerms.onerow <- onemodel.tidy.parametricTerms %>% tidyr::pivot_wider(names_from = term,
                                                                                                 values_from = all_of(var.parametricTerms),
                                                                                                 names_glue = "{term}.{.value}")
  } else {
    onemodel.tidy.parametricTerms.onerow <- onemodel.tidy.parametricTerms
  }
  
  if (nrow(onemodel.glance) >0) {
    onemodel.glance.onerow <- onemodel.glance %>%  tidyr::pivot_wider(names_from = term, 
                                                                      values_from = all_of(var.model),
                                                                      names_glue = "{term}.{.value}")
  } else {
    onemodel.glance.onerow <- onemodel.glance
  }
  
  ## bind together:
  onemodel.onerow <- bind_cols_check_emptyTibble(onemodel.tidy.smoothTerms.onerow, 
                                                 onemodel.tidy.parametricTerms.onerow)
  onemodel.onerow <- bind_cols_check_emptyTibble(onemodel.onerow, 
                                                 onemodel.glance.onerow)
  
  # add a column of element ids:
  colnames.temp <- colnames(onemodel.onerow)
  onemodel.onerow <- onemodel.onerow %>% tibble::add_column(element_id = i_element-1, .before = colnames.temp[1])   # add as the first column
  
  # change from tibble to a data.frame:
  onemodel.onerow <- onemodel.onerow %>% as.data.frame()  
  
  onemodel.onerow
  
}

#' Generate expected results for lm
#' 
#' @param fn.phenotypes A character, the filename of the csv file, usually the n50_cohort.csv
#' @param fn.h5 A character, the filename of the h5 file, usually the n50_fixels.h5
#' @param idx.fixel.lm An integer (>=1), the fixel index that is used for this lm testing
#' @param num.set.seed An integer, which will be used by `set.seed()`
#' @return expected.results A list of the expected results
helper_generate_expect_lm <- function(fn.phenotypes,
                                      fn.h5,
                                      idx.fixel.lm = 1,
                                      num.set.seed = 5) {
  ## basic setups:
  phenotypes <- read.csv(fn.phenotypes)
  nsubj <- nrow(phenotypes)
  # ordered factor:
  phenotypes$oSex <- ordered(phenotypes$sex, levels = c("F", "M"))  # ordered factor, "F" as reference group
  
  # set the seed:
  set.seed(num.set.seed)
  
  
  ## start to calculate:
  expected.results <- list()
  
  # load data for lm:
  # not to load from .txt - as it'll be off with that in .h5 (2e-6), making results a little off too (e.g. <2e-6), so cannot tell if the results are perfectly matched or not
  # fd.simu <- read.csv(paste0("data/n50_fixels_FD_idx-", toString(idx.fixel.lm), ".txt"))
  # fd.simu <- fd.simu$x
  
  h5f <- rhdf5::H5Fopen(fn.h5)
  h5d <- h5f$scalars$FD$values  # enter the dataset
  
  fd.simu <- h5d[idx.fixel.lm,]
  
  h5closeAll()
  
  
  data <- phenotypes
  data$FD <- fd.simu
  
  # different scenario:
  thename <- "age"
  formula <- FD ~ age
  dfout <- calcu_stat_lm(formula, data, idx.fixel.lm)
  expected.results[[thename]] <- dfout
  
  thename <- "age_sex"
  formula <- FD ~ age + sex
  dfout <- calcu_stat_lm(formula, data, idx.fixel.lm)
  expected.results[[thename]] <- dfout
  
  thename <- "age_factorA"
  formula <- FD ~ age + factorA
  dfout <- calcu_stat_lm(formula, data, idx.fixel.lm)
  expected.results[[thename]] <- dfout
  
  thename <- "age_factorB"
  formula <- FD ~ age + factorB
  dfout <- calcu_stat_lm(formula, data, idx.fixel.lm)
  expected.results[[thename]] <- dfout
  
  # phenotypes with NA:
  phenotypes_wNA <- phenotypes
  phenotypes_wNA$age[1] <- NA
  data_wNA <- phenotypes_wNA
  data_wNA$FD <- fd.simu
  
  thename <- "age_phenotypeswNA_na.action-na.omit"
  formula <- FD ~ age
  dfout <- calcu_stat_lm(formula, data_wNA, idx.fixel.lm, 
                         na.action = "na.omit")
  expected.results[[thename]] <- dfout
  
  # weights in lm() is different:
  weights.random <- abs(rnorm(nrow(phenotypes)))  # this requires `set.seed()`
  thename <- "age_weights-random"
  formula <- FD ~ age
  dfout <- calcu_stat_lm(formula, data, idx.fixel.lm, 
                         weights = weights.random)
  expected.results[[thename]] <- dfout
  
  # return:
  return(expected.results)
}

#' Generate expected results for gam
#' 
#' @param fn.phenotypes A character, the filename of the csv file, usually the n50_cohort.csv
#' @param fn.h5 A character, the filename of the h5 file, usually the n50_fixels.h5
#' @param idx.fixel.gam An integer (>=1), the fixel index that is used for this gam testing
#' @param num.set.seed An integer, which will be used by `set.seed()`
#' @return expected.results A list of the expected results
helper_generate_expect_gam <- function(fn.phenotypes,
                                       fn.h5,
                                       idx.fixel.gam = 11,
                                       num.set.seed = 5) {
  ## basic setups:
  phenotypes <- read.csv(fn.phenotypes)
  nsubj <- nrow(phenotypes)
  # ordered factor:
  phenotypes$oSex <- ordered(phenotypes$sex, levels = c("F", "M"))  # ordered factor, "F" as reference group
  
  # multiple levels:
  phenotypes$oMultiLevels <- c(rep("A",15),
                               rep("B",15),
                               rep("C",20))
  phenotypes$oMultiLevels <- ordered(phenotypes$oMultiLevels, levels = c("A","B","C"))  # ordered factor, "A" as reference group
  
  # set the seed:
  set.seed(num.set.seed)
  
  
  ## start to calculate:
  expected.results <- list()
  
  # load data for gam:
  h5f <- rhdf5::H5Fopen(fn.h5)
  h5d <- h5f$scalars$FD$values  # enter the dataset
  
  fd.simu <- h5d[idx.fixel.gam,]
  
  h5closeAll()
  
  data <- phenotypes
  data$FD <- fd.simu
  
  # different scenario:
  thename <- "s-age_sex"
  formula <- FD ~ s(age) + sex
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age"
  formula <- FD ~ s(age)
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "age"
  formula <- FD ~ age
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)   # compared to lm()'s results with same formula & data, random selected stat; 2022.5.23. Chenying
  expected.results[[thename]] <- dfout
  
  thename <- "s-age_s-factorA"
  formula <- FD ~ s(age) + s(factorA)
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age-k-4_sex"
  formula <- FD ~ s(age, k=4) + sex
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age-fx-T_sex"
  formula <- FD ~ s(age, fx=TRUE) + sex
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age-bs-cr_sex"
  formula <- FD ~ s(age, bs="cr") + sex
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age_sex_method-REML"
  formula <- FD ~ s(age) + sex
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam, method="REML")   # MAKE SURE THIS IS IN! TODO
  expected.results[[thename]] <- dfout
  
  thename <- "factorB_s-age_s-factorA"
  formula <- FD ~ factorB + s(age) + s(factorA)
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "factorB_s-factorA"
  formula <- FD ~ factorB + s(factorA)
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "factorB_s-age"
  formula <- FD ~ factorB + s(age)
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age-k-4"
  formula <- FD ~ s(age, k=4)
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "te-age"
  formula <- FD ~ te(age)
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age-factorA-fx-F-bs-tpcr"
  formula <- FD ~ s(age, factorA, fx = FALSE, bs = c("tp", "cr")) 
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age-factorA-k-4-bs-tptp"
  formula <- FD ~ s(age, factorA, k=4, bs = c("tp", "tp"))
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "ti-age-fx-F-bs-cr"
  formula <- FD ~ ti(age, fx = FALSE, bs = c("cr"))
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "ti-age-factorA-fx-T-bs-crtp"
  formula <- FD ~ ti(age, factorA, fx = TRUE, bs = c("cr", "tp")) 
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "oSex_s-age-k-4-fx-T_s-age-byoSex-fx-T_factorB"
  formula <- FD ~ oSex + s(age,k=4, fx=TRUE) + s(age, by=oSex, fx=TRUE) + factorB
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "oMultiLevels_s-age-k-4-fx-T_s-age-byoMultiLevels-fx-T_factorB"
  formula <- FD ~ oMultiLevels + s(age,k=4, fx=TRUE) + s(age, by=oMultiLevels, fx=TRUE) + factorB 
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "oSex_s-age-k-4-fx-T_factorB"
  formula <- FD ~ oSex + s(age,k=4, fx=TRUE) + factorB
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "oSex_s-age-byoSex-fx-T_factorB"
  formula <- FD ~ oSex + s(age, by=oSex, fx=TRUE) + factorB
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "s-age-k-4-fx-T_s-age-byoSex-fx-T_factorB"
  formula <- FD ~ s(age,k=4, fx=TRUE) + s(age, by=oSex, fx=TRUE) + factorB
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "ti-age-fx-T_ti-factorB-fx-T_ti-age-factorB-fx-T_factorA"
  formula <- FD ~ ti(age, fx=TRUE) + ti(factorB, fx=TRUE) + ti(age, factorB, fx=TRUE) + factorA
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "ti-age-fx-T_ti-factorB-fx-T_factorA"
  formula <- FD ~ ti(age, fx=TRUE) + ti(factorB, fx=TRUE)+ factorA
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  thename <- "intercept"   # not to replicate the process here.... otherwise dfout's class is "tbl_df"...
  formula <- FD ~ 1
  dfout <- calcu_stat_gam(formula, data, idx.fixel.gam)
  expected.results[[thename]] <- dfout
  
  return(expected.results)
}



#' Generate expected values for accessors, without using ModelArray functions
#' 
#' @param fn.h5 A character, the filename of the h5 file
#' @return the matrix of the scalar values
helper_generate_expect_accessors <- function(fn.h5) {
  # load the data:
  h5f <- rhdf5::H5Fopen(fn.h5)
  scalar.value.full <- h5f$scalars$FD$values  # enter the dataset
  # ^ has been realized as matrix, taking ~70MB of memory if using n50_fixels.h5 file
  
  h5closeAll()
  
  return(scalar.value.full)
  
  
}