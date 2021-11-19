#' print the additional arguments settings
#' @param FUN The function, e.g. mgcv::gam, without "()"
#' @param argu_name The argument name of the function
#' @param dots: list of additional arguments
#' @param message_default The message for default 
#' @param message_usr_input The message describing user's input

printAdditionalArgu <- function(FUN, argu_name, dots, message_default = NULL, message_usr_input = NULL) {
  dots_names <- names(dots)
  if (argu_name %in% dots_names) {
    if (is.null(message_default)) {
      message_default <- invisible(eval(formals(FUN)[[argu_name]]))
    }
    if (is.null(message_default)) {  # if the default is NULL:
      message_default <- "NULL"
    } 
    
    if (is.null(message_usr_input)) {
      m1 <- paste0(argu_name, " = ", dots[[argu_name]], " (default: ", message_default, ")") %>% crayon::black() %>% cat()    # or, %>% message()  
    } else {   # specified the message:
      m1 <- paste0(argu_name, " = ", message_usr_input, " (default: ", message_default, ")") %>% crayon::black() %>% cat()
    }
    

  } else {
    m1<-paste0(argu_name, ": default")%>% crayon::black() %>% cat() 
    
  }

  cat(m1, "\n")
}



#' Run a linear model at each fixel location
#'
#' @param formula Formula (passed to `lm()`)
#' @param data FixelArray dataset
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param verbose Print progress messages
#' @param idx A vector of fixel IDs to subset
#' @param pbar Print progress bar
#' @param n_cores The number of cores to run on
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.old.lm <- function(formula, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, ...){
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  ### display additional arguments:
  dots <- list(...)
  dots_names <- names(dots)
  
  FUN <- stats::lm
  
  # subset:
  m1 <- "no default"  
  printAdditionalArgu(FUN, "subset", dots, m1)
  
  # weights:
  m1 <- "no default"
  if ("weights" %in% dots_names) {  # if user provides weights
    m_usr_input <- paste0( class(dots$weights), " with length of ", length(dots$weights)) # message describing usr's input; cannot use dim on c(1,2) 
  } else {
    m_usr_input <- NULL
  }
  printAdditionalArgu(FUN, "weights", dots, m1, m_usr_input)
  
  # na.action:
  m1 <- "no default"
  printAdditionalArgu(FUN, "na.action", dots, m1)
  
  # method:
  printAdditionalArgu(FUN, "method", dots)  # default: "qr"
  
  # model:
  m1 <- invisible(eval(formals(FUN)[["model"]])) %>% as.character()   # default: [logical] TRUE
  printAdditionalArgu(FUN, "model", dots, m1)  
  
  # x:
  m1 <- invisible(eval(formals(FUN)[["x"]])) %>% as.character()   # default: [logical] FALSE
  printAdditionalArgu(FUN, "x", dots, m1)
  
  # y:
  m1 <- invisible(eval(formals(FUN)[["y"]])) %>% as.character()   # default: [logical] FALSE
  printAdditionalArgu(FUN, "y", dots, m1)
  
  # qr:
  m1 <- invisible(eval(formals(FUN)[["qr"]])) %>% as.character()   # default: [logical] TRUE
  printAdditionalArgu(FUN, "qr", dots, m1)
  
  # singular.ok:
  m1 <- invisible(eval(formals(FUN)[["singular.ok"]])) %>% as.character()   # default: [logical] TRUE
  printAdditionalArgu(FUN, "singular.ok", dots, m1)
  
  # contrasts:
  printAdditionalArgu(FUN, "contrasts", dots)   # default: NULL
  
  # offset:
  m1 <- "no default"   # there is no default
  printAdditionalArgu(FUN, "offset", dots, m1)
  
  ### start the process:
  if(verbose){
    message(glue::glue("Fitting fixel-wise linear models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        stats::lm(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, mc.cores = n_cores, ...)
      
    } else {
      
      fits <- parallel::mclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        stats::lm(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
        # NOTES: see current FixelArray.lm, after: remove tibble information and turn into numeric function
        
        
        # lm(formula, data = dat, ...) %>%
        #   broom::glance() %>%
        #   print()
        
      }, mc.cores = n_cores, ...)
      
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        stats::lm(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, ...)
      
    }
    else {
      
      fits <- lapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        stats::lm(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, ...)
      
    }
  }
  
  df_out <- do.call(rbind, fits)
  
  # if(write){
  #   WriteResult(data, df_out, glue::glue("scalars/{scalar}/results"))
  # }
  
  df_out
  
}

#' Run a linear model at each fixel location, write out each result just after the model fitting 
#' For p-value corrections (arguments correct.p.value.*), supported methods include all methods in `p.adjust.methods` except "none". Can be more than one method. Turn it off by setting to "none".
#'
#' @param formula Formula (passed to `lm()`)
#' @param data FixelArray class
#' @param phenotypes The cohort matrix with covariates to be added to the model  
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param fixel.subset The subset of fixel ids you want to run. Integers. First id starts from 1.
#' @param full.outputs Whether to return full set of outputs (TRUE or FALSE). If FALSE, it will only return those listed in var.terms and var.model; if TRUE, arguments var.terms and var.model will be ignored.
#' @param var.terms The list of variables to save for terms (got from lm %>% tidy())
#' @param var.model The list of variables to save for the model (got from lm %>% glance())
#' @param correct.p.value.terms To perform and add a column for p.value correction for each term. 
#' @param correct.p.value.model To perform and add a column for p.value correction for the model. 
#' @param verbose Print verbose message or not
#' @param pbar Print progress bar
#' @param n_cores The number of cores to run on
#' @return Tibble with the summarized model statistics for all fixels requested
#' @import doParallel
#' @import tibble
#' @export

FixelArray.lm <- function(formula, data, phenotypes, scalar, fixel.subset = NULL, full.outputs = FALSE, 
                              var.terms = c("estimate", "statistic", "p.value"), 
                              var.model = c("adj.r.squared", "p.value"), 
                          correct.p.value.terms = "none", correct.p.value.model = "none",
                              verbose = TRUE, pbar = TRUE, n_cores = 1, ...) {
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  # checker for min and max of fixel.subset; and whether elements are integer
  if (min(fixel.subset) < 1) {
    stop("Minimal value in fixel.subset should >= 1")
  }
  if (max(fixel.subset) > nrow(scalars(data)[[scalar]])) {
    stop(paste0("Maximal value in fixel.subset should <= number of fixels = "), as.character(nrow(scalars(data)[[scalar]])))
  }
  if (class(fixel.subset) != "integer") {
    stop("Please enter integers for fixel.subset!")
  }
  
  
  ### display additional arguments:
  dots <- list(...)
  dots_names <- names(dots)
  
  FUN <- stats::lm
  
  # subset:
  m1 <- "no default"  
  printAdditionalArgu(FUN, "subset", dots, m1)
  
  # weights:
  m1 <- "no default"
  if ("weights" %in% dots_names) {  # if user provides weights
    m_usr_input <- paste0( class(dots$weights), " with length of ", length(dots$weights)) # message describing usr's input; cannot use dim on c(1,2) 
  } else {
    m_usr_input <- NULL
  }
  printAdditionalArgu(FUN, "weights", dots, m1, m_usr_input)
  
  # na.action:
  m1 <- "no default"
  printAdditionalArgu(FUN, "na.action", dots, m1)
  
  # method:
  printAdditionalArgu(FUN, "method", dots)  # default: "qr"
  
  # model:
  m1 <- invisible(eval(formals(FUN)[["model"]])) %>% as.character()   # default: [logical] TRUE
  printAdditionalArgu(FUN, "model", dots, m1)  
  
  # x:
  m1 <- invisible(eval(formals(FUN)[["x"]])) %>% as.character()   # default: [logical] FALSE
  printAdditionalArgu(FUN, "x", dots, m1)
  
  # y:
  m1 <- invisible(eval(formals(FUN)[["y"]])) %>% as.character()   # default: [logical] FALSE
  printAdditionalArgu(FUN, "y", dots, m1)
  
  # qr:
  m1 <- invisible(eval(formals(FUN)[["qr"]])) %>% as.character()   # default: [logical] TRUE
  printAdditionalArgu(FUN, "qr", dots, m1)
  
  # singular.ok:
  m1 <- invisible(eval(formals(FUN)[["singular.ok"]])) %>% as.character()   # default: [logical] TRUE
  printAdditionalArgu(FUN, "singular.ok", dots, m1)
  
  # contrasts:
  printAdditionalArgu(FUN, "contrasts", dots)   # default: NULL
  
  # offset:
  m1 <- "no default"   # there is no default
  printAdditionalArgu(FUN, "offset", dots, m1)
  
  
  
  ### other setups:
  if (full.outputs == TRUE) {   # full set of outputs
    var.terms = c("estimate","std.error","statistic","p.value")
    var.model = c("r.squared", "adj.r.squared", "sigma", "statistic", "p.value", "df", "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs")
  }
  
  # check for p.value correction:
  p.adjust.methods.full <- p.adjust.methods[ p.adjust.methods != "none" ]
    # check for terms:
  if ( all(correct.p.value.terms == "none") == FALSE) {    # any element is not "none"
    checker.method.in <- correct.p.value.terms %in% p.adjust.methods.full
    if ( all(checker.method.in) == FALSE) {   # not all "TRUE"
      stop(paste0("Some of elements in correct.p.value.terms are not valid, so not to perform terms's p.value corrections. Valid inputs are: ", paste(p.adjust.methods.full, collapse = ',')))
    } 
    
    if ("p.value" %in% var.terms == FALSE) {  # not in the list | # check whether there is "p.value" in var.terms'
      warning(paste0("p.value was not included in var.terms, so not to perform terms's p.value corrections"))
    }
  }
  
    # check for model:
  if ( all(correct.p.value.model == "none") == FALSE ) {    # any element is not "none"
    checker.method.in <- correct.p.value.model %in% p.adjust.methods.full
    if ( all(checker.method.in) == FALSE) {   # not all "TRUE"
      stop(paste0("Some of elements in correct.p.value.model are not valid, so not to perform model's p.value corrections. Valid inputs are: ", paste(p.adjust.methods.full, collapse = ',')))
    } 
    
    if ("p.value" %in% var.model == FALSE) {  # not in the list | # check whether there is "p.value" in var.model'
      warning(paste0("p.value was not included in var.model, so not to perform model's p.value corrections"))
    }
  }
  

  
  ### start the process:
  if(verbose){
    message(glue::glue("Fitting fixel-wise linear models for {scalar}", ))
    message(glue::glue("initiating....", ))
  }
  

  
  # initiate: get the example of one fixel and get the column names
  outputs_initiator <- analyseOneFixel.lm(i_fixel=1, formula, data, phenotypes, scalar, 
                                     var.terms, var.model, 
                                     flag_initiate = TRUE, 
                                     ...)
  column_names <- outputs_initiator$column_names
  list.terms <- outputs_initiator$list.terms

  # loop (by condition of pbar and n_cores)
  if(verbose){
    message(glue::glue("looping across fixels....", ))
  }
  
  # is it a multicore process?
  flag_initiate <- FALSE
  if(n_cores > 1){
    
    if (pbar) {
      
      fits <- pbmcapply::pbmclapply(fixel.subset,   # a list of i_fixel
                                    analyseOneFixel.lm,  # the function
                                    mc.cores = n_cores,
                                    formula, data, phenotypes, scalar,
                                    var.terms, var.model,
                                    flag_initiate = FALSE,
                                    ...)
      
    } else {
      
      # foreach::foreach
      
      fits <- parallel::mclapply(fixel.subset,   # a list of i_fixel 
                                 analyseOneFixel.lm,  # the function
                                 mc.cores = n_cores,
                                 formula, data, phenotypes, scalar,
                                 var.terms, var.model,
                                 flag_initiate = FALSE,
                                 ...)
      
    }
  } else  {  # n_cores ==1, not multi-core
    
    if (pbar) {
      
      fits <- pbapply::pblapply(fixel.subset,   # a list of i_fixel
                                analyseOneFixel.lm,  # the function
                                formula, data, phenotypes, scalar,
                                var.terms, var.model,
                                flag_initiate = FALSE,
                                ...)
      
    } else {
      
      fits <- lapply(fixel.subset,   # a list of i_fixel
                     analyseOneFixel.lm,  # the function
                     formula, data, phenotypes, scalar,
                     var.terms, var.model,
                     flag_initiate = FALSE,
                     ...)
    }
  }  
    

  
  df_out <- do.call(rbind, fits)    
  df_out <- as.data.frame(df_out)    # turn into data.frame
  colnames(df_out) <- column_names     # add column names
  

  # Add corrections of p.values:
  
    # loop over elements in correct.p.value.model
      # if == "none": do nothing
        # else, %in% default methods in p.adjust --> all TRUE? if not, error: not support | checker at beginning of this function
          # else, "p.value" %in% var.model == FALSE: warning: nothing to correct | checker at beginning of this function
            # else, iterate
  
  # add correction of p.values: for terms
  if ( all(correct.p.value.terms == "none") ) {    # all() is to accormodate for multiple elements in correct.p.value.terms: if one of is not "none", FALSE
    # do nothing
    
  } else {
    if ("p.value" %in% var.terms == TRUE) {   # check whether there is "p.value" in var.terms' | if FALSE: print warning (see beginning of this function)
      
      for (methodstr in correct.p.value.terms) {
        
        for (tempstr in list.terms) {
          tempstr.raw <- paste0(tempstr, ".p.value")
          tempstr.corrected <- paste0(tempstr.raw, ".", methodstr)
          temp.corrected <- p.adjust(df_out[[tempstr.raw]], method = methodstr)
          df_out <- df_out %>% tibble::add_column( "{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
        }
        
      }
      
    }
    
  }
  

  # if ( correct.p.value.terms %in% c("fdr", "bonferroni") ) { 
  #   
  #   if ("p.value" %in% var.terms == FALSE) {  # not in the list | # check whether there is "p.value" in var.terms
  #     warning(paste0("p.value was not included in var.terms, so not to perform terms' p.value corrections"))
  #     
  #   } else {
  #     
  #     for (tempstr in list.terms) {
  #       tempstr.raw <- paste0(tempstr, ".p.value")
  #       tempstr.corrected <- paste0(tempstr.raw, ".", correct.p.value.terms)
  #       temp.corrected <- p.adjust(df_out[[tempstr.raw]], method = correct.p.value.terms)
  #       df_out <- df_out %>% tibble::add_column( "{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
  #     }
  #   }
  #   
  # } else if (correct.p.value.terms == "none") {
  #   # do nothing
  # } else {
  #   warning(paste0("correct.p.value.terms = ", toString(correct.p.value.terms), " is not supported!"))
  # }
  
  
  # add correction of p.values: for the model
  if (  all(correct.p.value.model == "none") ) {
    # do nothing
    
  } else {
    if ("p.value" %in% var.model == TRUE) {   # check whether there is "p.value" in var.model' | if FALSE: print warning (see beginning of this function)

      for (methodstr in correct.p.value.model) {
        
        tempstr.raw <- "model.p.value"
        tempstr.corrected <- paste0(tempstr.raw, ".", methodstr)
        temp.corrected <- p.adjust(df_out[[tempstr.raw]], method = methodstr)
        df_out <- df_out %>% tibble::add_column( "{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
        
      }
        
    }
      
  }
 
  # 
  # if (correct.p.value.model %in% c("fdr", "bonferroni")) {
  #   
  #   if ("p.value" %in% var.model == FALSE) {  # not in the list | # check whether there is "p.value" in var.model
  #     warning(paste0("p.value was not included in var.model, so not to perform model's p.value corrections"))
  #     
  #   } else {
  #     
  #     tempstr.raw <- "model.p.value"
  #     tempstr.corrected <- paste0(tempstr.raw, ".", correct.p.value.model)
  #     temp.corrected <- p.adjust(df_out[[tempstr.raw]], method = correct.p.value.model)
  #     df_out <- df_out %>% tibble::add_column( "{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
  #     
  #   }
  #   
  # } else if (correct.p.value.model == "none") {
  #   # do nothing
  # } else {
  #   warning(paste0('correct.p.value.model = "', toString(correct.p.value.model), '" is not supported!'))
  # }
  
  
  
  df_out   # return

}



#' Run a GAM model at each fixel location
#' For p-value corrections (arguments correct.p.value.*), supported methods include all methods in `p.adjust.methods` except "none". Can be more than one method. Turn it off by setting to "none".
#' @param formula Formula (passed to `mgcv::gam()`)
#' @param data FixelArray class
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param fixel.subset The subset of fixel ids you want to run. Integers. First id starts from 1.
#' @param full.outputs Whether to return full set of outputs (TRUE or FALSE). If FALSE, it will only return those listed in var.terms and var.model; if TRUE, arguments var.terms and var.model will be ignored.
#' @param var.smoothTerms The list of variables to save for smooth terms (got from gam %>% tidy(parametric = FALSE)). Example smooth term: age in formula "outcome ~ s(age)".
#' @param var.parametricTerms The list of variables to save for parametric terms (got from gam %>% tidy(parametric = TRUE)). Example parametric term: sex in formula "outcome ~ s(age) + sex"
#' @param var.model The list of variables to save for the model (got from lm %>% glance())
#' @param correct.p.value.smoothTerms To perform and add a column for p.value correction for each smooth term. 
#' @param correct.p.value.parametricTerms To perform and add a column for p.value correction for each parametric term. 
#' @param correct.p.value.model To perform and add a column for p.value correction for the model. 
#' @param verbose Print progress messages
#' @param pbar Print progress bar
#' @param n_cores The number of cores to run on
#' @return Tibble with the summarized model statistics for all fixels requested
#' @import doParallel
#' @import tibble
#' @export

FixelArray.gam <- function(formula, data, phenotypes, scalar, fixel.subset = NULL, full.outputs = FALSE, 
                              var.parametricTerms = c("statistic","p.value"),
                              var.smoothTerms = c("estimate", "statistic", "p.value"),
                              var.model = c("adj.r.squared", "p.value"), 
                              correct.p.value.smoothTerms = "none", correct.p.value.parametricTerms = "none", correct.p.value.model = "none",
                              verbose = TRUE, pbar = TRUE, n_cores = 1, ...){
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  
  if(verbose){
    message(glue::glue("Fitting fixel-wise GAM for {scalar}", ))
  }
  
  
  # turn to explicit ones, for checking input arguments are valid
  dots <- list(...)    
  dots_names <- names(dots)
  
  # if ("family" %in% dots_names) {
  #   m <- print(eval(formals(mgcv::gam)$family))
  #   message(paste0("family = ", dots$family, " (default: ", "Family: ", m$family, "; Link function: ", m$link, ")"))    #  eval(formals(mgcv::gam)$family) %>% print()
  # } else {
  #   message("family: default")
  # }
  
  FUN <- mgcv::gam
  
  m <- invisible(eval(formals(FUN)$family))    # should not use message(), but print() --> but will print out or invisible()
  m1 <- paste0("Family: ", m$family, "; Link function: ", m$link)
  printAdditionalArgu(FUN, "family", dots, m1)
  
  
  # eval(formals(mgcv::gam)$data) # return the default setting of argument "data"
  
  
  
  
}



#' Run a t.test at each fixel location
#'
#' @param formula Formula (passed to `lm()`)
#' @param data FixelArray dataset
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param subset A vector of fixel IDs to subset
#' @param verbose Print progress messages
#' @param n_cores The number of cores to run on
#' @param pbar Print progress bar
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.t.test <- function(formula, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, ...){
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  
  if(verbose){
    message(glue::glue("Running t.test models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        t.test(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, mc.cores = n_cores, ...)
      
    } else {
      
      # if it's Linux or Mac:
      system_name <- Sys.info()['sysname']
      if (system_name == "Linux" || system_name == "Darwin") { # Darwin = OSX = Mac

        fits <- parallel::mclapply(ids, function(i, ...){
          
          values <- scalars(data)[[scalar]][i,]
          # values <- data@scalars[[scalar]][i,]
          dat <- phenotypes
          dat[[scalar]] <- values
          
          t.test(formula, data = dat, ...) %>%
            broom::tidy() %>%
            dplyr::mutate(fixel_id = i-1)
          
          
        }, mc.cores = n_cores, ...)
      } else if (system_name == "Windows") {   # 7/27/2021: there is still error for Windows system.... current error: see below
        # cl <- makeCluster(getOption("cl.cores", n_cores))
        
        # message(slotNames(data))  # identified slots for data here..
        
        cl <- makeCluster(n_cores)
        
        clusterEvalQ(cl, {"FixelArray"}) 
        
        fits <- parallel::parLapply(cl, ids, function(i, ...){
          # source("FixelArray.R")
          # values <- scalars(data)[[scalar]][i,]   # ERROR: could not find function "scalars"
          values <- data@scalars[[scalar]][i,]   # ERROR: object of type 'S4' is not subsettable
          dat <- phenotypes
          dat[[scalar]] <- values
          
          t.test(formula, data = dat, ...) %>%
            broom::tidy() %>%
            dplyr::mutate(fixel_id = i-1)
          
          
        }, mc.cores = n_cores, ...)
        
        stopCluster(cl)   # for windows
      }
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        t.test(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, ...)
      
    }
    else {
      
      fits <- lapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        t.test(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, ...)
      
    }
  }
  
  df_out <- do.call(rbind, fits)
  
  # if(write){
  #   WriteResult(data, df_out, glue::glue("scalars/{scalar}/results"))
  # }
  
  df_out
  
}

#' Run a GAMM4 model at each fixel location
#'
#' @param formula Formula (passed to `gamm4()`)
#' @param data FixelArray dataset
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param subset A vector of fixel IDs to subset
#' @param verbose Print progress messages
#' @param n_cores The number of cores to run on
#' @param pbar Print progress bar
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.gamm4 <- function(formula, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, ...){
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  
  if(verbose){
    message(glue::glue("Fitting fixel-wise linear models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        gamm4::gamm4(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, mc.cores = n_cores, ...)
      
    } else {
      
      fits <- parallel::mclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        gamm4::gamm4(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
        
      }, mc.cores = n_cores, ...)
      
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        gamm4::gamm4(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, ...)
      
    }
    else {
      
      fits <- lapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        gamm4::gamm4(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, ...)
      
    }
  }
  
  df_out <- do.call(rbind, fits)
  
  # if(write){
  #   WriteResult(data, df_out, glue::glue("scalars/{scalar}/results"))
  # }
  
  df_out
  
}



#' Run a user-defined function on each fixel location
#'
#' @param formula Formula (passed to `lm()`)
#' @param FUN User-defined modelling function; must return a 1-row named vector or data.frame
#' @param data FixelArray dataset
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param subset A vector of fixel IDs to subset
#' @param verbose Print progress messages
#' @param n_cores The number of cores to run on
#' @param pbar Print progress bar
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.model <- function(formula, FUN, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, ...){
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  

  FUN <- match.fun(FUN)
  
  if(verbose){
    message(glue::glue("Fitting fixel-wise linear models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(x, ...){
        
        values <- scalars(data)[[scalar]][x,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        FUN(formula, data = dat, ...) %>%
          broom::tidy()
        
      }, mc.cores = n_cores, ...)
      
    } else {
      
      fits <- parallel::mclapply(ids, function(x, ...){
        
        values <- scalars(data)[[scalar]][x,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        FUN(formula, data = dat, ...) %>%
          broom::tidy()
        
        
      }, mc.cores = n_cores, ...)
      
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(x, ...){
        
        values <- scalars(data)[[scalar]][x,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        FUN(formula, data = dat, ...) %>%
          broom::tidy()
        
      }, ...)
      
    }
    else {
      
      fits <- lapply(ids, function(x, ...){
        
        values <- scalars(data)[[scalar]][x,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        FUN(formula, data = dat, ...) %>%
          broom::tidy()
        
      }, ...)
      
    }
  }
  
  return(do.call(rbind, fits))
  
}