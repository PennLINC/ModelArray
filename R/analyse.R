#' print the additional arguments settings
#' @param FUN The function, e.g. mgcv::gam, without "()"
#' @param argu_name The argument name of the function
#' @param dots: list of additional arguments
#' @param message_default The message for default 
#' @param message_usr_input The message describing user's input
#' @import crayon
#' @import dplyr
#' 
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


#' Check if the list of p-value correction methods are valid for a specific type of term/model. 
#' Can be used for any statistical model. As long as the p.value to be correct is named as "p.value".
#' 
#' @param correct.list The list of correction methods for this type of term/model
#' @param name.correct.list The name of the list of correction methods for this type of term/model
#' @param var.list The list of statistics to be saved for this type of term/model
#' @param name.var.list The name of the list of statistics to be saved for this type of term/model

check_validity_correctPValue <- function(correct.list, name.correct.list, 
                                         var.list, name.var.list) {
  p.adjust.methods.full <- p.adjust.methods[ p.adjust.methods != "none" ]
  if ( all(correct.list == "none") == FALSE) {    # any element is not "none"
    checker.method.in <- correct.list %in% p.adjust.methods.full
    if ( all(checker.method.in) == FALSE) {   # not all "TRUE"
      stop(paste0("Some of elements in ",name.correct.list," are not valid. Valid inputs are: ", paste(p.adjust.methods.full, collapse = ', ')))
    } 
    
    if ("p.value" %in% var.list == FALSE) {  # not in the list | # check whether there is "p.value" in var.list
      warning(paste0("p.value was not included in ",name.var.list,", so not to perform its p.value corrections"))   # TODO: why this warning comes out after FixelArray.aModel is done?
    }
  }
}


#' Print out important arguments in smooth terms s() in mgcv::gam() formula
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/s
#' TODO: finish the description
#' @param ofInterest got via: gam.formula.breakdown <- mgcv::interpret.gam(formula); ofInterest <- gam.formula.breakdown$smooth.spec[[i]]
#' @import mgcv
#' @import dplyr
#' @import crayon
#' 
checker_gam_s <- function(ofInterest) {
  FUN <- mgcv::s
  
  paste0(ofInterest$label, ": ") %>% crayon::black() %>% cat()
  
  ### k (or bs.dim):   # could be multiple values
  m1 <- invisible(eval(formals(FUN)[["k"]]))  # default
  m2 <- ofInterest$bs.dim   # could be a list of multiple values
  if ((length(unique(m2)) == 1) & (m1 %in% unique(m2))) {  # default
    msg_k <- " (default)"
  } else {
    msg_k <- ""
  }
  
  paste0("  k = ", 
         paste(as.character(m2), collapse = ", "), msg_k, "; ") %>% crayon::black() %>% cat()
  
  ### fx:
  m1 <- invisible(eval(formals(FUN)[["fx"]])) %>% as.character()  # default
  m2 <- ofInterest$fixed   # actual
  
  if (as.character(!as.logical(m1)) %in% m2) {  # there is an opposite logical value in m2 (actual)
    msg_fx <- ""
  } else {
    msg_fx <- " (default)"
  }
  
  paste0("  fx = ", toString(ofInterest$fixed), msg_fx, "; ") %>% crayon::black() %>% cat()
  
  ### bs: 
  m1 <- invisible(eval(formals(FUN)[["bs"]]))   # default
  # actual:
  mybs <- gsub(".smooth.spec", "",class(ofInterest))  # if there are multiple elements in bs, mybs will be a list
  if ((length(unique(mybs)) == 1) & (m1 %in% unique(mybs))) {  # default
    msg_bs <- " (default)"
  } else {
    msg_bs <- ""
  }
  
  paste0("  bs = ",
         paste(as.character(mybs), collapse = ", "), msg_bs) %>% crayon::black() %>% cat()
  
  ### TODO: add "by" - interaction terms | if the factor after by is   not listed in the "label"
  
  
  cat("\n")
  

}

#' Print out important arguments in smooth term te() or ti() or t2() in mgcv::gam() formula
#' Why a separate function is needed for t(), cannot using s(): in ofInterest, "fx" is "fx" for t(), but "fixed" for s() - so they are different.
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/te or /t2()
#' TODO: finish the description
#' @param FUN could be mgcv::te(), ti() or t2()
#' @param ofInterest got via: gam.formula.breakdown <- mgcv::interpret.gam(formula); ofInterest <- gam.formula.breakdown$smooth.spec[[i]]
#' @import mgcv
#' @import dplyr
#' @import crayon
#' 
checker_gam_t <- function(FUN, ofInterest) {
  paste0(ofInterest$label, ": ") %>% crayon::black() %>% cat()
  
  # t()'s interpret.gam()'s smooth.spec does not have k or bs.dim as s() does; may provided in ofInterest$margin[[i-xx]]$bs.dim but not fully clear

  ### fx:
  m1 <- invisible(eval(formals(FUN)[["fx"]])) %>% as.character()  # default
  m2 <- ofInterest$fx   # actual

  if (as.character(!as.logical(m1)) %in% m2) {  # there is an opposite logical value in m2 (actual)
    msg_fx <- ""
  } else {
    msg_fx <- " (default)"
  }
  
  paste0("  fx = ", toString(ofInterest$fx), msg_fx, "; ") %>% crayon::black() %>% cat()
  
  ### bs:  # also different way of extracting from s()'s
  m1 <- invisible(eval(formals(FUN)[["bs"]])) %>% as.character()   # default

  mybs <- list()  # actual, as a list
  for (i in 1:length(ofInterest$margin)) {
    temp <- gsub(".smooth.spec", "", ofInterest$margin[[i]] %>% class() )
    mybs[i] <- temp
  }

  # then check if all elements are default value of bs
  if ((length(unique(mybs)) == 1) & (m1 %in% unique(mybs))) {   # all elements are the same, and = default
    msg_bs <- " (default)"
  } else {
    msg_bs <- ""
  }
  # print out
  paste0("  bs = ", 
         paste(as.character(mybs), collapse = ", "), msg_bs) %>% crayon::black() %>% cat()
  
  ### TODO: add "by" - interaction terms | first check if "by" is accepted by t()
  
  cat("\n")
}

#' A checker for formula in gam for FixelArray.gam()
#' TODO: finish the description
#' @import mgcv
#' @import dplyr
#' @import rlang
checker_gam_formula <- function(formula, gam.formula.breakdown, onemodel) {
  if (length(gam.formula.breakdown$smooth.spec) != 0) {   # if there is smooth term
    list_smooth_terms <- character(length(gam.formula.breakdown$smooth.spec))
    for (i_smoothTerm in 1:length(gam.formula.breakdown$smooth.spec)) {
      ofInterest <- gam.formula.breakdown$smooth.spec[[i_smoothTerm]]
      list_smooth_terms[i_smoothTerm] <- ofInterest$label
      
      # check what class of smooth term (s or ti or ?); then call the corresponding function - throw out the message for important argument for this class
      # ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/smooth.terms
      smooth.class <- strsplit(ofInterest$label, "[(]")[[1]][1]
      if (smooth.class == "s") {
        checker_gam_s(ofInterest)
      } else if (smooth.class == "te") {
        checker_gam_t(mgcv::te, ofInterest)
      } else if (smooth.class == "ti") {
        checker_gam_t(mgcv::ti, ofInterest)
      } else if (smooth.class == "t2") {
        checker_gam_t(mgcv::t2, ofInterest)
      } else {
        stop(paste0("invalid smooth class for term ", ofInterest$label))
      }

    }
  } else {   # no smooth term
    message("Warning: there is no smooth term in the requested formula")
    
  }
  
  ### duplicated terms: 
  # duplicated terms: s(age*factorA)+s(age) --> two s(age); s(age+factorA) + s(age) --> two s(age)
  # use gam.formula.breakdown; compare the smooth term names
  if (length(gam.formula.breakdown$smooth.spec) != 0) {   # if there is smooth term
    if (length(unique(list_smooth_terms)) < length(list_smooth_terms)) {  # there is duplicated one
      stop(paste0("There are duplicated smooth terms: ", paste(list_smooth_terms, collapse = ", ")))
      # tip: s(age*factorA)+s(age); s(age+factorA) + s(age) can also cause this error
    }
  }
  
  # TODO: duplicated terms: s(age) + ti(age); 
  
  ### if there is interactions in smooth term, check which formula it belongs to
  invalid_patternInteract_in_s <- c("*", "+",",")   # invalid: s(a*b), s(a+b), s(a,b)  # however s(a*b) or s(a+b) will become s(a) after brooming...
  invalid_patternInteract_in_t <- c("*", "+",":")   # invalid: t(a*b), t(a+b), t(a, by=b) --> t(a):b
  invalid_smoothClass_interact <- c("te", "t2")
  
  valid_patternInteract_s <- c(":")   # valid: s(a):b <-- got from s(a,by=b)
  valid_patternInteract_t <- c(",")   # valid: t(a,b)
  
  ## check via formula as a string
  # any "*":
  str_var <- as.character(formula)[3]
  if (grepl("*", str_var, fixed=TRUE)) {  # any "*" in the RHS of the formula
    stop("The right hand side of formula should not contain *")
  }
  #strsplit(str_var, "[()]")# split by ( and )  # however not 
  # TODO: any "+":
  for (smooth.class in c("s","ti","te","t2")) {
    # 1. find the smooth term from str_var, including smooth.class   # trying out str_match and strsplit.....
    # 2. fun_call <- quo(<smooth_term>)   # e.g. mgcv::s(age+factorA, k=4)
    # 3. user_args <- call_args(call_standardise(fun_call))  # and get the string for the first argument, age + factorA
  }
  
  ## check via term name after model fitting:
  # fit for one fixel, get the summarized stat:
  onemodel.tidy.smoothTerms <- onemodel %>% broom::tidy(parametric = FALSE)
  onemodel.tidy.parametricTerms <- onemodel %>% broom::tidy(parametric = TRUE)
  onemodel.glance <- onemodel %>% broom::glance()
  onemodel.summary <- onemodel %>% summary()
  
  if (nrow(onemodel.tidy.smoothTerms) != 0) {   # there is smooth term(s)
    list_smoothTerms_name <- onemodel.tidy.smoothTerms$term   # list of names of smooth terms
    
    num.interact.term <- 0   # any finding of more than one smooth term has interaction term; or a smooth term contains more than two variable to be interacted with   # e.g. two ","
    for (i_smoothTerm in 1:length(list_smoothTerms_name)) {
      smoothTerm_name <- list_smoothTerms_name[i_smoothTerm]
      smooth.class <- strsplit(smoothTerm_name, "[(]")[[1]][1]
        
      print(smoothTerm_name)
      
      ## check if there is any sign of WRONG pattern of interaction - via term name after model fitting:
      if (smooth.class == "s") {
        for (invalid_pattern in invalid_patternInteract_in_s) {
          if (grepl(invalid_pattern, smoothTerm_name, fixed = TRUE)) {   # contains the wrong pattern of interaction
            stop(paste0(smoothTerm_name, " contains the invalid interaction pattern for ", smooth.class, "(): one of ", paste(invalid_patternInteract_in_s, collapse = ', ')))
          }  
        }
        # now, it's s() and there is no invalid pattern of interaction; check if there is valid one:
        for (valid_pattern in valid_patternInteract_s) {
          num.interact.term <- num.interact.term + 1
        }
      } else if (smooth.class == "ti") {
        for (invalid_pattern in invalid_patternInteract_in_t) {
          if (grepl(invalid_pattern, smoothTerm_name, fixed = TRUE)) {   # contains the wrong pattern of interaction
            stop(paste0(smoothTerm_name, " contains the invalid interaction pattern for ", smooth.class, "(): one of ", paste(invalid_patternInteract_in_t, collapse = ', ')))
          }  
        }
        # now, it's t() and there is no invalid pattern of interaction; check if there is valid one:
        for (valid_pattern in valid_patternInteract_t) {
          num.interact.term <- num.interact.term + 1
        }
      }
      # after above checking, there is no WRONG pattern; but 
      
      
      
      if (num.interact.term >1) {
        stop("there is more than one interaction term! Either there is more than one smooth term with interaction component, or there is smooth term that contains more than two variables to be interacted with")
      }
    }
    
  } else {
    num.interact.term <- 0
  }
  
  
}

#' Generate GAM formula with interaction term: factor-smooth interaction
#' Example 
#' `factor.var` and `smooth.var` should come from data.frame `phenotypes`
#' TODO: finish description
#' @param response.var character class, the variable name for response
#' @param factor.var character class, the variable name for factor. It should be an ordered factor. If not, it will generate it as a new column in `phenotypes`, which requires `reference.group`.
#' @param smooth.var character class, the variable name in smooth term as main effect
#' @param phenotypes data.frame class, the cohort matrix with covariates to be added to the model 
#' @param reference.group character class, the reference group for ordered factor of `factor.var`; required when `factor.var` in `phenotypes` is not an ordered factor. 
#' @param prefix.ordered.factor character class, the prefix for ordered factor; required when `factor.var` in `phenotypes` is not an ordered factor.
#' @param fx TRUE or FALSE, to be used in smooth term s(). Recommend TRUE.
#' @param k integer, to be used in smooth term s(). Default is -1 as in mgcv::s()
#' @return a list, including: 1) formula generated; 2) data.frame phenotypes - updated if argument factor.var is not an ordered factor
#' @import mgcv
#' 
generator_gamFormula_factorXsmooth <- function(response.var, factor.var, smooth.var, phenotypes, 
                                               reference.group = NULL, prefix.ordered.factor = "o",
                                               fx=TRUE, k=-1) {
  class.factor.var <- class(phenotypes[[factor.var]])
  if (  !( (length(class.factor.var) == 2) & (class.factor.var[1] == "ordered") & (class.factor.var[2] == "factor")  )  ) {   # class is not c("ordered", "factor")

    message("input `factor.var` is not an ordered factor; will generate one in `phenotypes` which will be returned")
    if (is.null(reference.group)) {
      stop("requires a reference.group to generate the ordered factor")
    }
    
    # name of the ordered factor:
    unordered.factor.var <- factor.var
    factor.var <- paste0(prefix.ordered.factor, unordered.factor.var)
    
    message(paste0("the ordered factor will be named as: ", factor.var))
    
    # check if factor.var already exists:
    if (factor.var %in% colnames(phenotypes)) {
      stop(paste0("a column with the same name '", factor.var,"' already exists in `phenotypes` data frame! Please change the `prefix.ordered.factor`!"))
    }
    
    # add the column to phenotypes:
    list.groups <- unique(phenotypes[[unordered.factor.var]])
    list.groups <- list.groups[list.groups != reference.group]# temporarily drop the reference group
    list.groups <- c(reference.group, list.groups)   # add as first
    phenotypes[[factor.var]] <- ordered(phenotypes[[unordered.factor.var]], levels = list.groups)   # the first element in the list.groups would be the reference group

  }
  
  # generate the formula:
  formula <- paste0(response.var, "~", factor.var, "+")
  formula <- paste0(formula, "s(", smooth.var, ",k=",toString(k),",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "s(", smooth.var, ",k=",toString(k),",by=", factor.var,",fx=", toString(fx), ")")
  
  # return
  formula <- as.formula(formula)  # when printing formula, it could be shown in more than one lines...
  toReturn = list(formula = formula,
                  phenotypes = phenotypes)
  return(toReturn)
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
  var.terms.full = c("estimate","std.error","statistic","p.value")
  var.model.full = c("r.squared", "adj.r.squared", "sigma", "statistic", "p.value", "df", "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs")
  if (full.outputs == TRUE) {   # full set of outputs
    var.terms = var.terms.full
    var.model = var.model.full
  }
  # check on validity of list of vars:
  var.terms <- var.terms[!duplicated(var.terms)]  # remove duplicated element(s)
  var.model <- var.model[!duplicated(var.model)]
  for (var in var.terms) {
    if (!(var %in% var.terms.full)) {
      stop(paste0(var, " is not valid for var.terms!"))
    }
  }
  for (var in var.model) {
    if (!(var %in% var.model.full)) {
      stop(paste0(var, " is not valid for var.model!"))
    }
  }


  # check for p.value correction:
    # check for terms:
  check_validity_correctPValue(correct.p.value.terms, "correct.p.value.terms",
                              var.terms, "var.terms")
    # check for model:
  check_validity_correctPValue(correct.p.value.model, "correct.p.value.model",
                              var.model, "var.model")


  
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
 
  
  


  df_out   # return

}



#' Run a GAM model at each fixel location
#' For p-value corrections (arguments correct.p.value.*), supported methods include all methods in `p.adjust.methods` except "none". Can be more than one method. Turn it off by setting to "none".
#' Please notice that there is no p.value for the model, so no "correct.p.value.model" for GAM model.
#' @param formula Formula (passed to `mgcv::gam()`)
#' @param data FixelArray class
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param fixel.subset The subset of fixel ids you want to run. Integers. First id starts from 1.
#' @param full.outputs Whether to return full set of outputs (TRUE or FALSE). If FALSE, it will only return those listed in var.terms and var.model; if TRUE, arguments var.terms and var.model will be ignored.
#' @param var.smoothTerms The list of variables to save for smooth terms (got from gam %>% tidy(parametric = FALSE)). Example smooth term: age in formula "outcome ~ s(age)".
#' @param var.parametricTerms The list of variables to save for parametric terms (got from gam %>% tidy(parametric = TRUE)). Example parametric term: sex in formula "outcome ~ s(age) + sex"
#' @param var.model The list of variables to save for the model (got from lm %>% glance())
#' @param eff.size.term.index The i-th term of the formula's right hand side as the term of interest for effect size. Positive integer or integer list. Usually term of interest is smooth term.
#' @param correct.p.value.smoothTerms To perform and add a column for p.value correction for each smooth term. 
#' @param correct.p.value.parametricTerms To perform and add a column for p.value correction for each parametric term. 
#' @param verbose Print progress messages
#' @param pbar Print progress bar
#' @param n_cores The number of cores to run on
#' @return Tibble with the summarized model statistics for all fixels requested
#' @import doParallel
#' @import tibble
#' @import mgcv
#' @export

FixelArray.gam <- function(formula, data, phenotypes, scalar, fixel.subset = NULL, full.outputs = FALSE, 
                              var.smoothTerms = c("statistic","p.value"),
                              var.parametricTerms = c("estimate", "statistic", "p.value"),
                              var.model = c("dev.expl"), 
                              eff.size.term.index = NULL,
                              correct.p.value.smoothTerms = "none", correct.p.value.parametricTerms = "none",
                              verbose = TRUE, pbar = TRUE, n_cores = 1, ...){
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
  
  ### TODO: print additional arguments in smooth term (s(), te(), etc) - only displaying the important arguments
  # check if the formula is valid in terms of mgcv::gam()
  tryCatch(
    {
      # try
      gam.formula.breakdown <- mgcv::interpret.gam(formula)  # if error, it means the formula is not valid in terms of mgcv::gam()
    },
    error = function(cond) {
      stop(paste0("The formula is not valid for mgcv::gam()! Please check and revise."))
    }
  )
  
  # to check formula, we need to fit one fixel:
  values <- scalars(data)[[scalar]][1,]
  dat <- phenotypes
  dat[[scalar]] <- values
  onemodel <- mgcv::gam(formula = formula, data = dat)
    
  #checker_gam_formula(formula, gam.formula.breakdown, onemodel)
    
  # what smooth? s or te or?
  # additional arguments in the smooth term, and are they valid for this specific term type?
  
  ### display additional arguments:
  dots <- list(...)    
  dots_names <- names(dots)
  
  FUN <- mgcv::gam
  
  # family:
  m <- invisible(eval(formals(FUN)$family))    # should not use message(), but print() --> but will print out or invisible()
  m1 <- paste0("Family: ", m$family, "; Link function: ", m$link)
  printAdditionalArgu(FUN, "family", dots, m1)
  
  
  # eval(formals(mgcv::gam)$data) # return the default setting of argument "data"
  
  
  # TODO: finish this part: display additional arguments
  
  # TODO: check if fx=FALSE; if so, add edf to the list of var + warning: fx=TRUE is recommended
  
  # when full.outputs = TRUE:
  var.smoothTerms.full <- c("edf","ref.df","statistic","p.value")
  var.parametricTerms.full <- c("estimate", "std.error","statistic","p.value")
  var.model.full <- c("adj.r.squared","dev.expl", "sp.criterion", "scale",
                      "df", "logLik","AIC", "BIC", "deviance", "df.residual", "nobs")

  if (full.outputs == TRUE) {   # full set of outputs
    var.smoothTerms <- var.smoothTerms.full
    var.parametricTerms <- var.parametricTerms.full
    var.model <- var.model.full
    
    # TODO: also add exiting smooth terms' index to eff.size.term.index
  }

  ### check on validity of arguments: var.term and var.model 
  var.smoothTerms <- var.smoothTerms[!duplicated(var.smoothTerms)]  # remove duplicated element(s)
  var.parametricTerms <- var.parametricTerms[!duplicated(var.parametricTerms)]
  var.model <- var.model[!duplicated(var.model)]
  
  # check if all var* are empty:
  if (length(var.smoothTerms)==0 & length(var.parametricTerms) ==0 & length(var.model) == 0) {
    stop("All var* [var.smoothTerms, var.parametricTerms, var.model] are empty!")
  }

  # check if every var is valid:
  for (var in var.smoothTerms) {   
    if (!(var %in% var.smoothTerms.full)) {
      stop(paste0(var, " is not valid for var.smoothTerms!"))
    }
  }
  for (var in var.parametricTerms) {
    if (!(var %in% var.parametricTerms.full)) {
      stop(paste0(var, " is not valid for var.parametricTerms!"))
    }
  }
  for (var in var.model) {
    if (!(var %in% var.model.full)) {
      stop(paste0(var, " is not valid!"))
    }
  }


  var.model.orig <- var.model
  if (!is.null(eff.size.term.index)) {    # eff.size is not null --> requested

    # check if the term index is valid:
    if (min(eff.size.term.index)<=0) {# any of not positive | can't really check if it's integer as is.integer(1) is FALSE...
      stop(paste0("There is element(s) in eff.size.term.index <= 0. It should be a (list of) positive integer!"))
    }

    terms.full.formula <- terms(formula, keep.order = TRUE)   # not the re-order the terms | see: https://rdrr.io/r/stats/terms.formula.html
    if (max(eff.size.term.index) > length(labels(terms.full.formula))) {  # if max is more than the number of terms on RHS of formula
      stop(paste0("Largest index in eff.size.term.index is more than the term number on the right hand side of formula!"))
    }
    
    # check how many variables on RHS; if no (but intercept, i.e. xx ~ 1), stop
    if (length(labels(terms.full.formula)) ==0) {
      stop(paste0("To analyze effect size but there is no variable (except intercept 1) on right hand side of formula! Please provide at least one valid variable."))
    }
    
    # print warning:
    message("will get effect size (eff.size) so the execution time will be longer.")
    # add adj.r.squared into var.model
    if (!("adj.r.squared" %in% var.model)) {
      var.model <- c(var.model, "adj.r.squared")
    }
  }

  ### check on arguments: p-values correction methods
  # check for smoothTerms:
  check_validity_correctPValue(correct.p.value.smoothTerms, "correct.p.value.smoothTerms",
                              var.smoothTerms, "var.smoothTerms")
  # check for parametricTerms:
  check_validity_correctPValue(correct.p.value.parametricTerms, "correct.p.value.parametricTerms",
                              var.parametricTerms, "var.parametricTerms")
  
  
  ### run
  # start the process:
  if(verbose){
    message(glue::glue("Fitting fixel-wise GAMs for {scalar}", ))
    message(glue::glue("initiating....", ))
  }

  # initiate: get the example of one fixel and get the column names
  outputs_initiator <- analyseOneFixel.gam(i_fixel=1, formula, data, phenotypes, scalar,
                                          var.smoothTerms, var.parametricTerms, var.model,
                                          flag_initiate = TRUE,
                                          ...)
  column_names <- outputs_initiator$column_names
  list.smoothTerms = outputs_initiator$list.smoothTerms
  list.parametricTerms = outputs_initiator$list.parametricTerms

  # loop (by condition of pbar and n_cores)
  if(verbose){
    message(glue::glue("looping across fixels....", ))
  }

  # is it a multicore process?
  flag_initiate <- FALSE
  if(n_cores > 1){
    
    if (pbar) {
      
      fits <- pbmcapply::pbmclapply(fixel.subset,   # a list of i_fixel
                                    analyseOneFixel.gam,  # the function
                                    mc.cores = n_cores,
                                    formula, data, phenotypes, scalar,
                                    var.smoothTerms, var.parametricTerms, var.model,
                                    flag_initiate = FALSE,
                                    ...)
      
    } else {
      
      # foreach::foreach
      
      fits <- parallel::mclapply(fixel.subset,   # a list of i_fixel 
                                 analyseOneFixel.gam,  # the function
                                 mc.cores = n_cores,
                                 formula, data, phenotypes, scalar,
                                 var.smoothTerms, var.parametricTerms, var.model,
                                 flag_initiate = FALSE,
                                 ...)
      
    }
  } else  {  # n_cores ==1, not multi-core
    
    if (pbar) {
      
      fits <- pbapply::pblapply(fixel.subset,   # a list of i_fixel
                                analyseOneFixel.gam,  # the function
                                formula, data, phenotypes, scalar,
                                var.smoothTerms, var.parametricTerms, var.model,
                                flag_initiate = FALSE,
                                ...)
      
    } else {
      
      fits <- lapply(fixel.subset,   # a list of i_fixel
                     analyseOneFixel.gam,  # the function
                     formula, data, phenotypes, scalar,
                     var.smoothTerms, var.parametricTerms, var.model,
                     flag_initiate = FALSE,
                     ...)
    }
  }  


  df_out <- do.call(rbind, fits)    
  df_out <- as.data.frame(df_out)    # turn into data.frame
  colnames(df_out) <- column_names     # add column names
  


  ### get the effect size for smooth terms:
  if (!is.null(eff.size.term.index)) {   # if eff.size is requested
    message("Getting the effect size: running the reduced model...")
  
    # list of term of interest for eff.size:
    eff.size.term.fullFormat.list <- labels(terms.full.formula)[eff.size.term.index]  # the term for effect size, in full format
    # get the short version:
    eff.size.term.shortFormat.list <- list()
    for (eff.size.term.fullFormat in eff.size.term.fullFormat.list) {
      
      temp <- strsplit(eff.size.term.fullFormat, "[(,)]")[[1]]   # format: s(age, k=xxx) --> s, age, k=xxxx
      eff.size.term.shortFormat <- paste0(temp[1], "_",temp[2])   # remove optional arguments in s(), replace () with _: get e.g. s_age
      eff.size.term.shortFormat.list <- append(eff.size.term.shortFormat.list, eff.size.term.shortFormat)
      # TODO: ask Bart: if it's appropriate to have s(age, k=3) + s(age, k=4). If not, throw out a warning saying it's not good. If not, it's okay to use shortFormat as column name of eff.size
    }
    
    # loop of each eff.size.term.index (i.e. each term of interest)
    for (i.eff.size.term in 1:length(eff.size.term.fullFormat.list)) {
      idx.eff.size.term <- eff.size.term.index[i.eff.size.term]   # index
      eff.size.term.fullFormat <- eff.size.term.fullFormat.list[i.eff.size.term]
      eff.size.term.shortFormat <- eff.size.term.shortFormat.list[[i.eff.size.term]][1]   # it's nested
      
      # get the formula of reduced model
        # check if there is only one term (after removing it in reduced model, there is no term but intercept in the formula...)
      if (length(labels(terms.full.formula)) ==1) {
        temp <- toString(formula) %>% strsplit("[, ]")  
        reduced.formula <- as.formula(paste0(temp[[1]][3], "~1"))
      } else {
        reduced.formula <- formula(drop.terms(terms.full.formula, idx.eff.size.term, keep.response = TRUE))  # index on RHS of formula -> change to class of formula (as.formula does not work)
      }
      
      message(paste0("* Getting effect size for term: ", eff.size.term.fullFormat, " via reduced model as below","; will show up as ", eff.size.term.shortFormat, " in final dataframe"))   # NOTES: it would be great to figure out how to paste and print formula without format changed. May try out stringf?
      print(reduced.formula)

      # var* for reduced model: only adjusted r sq is enough
      # initiate:
      reduced.model.outputs_initiator <- analyseOneFixel.gam(i_fixel=1, reduced.formula, data, phenotypes, scalar,
                                            var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                                            flag_initiate = TRUE,
                                            ...)
      reduced.model.column_names <- reduced.model.outputs_initiator$column_names
  
      # run on reduced model, get the adj r sq of reduced model
      if(n_cores > 1){
      
        if (pbar) {
          
          reduced.model.fits <- pbmcapply::pbmclapply(fixel.subset,   # a list of i_fixel
                                        analyseOneFixel.gam,  # the function
                                        mc.cores = n_cores,
                                        reduced.formula, data, phenotypes, scalar,
                                        var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                                        flag_initiate = FALSE,
                                        ...)
          
        } else {
          
          # foreach::foreach
          
          reduced.model.fits <- parallel::mclapply(fixel.subset,   # a list of i_fixel 
                                     analyseOneFixel.gam,  # the function
                                     mc.cores = n_cores,
                                     reduced.formula, data, phenotypes, scalar,
                                     var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                                     flag_initiate = FALSE,
                                     ...)
          
        }
      } else  {  # n_cores ==1, not multi-core
        
        if (pbar) {
          
          reduced.model.fits <- pbapply::pblapply(fixel.subset,   # a list of i_fixel
                                    analyseOneFixel.gam,  # the function
                                    reduced.formula, data, phenotypes, scalar,
                                    var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                                    flag_initiate = FALSE,
                                    ...)
          
        } else {
          
          reduced.model.fits <- lapply(fixel.subset,   # a list of i_fixel
                         analyseOneFixel.gam,  # the function
                         reduced.formula, data, phenotypes, scalar,
                         var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                         flag_initiate = FALSE,
                         ...)
        }
      }  # end of loop for calculating reduced model across fixels
  
      reduced.model.df_out <- do.call(rbind, reduced.model.fits)
      reduced.model.df_out <- as.data.frame(reduced.model.df_out)
      colnames(reduced.model.df_out) <- reduced.model.column_names
      
      # rename "adj.r.squared" as "redModel.adj.r.squared" before merging into df_out
      names(reduced.model.df_out)[names(reduced.model.df_out) == 'model.adj.r.squared'] <- "redModel.adj.r.squared"
        
      # combine new df_out to original one:
      df_out <- merge(df_out, reduced.model.df_out, by = "fixel_id")
      
      # calculate the eff.size, add to the df_out:
      df_out <- df_out %>% dplyr::mutate("{eff.size.term.shortFormat}.eff.size" := model.adj.r.squared - redModel.adj.r.squared)

      # remove column of redModel
      df_out <- df_out %>% subset(select = -c(redModel.adj.r.squared))
    }  # end of for loop across term of interest for effect size

    
  
    # if adjusted r sq is not requested (see var.model.orig), remove it:
    if (!("adj.r.squared" %in% var.model.orig)) {
      df_out <- df_out %>% subset(select = -c(model.adj.r.squared))
    }
    
  }   # end of if: requesting eff.size
  
  
  
  
  
  ### correct p values
  # add correction of p.values: for smoothTerms
  if ( all(correct.p.value.smoothTerms == "none") ) {    # all() is to accormodate for multiple elements in correct.p.value.smoothTerms: if one of is not "none", FALSE
    # do nothing
    
  } else {
    if ("p.value" %in% var.smoothTerms == TRUE) {   # check whether there is "p.value" in var.smoothTerms' | if FALSE: print warning (see beginning of this function)
      
      for (methodstr in correct.p.value.smoothTerms) {
        
        for (tempstr in list.smoothTerms) {
          tempstr.raw <- paste0(tempstr, ".p.value")
          tempstr.corrected <- paste0(tempstr.raw, ".", methodstr)
          temp.corrected <- p.adjust(df_out[[tempstr.raw]], method = methodstr)
          df_out <- df_out %>% tibble::add_column( "{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
        }
        
      }
      
    }
    
  }

  # add correction of p.values for parametricTerms
  if ( all(correct.p.value.parametricTerms == "none") ) {    # all() is to accormodate for multiple elements in correct.p.value.parametricTerms: if one of is not "none", FALSE
    # do nothing
    
  } else {
    if ("p.value" %in% var.parametricTerms == TRUE) {   # check whether there is "p.value" in var.parametricTerms' | if FALSE: print warning (see beginning of this function)
      
      for (methodstr in correct.p.value.parametricTerms) {
        
        for (tempstr in list.parametricTerms) {
          tempstr.raw <- paste0(tempstr, ".p.value")
          tempstr.corrected <- paste0(tempstr.raw, ".", methodstr)
          temp.corrected <- p.adjust(df_out[[tempstr.raw]], method = methodstr)
          df_out <- df_out %>% tibble::add_column( "{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
        }
        
      }
      
    }
    
  }



  ### return
  df_out
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