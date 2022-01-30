#' check if an object in .h5 exists
#' @param fn_h5 filename of the .h5 file
#' @param group_name full directory of this object in .h5 name
#' @param object_name name of the object, should be a string without "/"
#' @noRd
#' @importFrom rhdf5 h5ls h5closeAll
#' @importFrom dplyr filter
#' @importFrom rlang .data
flagObjectExistInh5 <- function(fn_h5, group_name="/results",object_name="myAnalysis") {
  
  rhdf5::h5closeAll()
  
  h5 <- rhdf5::h5ls(fn_h5)
  
  h5 %>% 
    dplyr::filter(.data$group==group_name & .data$name==object_name) %>% 
    nrow() -> h5.nrow
  
  if (h5.nrow==0) {
    object_exists <- FALSE
  } else {
    object_exists <- TRUE
  }
  
  object_exists
}



#' check if h5 group "results" exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @noRd
#' @importFrom rhdf5 h5ls h5closeAll
#' @importFrom dplyr filter
#' @importFrom rlang .data
flagResultsGroupExistInh5 <- function(fn_h5) {
  
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  
  h5 %>% 
    dplyr::filter(.data$group=="/" & .data$name=="results") %>% 
    nrow() -> h5.nrow
  
  if (h5.nrow==0) {
    object_exists <- FALSE
  } else {
    object_exists <- TRUE
  }
  
  object_exists
}

#' check if a subfolder of results exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @param analysis_name The subfolder name in "results" in .h5 file 
#' @noRd
#' @importFrom rhdf5 h5ls h5closeAll
#' @importFrom dplyr filter
#' @importFrom rlang .data
flagAnalysisExistInh5 <- function(fn_h5, analysis_name) {
  
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  
  h5 %>% 
    dplyr::filter(.data$group=="/results" & .data$name==analysis_name) %>% 
    nrow() -> h5.nrow
  
  if (h5.nrow==0) {
    object_exists <- FALSE
  } else {
    object_exists <- TRUE
  }
  
  object_exists
}


#' print the additional arguments settings
#' @param FUN The function, e.g. mgcv::gam, without "()"
#' @param argu_name The argument name of the function
#' @param dots: list of additional arguments
#' @param message_default The message for default 
#' @param message_usr_input The message describing user's input
#' @importFrom crayon black
#' @importFrom dplyr %>%
#' @noRd
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
#' @importFrom stats p.adjust.methods
#' @noRd
check_validity_correctPValue <- function(correct.list, name.correct.list, 
                                         var.list, name.var.list) {
  p.adjust.methods.full <- stats::p.adjust.methods[ stats::p.adjust.methods != "none" ]
  if ( all(correct.list == "none") == FALSE) {    # any element is not "none"
    checker.method.in <- correct.list %in% p.adjust.methods.full
    if ( all(checker.method.in) == FALSE) {   # not all "TRUE"
      stop(paste0("Some of elements in ",name.correct.list," are not valid. Valid inputs are: ", paste(p.adjust.methods.full, collapse = ', ')))
    } 
    
    if ("p.value" %in% var.list == FALSE) {  # not in the list | # check whether there is "p.value" in var.list
      warning(paste0("p.value was not included in ",name.var.list,", so not to perform its p.value corrections"))   # TODO: why this warning comes out after ModelArray.aModel is done?
    }
  }
}


#' Print out important arguments in smooth terms s() in mgcv::gam() formula
#' 
#' @details
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/s
#' 
#' @param ofInterest got via: gam.formula.breakdown <- mgcv::interpret.gam(formula); ofInterest <- gam.formula.breakdown$smooth.spec[[i]]
#' @importFrom mgcv s
#' @importFrom dplyr %>%
#' @importFrom crayon black
#' @noRd
#' 
checker_gam_s <- function(ofInterest) {
  FUN <- mgcv::s
  
  # add "by=?" back to term name:
  if (ofInterest$by == "NA") {
    term_name <- ofInterest$label
  } else {
    term_name <- paste0(substr(ofInterest$label, 1, nchar(ofInterest$label)-1), ", by=",ofInterest$by,")")
  }
  
  paste0(term_name, ": ") %>% crayon::black() %>% cat()
  
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
  
  cat("\n")
  
  
}

#' Print out important arguments in smooth term te() or ti() or t2() in mgcv::gam() formula
#' 
#' @details
#' Why a separate function is needed for t(), cannot using s(): in ofInterest, "fx" is "fx" for t(), but "fixed" for s() - so they are different.
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/te or /t2()
#' 
#' @param FUN could be mgcv::te(), ti() or t2()
#' @param ofInterest got via: gam.formula.breakdown <- mgcv::interpret.gam(formula); ofInterest <- gam.formula.breakdown$smooth.spec[[i]]
#' @importFrom mgcv te ti t2
#' @importFrom dplyr %>%
#' @importFrom crayon black
#' @noRd
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
  
  cat("\n")
}

#' A checker for formula in gam for ModelArray.gam()
#' @param formula The formula
#' @param gam.formula.breakdown Got from mgcv::interpret.gam(formula)
#' @param onemodel The model of one element got from mgcv::gam()
#' @importFrom mgcv s ti t2 te
#' @importFrom dplyr %>%
#' @importFrom crayon black
#' @noRd
#' 
checker_gam_formula <- function(formula, gam.formula.breakdown, onemodel=NULL) {
  # print out the formula:
  temp <- formula %>% as.character() 
  str_formula <- paste0(temp[2], " ~ ", temp[3])
  
  m1 <- paste0("The formula requested: ", str_formula) %>% crayon::black() %>% cat()
  cat(m1, "\n")
  
  
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
  
  
  # # fit for one element, get the summarized stat:
  # onemodel.tidy.smoothTerms <- onemodel %>% broom::tidy(parametric = FALSE)
  # onemodel.tidy.parametricTerms <- onemodel %>% broom::tidy(parametric = TRUE)
  # onemodel.glance <- onemodel %>% broom::glance()
  # onemodel.summary <- onemodel %>% summary()
  
  
}

#' Generate GAM formula with factor-smooth interaction
#' 
#' @description 
#' This function will generate a formula in the following format: \code{y ~ orderedFactor + s(x) + s(x, by=orderedFactor)},
#' where \code{y} is \code{response.var}, \code{x} is \code{smooth.var}, and \code{orderedFactor} is \code{factor.var} - see \code{factor.var} for more.
#' The formula generated could be further modified, e.g. adding covariates.
#' 
#' @param response.var character class, the variable name for response
#' @param factor.var character class, the variable name for factor. It should be an ordered factor. If not, it will generate it as a new column in `phenotypes`, which requires `reference.group`.
#' @param smooth.var character class, the variable name in smooth term as main effect
#' @param phenotypes data.frame class, the cohort matrix with columns of independent variables (including \code{factor.var} and \code{smooth.var}) to be added to the model 
#' @param reference.group character class, the reference group for ordered factor of `factor.var`; required when `factor.var` in `phenotypes` is not an ordered factor. 
#' @param prefix.ordered.factor character class, the prefix for ordered factor; required when `factor.var` in `phenotypes` is not an ordered factor.
#' @param fx TRUE or FALSE, to be used in smooth term s(). Recommend TRUE.
#' @param k integer, to be used in smooth term including the interaction term. If NULL (no entry), will use default value as in mgcv::s()
#' @return a list, including: 1) formula generated; 2) data.frame phenotypes - updated if argument factor.var is not an ordered factor
#' @importFrom mgcv s
#' @importFrom stats as.formula
#' @export
#' 
generator_gamFormula_factorXsmooth <- function(response.var, factor.var, smooth.var, phenotypes, 
                                               reference.group = NULL, prefix.ordered.factor = "o",
                                               fx=TRUE, k=NULL) {
  class.factor.var <- class(phenotypes[[factor.var]])
  if (  !( (length(class.factor.var) == 2) & (class.factor.var[1] == "ordered") & (class.factor.var[2] == "factor")  )  ) {   # class is not c("ordered", "factor")
    
    message("input `factor.var` is not an ordered factor; will generate one in data.frame `phenotypes` which will be returned")
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
  
  if (is.null(k)) {
    k <- invisible(eval(formals(mgcv::s)[["k"]]))
  }
  
  # generate the formula:
  formula <- paste0(response.var, "~", factor.var, "+")
  formula <- paste0(formula, "s(", smooth.var, ",k=",toString(k),",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "s(", smooth.var, ",k=",toString(k),",by=", factor.var,",fx=", toString(fx), ")")
  
  # return
  formula <- stats::as.formula(formula)  # when printing formula, it could be shown in more than one lines...
  toReturn = list(formula = formula,
                  phenotypes = phenotypes)
  return(toReturn)
}



#' Generate GAM formula with continuous*continuous interaction
#' 
#' @description 
#' This function will generate a formula in the following format: \code{y ~ ti(x) + ti(z) + ti(x,z)},
#' where \code{y} is \code{response.var}, \code{x} is \code{cont1.var}, and \code{z} is \code{cont2.var}.
#' The formula generated could be further modified, e.g. adding covariates.
#' 
#' @param response.var character class, the variable name for response
#' @param cont1.var character class, the name of the first continuous variable
#' @param cont2.var character class, the name of the second continuous variable
#' @param fx TRUE or FALSE, to be used in smooth term s(). Recommend TRUE.
#' @param k integer, to be used in smooth term including the interaction term. If NULL (no entry), will use default value as in mgcv::s()
#' @return The formula generated
#' @importFrom mgcv ti
#' @importFrom stats as.formula
#' @export
#' 
generator_gamFormula_continuousInteraction <- function(response.var, cont1.var, cont2.var,
                                                       fx=TRUE, k=NULL) {
  
  if (is.null(k)) {
    k <- invisible(eval(formals(mgcv::ti)[["k"]]))
  }
  
  formula <- paste0(response.var, "~")
  formula <- paste0(formula, "ti(", cont1.var, ",k=",toString(k),",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "ti(", cont2.var, ",k=",toString(k),",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "ti(", cont1.var, ",", cont2.var,",k=",toString(k),",fx=", toString(fx), ")")
  
  formula <- stats::as.formula(formula)
  return(formula)
}


#' @param fullmodel The full model object returned by e.g. mgcv::gam() with full formula
#' @param redmodel The reduced model object returned by e.g. mgcv::gam() with reduced formula (i.e. the term of interest has been removed)
#' @return A list including: partial R squared, sse for full and reduced models
#' 
#' @details When calculating reduced model `redmodel`, we recommend to use `fullmodel$model` (i.e. the data used in full model, after excluding subjects with NA) as the input data.frame, so that the list of subjects used is consistent between full and reduced model. 
#' This is mainly for a rare case, where e.g. first subject has missing age (NA), but s(age) is of interest for changed.rsq; then full model will not include it, and if still using the original data.frame, reduced model will include this subject...
partialRsq <- function(fullmodel, redmodel) {
  # calculating SSE: used observed y (i.e. excluding observations with NA), and fitted values, directly from model object
  
  sse.full <- sum( (fullmodel$y - fullmodel$fitted.values)^2 )
  sse.red <- sum( (redmodel$y - redmodel$fitted.values)^2 )
  
  partialRsq <- (sse.red - sse.full) / sse.red
  
  toReturn <- list(partialRsq = partialRsq,
                   sse.full = sse.full,
                   sse.red = sse.red)
  return(toReturn)
}

