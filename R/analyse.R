#' print the additional arguments settings
#' @param FUN The function, e.g. mgcv::gam, without "()"
#' @param argu_name The argument name of the function
#' @param dots: list of additional arguments
#' @param message_default The message for default 
#' @param message_usr_input The message describing user's input
#' @importFrom crayon black
#' @importFrom dplyr %>%
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
#' @importFrom stats p.adjust.methods

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
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/s
#' TODO: finish the description
#' @param ofInterest got via: gam.formula.breakdown <- mgcv::interpret.gam(formula); ofInterest <- gam.formula.breakdown$smooth.spec[[i]]
#' @importFrom mgcv s
#' @importFrom dplyr %>%
#' @importFrom crayon black
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
#' Why a separate function is needed for t(), cannot using s(): in ofInterest, "fx" is "fx" for t(), but "fixed" for s() - so they are different.
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/te or /t2()
#' TODO: finish the description
#' @param FUN could be mgcv::te(), ti() or t2()
#' @param ofInterest got via: gam.formula.breakdown <- mgcv::interpret.gam(formula); ofInterest <- gam.formula.breakdown$smooth.spec[[i]]
#' @importFrom mgcv te ti t2
#' @importFrom dplyr %>%
#' @importFrom crayon black
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
#' TODO: finish the description
#' @importFrom mgcv s ti t2 te
#' @importFrom dplyr %>%
#' @importFrom crayon black
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


#' Fit linear model for element-wise data
#' 
#' @description 
#' `ModelArray.lm` fits linear model (`stats::lm()`) for each of elements requested, and returns a tibble dataframe of requested model statistics.
#' 
#' @details 
#' You may request returning specific statistical variables by setting \code{var.*}, or you can get all by setting \code{full.outputs=TRUE}. 
#' Note that statistics covered by \code{full.outputs} or \code{var.*} are the ones from broom::tidy() and broom::glance() only, and do not include corrected p values.
#' List of acceptable statistic names for each of \code{var.*}:
#' \itemize{
#'   \item \code{var.terms}: c("estimate","std.error","statistic","p.value"); For interpretation please see `broom:tidy()`.
#'   \item \code{var.model}: c("r.squared", "adj.r.squared", "sigma", "statistic", "p.value", "df", "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs"); For interpretation please see `broom::glance()`.
#' }
#' For p-value corrections (arguments \code{correct.p.value.*}), supported methods include all methods in `p.adjust.methods` except "none". Can be more than one method. Turn it off by setting to "none".
#'
#' @param formula Formula (passed to `stats::lm()`)
#' @param data ModelArray class
#' @param phenotypes A data.frame of the cohort with columns of independent variables and covariates to be added to the model. It should contains a column called "source_file", and this column should match to that in \code{data}.
#' @param scalar A character. The name of the element-wise scalar to be analysed
#' @param element.subset A list of positive integers (min = 1, max = number of elements). The subset of elements you want to run. Default is `NULL`, i.e. requesting all elements in `data`.
#' @param full.outputs TRUE or FALSE, Whether to return full set of outputs. If FALSE, it will only return those listed in arguments \code{var.*}; if TRUE, arguments \code{var.*} will be ignored.
#' @param var.terms A list of characters. The list of variables to save for terms (got from `broom::tidy()`). See "Details" section for more.
#' @param var.model A list of characters. The list of variables to save for the model (got from `broom::glance()`). See "Details" section for more.
#' @param correct.p.value.terms A list of characters. To perform and add a column for p.value correction for each term. See "Details" section for more.
#' @param correct.p.value.model A list of characters. To perform and add a column for p.value correction for the model. See "Details" section for more.
#' @param verbose TRUE or FALSE, to print verbose message or not
#' @param pbar TRUE or FALSE, to print progress bar or not
#' @param n_cores Positive integer, The number of CPU cores to run with
#' @param ... Additional arguments for `stats::lm()`
#' @return Tibble with the summarized model statistics for all elements requested
#' @importFrom dplyr %>%
#' @import doParallel
#' @import tibble
#' @importFrom stats p.adjust lm
#' @importFrom glue glue
#' @export

ModelArray.lm <- function(formula, data, phenotypes, scalar, element.subset = NULL, full.outputs = FALSE, 
                              var.terms = c("estimate", "statistic", "p.value"), 
                              var.model = c("adj.r.squared", "p.value"), 
                          correct.p.value.terms = "none", correct.p.value.model = "none",
                              verbose = TRUE, pbar = TRUE, n_cores = 1, ...) {
  # data type assertions
  if(class(data) != "ModelArray") {
    stop("data's class is not ModelArray!")
  }
  
  ## element.subset:
  if (is.null(element.subset)) {  # request all elements
    num.element.total <- numElementsTotal(modelarray=data, scalar_name = scalar)
    element.subset <- 1:num.element.total
  }
  # checker for min and max of element.subset; and whether elements are integer
  if (min(element.subset) < 1) {
    stop("Minimal value in element.subset should >= 1")
  }
  if (max(element.subset) > nrow(scalars(data)[[scalar]])) {
    stop(paste0("Maximal value in element.subset should <= number of elements = "), as.character(nrow(scalars(data)[[scalar]])))
  }
  if (class(element.subset) != "integer") {
    stop("Please enter integers for element.subset!")
  }
  
  ### sanity check: whether they match: modelarray's source file list and phenotypes' source file list:
  sources.modelarray <- sources(data)[[scalar]]
  sources.phenotypes <- phenotypes[["source_file"]]
  if (is.null(sources.phenotypes)) {
    stop(paste0("Did not find column 'source_file' in argument 'phenotypes'. Please check!"))
  }
  
  ## length should be the same:
  if (length(sources.modelarray) != length(sources.phenotypes)) {
    stop(paste0("The length of source file list from phenotypes's column 'source_file' is not the same as that in ModelArray 'data'! Please check out! The latter one can be accessed by: sources(data)[[scalar]]"))
  }
  
  ## check if the list is unique:
  if (length(sources.modelarray) != length(unique(sources.modelarray))) {
    stop(paste0("The source files in ModelArray 'data' are not unique! Please check out! It can be accessed by: sources(data)[[scalar]]"))
  }
  if (length(sources.phenotypes) != length(unique(sources.phenotypes)) ) {
    stop(paste0("The source files from phenotypes's column 'source_file' are not unique! Please check out and remove the duplicated one!"))
  }
  
  if (identical(sources.modelarray, sources.phenotypes)) {
    # identical, pass
  } else {  # not identical (but length is the same):
    # check if two lists can be matched (i.e. no unmatched source filename)
    if ((all(sources.modelarray %in% sources.phenotypes)) & ((all(sources.phenotypes %in% sources.modelarray)))) {
      # can be matched, just the order is different. Use match() function:
      reorder_idx <- match(sources.modelarray,  # vector of values in the order we want
                           sources.phenotypes)   # vector to be reordered
      # apply to phenotypes:
      phenotypes <- phenotypes[reorder_idx, ]
      row.names(phenotypes) <- NULL # reset the row name, just to be safe for later adding scalar values... see ModelArray_paper/notebooks/test_match_sourceFiles.Rmd
      if (!identical(phenotypes[["source_file"]], sources.modelarray)) {
        stop("matching source file names were not successful...")
      }
    } else {
      stop(paste0("phenotypes's column 'source_file' have different element(s) from source file list in ModelArray 'data'! Please check out! The latter one can be accessed by: sources(data)[[scalar]]"))
    }
    
    # stop(paste0("The source file list from phenotypes's column 'source_file' is not identical to that in ModelArray 'data'! Please check out! The latter one can be accessed by: sources(data)[[scalar]] "))
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
    message(glue::glue("Fitting element-wise linear models for {scalar}", ))
    message(glue::glue("initiating....", ))
  }
  

  
  # initiate: get the example of one element and get the column names
  outputs_initiator <- analyseOneElement.lm(i_element=1, formula, data, phenotypes, scalar, 
                                     var.terms, var.model, 
                                     flag_initiate = TRUE, 
                                     ...)
  column_names <- outputs_initiator$column_names
  list.terms <- outputs_initiator$list.terms

  # loop (by condition of pbar and n_cores)
  if(verbose){
    message(glue::glue("looping across elements....", ))
  }
  
  # is it a multicore process?
  flag_initiate <- FALSE
  if(n_cores > 1){
    
    if (pbar) {
      
      fits <- pbmcapply::pbmclapply(element.subset,   # a list of i_element
                                    analyseOneElement.lm,  # the function
                                    mc.cores = n_cores,
                                    formula, data, phenotypes, scalar,
                                    var.terms, var.model,
                                    flag_initiate = FALSE,
                                    ...)
      
    } else {
      
      # foreach::foreach
      
      fits <- parallel::mclapply(element.subset,   # a list of i_element 
                                 analyseOneElement.lm,  # the function
                                 mc.cores = n_cores,
                                 formula, data, phenotypes, scalar,
                                 var.terms, var.model,
                                 flag_initiate = FALSE,
                                 ...)
      
    }
  } else  {  # n_cores ==1, not multi-core
    
    if (pbar) {
      
      fits <- pbapply::pblapply(element.subset,   # a list of i_element
                                analyseOneElement.lm,  # the function
                                formula, data, phenotypes, scalar,
                                var.terms, var.model,
                                flag_initiate = FALSE,
                                ...)
      
    } else {
      
      fits <- lapply(element.subset,   # a list of i_element
                     analyseOneElement.lm,  # the function
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
          temp.corrected <- stats::p.adjust(df_out[[tempstr.raw]], method = methodstr)
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
        temp.corrected <- stats::p.adjust(df_out[[tempstr.raw]], method = methodstr)
        df_out <- df_out %>% tibble::add_column( "{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
        
      }
        
    }
      
  }
 
  
  


  df_out   # return

}



#' Run GAM for element-wise data
#' 
#' @description 
#' `ModelArray.gam` fits gam model for each of elements requested, and returns a tibble dataframe of requested model statistics.
#' 
#' @details 
#' You may request returning specific statistical variables by setting \code{var.*}, or you can get all by setting \code{full.outputs=TRUE}. 
#' Note that statistics covered by \code{full.outputs} or \code{var.*} are the ones from broom::tidy(), broom::glance(), and summary() only, and do not include effect size or corrected p values.
#' List of acceptable statistic names for each of \code{var.*}:
#' \itemize{
#'  \item \code{var.smoothTerms}: c("edf","ref.df","statistic","p.value"); For interpretation please see `broom::tidy(parametric=FALSE)`.
#'  \item \code{var.parametricTerms}: c("estimate", "std.error","statistic","p.value"); For interpretation please see `broom::tidy(parametric=TRUE)`.
#'  \item \code{var.model}: c("adj.r.squared","dev.expl", "sp.criterion", "scale", "df", "logLik","AIC", "BIC", "deviance", "df.residual", "nobs"); "adj.r.squared" is \code{r.sq} from `summary()`; "sp.criterion" is \code{sp.criterion} from `summary()`; For interpretation please see `broom::glance()` and `summary()`.
#' }
#' Regarding formula: So far these kinds of formula are tested:
#' \itemize{
#'   \item formula with smooth term, but without any interactions. Examples like \code{y ~ s(x) + orderedFactor}; \code{y ~ s(x) + s(z)}
#'   \item formula with interaction, but limited to only one interaction term, and in the formats of:
#'   \itemize{
#'       \item Formula #1: \code{y ~ orderedFactor + s(x) + s(x, by=orderedFactor) + other_covariate}, where \code{orderedFactor} should be discrete variables and generated by `ordered`. The interaction term will be displayed as "s_x_BYorderedFactor" in the column name in returned data.frame. You may use function `generator_gamFormula_factorXsmooth()` to generate one.
#'       \item Formula #2: \code{y ~ ti(x) + ti(z) + ti(x,z) + other_covariate}, where \code{x} and \code{z} should be continuous variables. The interaction term will be displayed as "ti_x_z" in the column name in the returned data.frame. You may use function `generator_gamFormula_continuousInteraction()` to generate one.
#'   }
#' }
#' Effect size is calculated by the difference between adjusted R squared of full model (formula requested) and that of reduced model (formula without the term requested)
#' \itemize{
#'   \item When requesting effect size, \code{fx} should be set as \code{TRUE}, so that degree of freedom is fixed.
#'   \item For formula with interactions, only formula in above formats are tested, and only effect size for interaction term is validated. The effect size for main effect (such as s(x) in Formula #1) may not "functionally" be its effect size, as the definition should be changed to reduced formula without both main effect and interaction term.
#' }
#' For p-value corrections (arguments \code{correct.p.value.*}), supported methods include all methods in `p.adjust.methods` except "none". Can be more than one method. Turn it off by setting to "none".
#' Please notice that different from `ModelArray.lm`, there is no p.value for the GAM model, so no "correct.p.value.model" for GAM model.
#' @param formula Formula (passed to `mgcv::gam()`)
#' @param data ModelArray class
#' @param phenotypes A data.frame of the cohort with columns of independent variables and covariates to be added to the model. It should contains a column called "source_file", and this column should match to that in \code{data}.
#' @param scalar A character. The name of the element-wise scalar to be analysed
#' @param element.subset A list of positive integers (min = 1, max = number of elements). The subset of elements you want to run. Default is `NULL`, i.e. requesting all elements in `data`.
#' @param full.outputs TRUE or FALSE, Whether to return full set of outputs. If FALSE, it will only return those listed in arguments \code{var.*}; if TRUE, arguments \code{var.*} will be ignored.
#' @param var.smoothTerms A list of characters. The list of variables to save for smooth terms (got from `broom::tidy(parametric = FALSE)`). Example smooth term: age in formula "outcome ~ s(age)". See "Details" section for more.
#' @param var.parametricTerms A list of characters. The list of variables to save for parametric terms (got from `broom::tidy(parametric = TRUE)`). Example parametric term: sex in formula "outcome ~ s(age) + sex". See "Details" section for more.
#' @param var.model A list of characters. The list of variables to save for the model (got from `broom::glance()` and `summary()`). See "Details" section for more.
#' @param eff.size.term.index A list of (one or several) positive integers. Each element in the list means the i-th term of the formula's right hand side as the term of interest for effect size. Effect size will be calculated for each of term requested. Positive integer or integer list. Usually term of interest is smooth term, or interaction term in models with interactions.
#' @param correct.p.value.smoothTerms A list of characters. To perform and add a column for p.value correction for each smooth term. See "Details" section for more.
#' @param correct.p.value.parametricTerms A list of characters. To perform and add a column for p.value correction for each parametric term. See "Details" section for more.
#' @param verbose TRUE or FALSE, to print verbose messages or not
#' @param pbar TRUE or FALSE, to print progress bar or not
#' @param n_cores Positive integer, The number of CPU cores to run with
#' @param ... Additional arguments for `mgcv::gam()`
#' @return Tibble with the summarized model statistics for all elements requested
#' @importFrom dplyr %>% mutate
#' @import doParallel
#' @import tibble
#' @import mgcv
#' @importFrom stats terms as.formula drop.terms p.adjust
#' @importFrom glue glue
#' @export

ModelArray.gam <- function(formula, data, phenotypes, scalar, element.subset = NULL, full.outputs = FALSE, 
                              var.smoothTerms = c("statistic","p.value"),
                              var.parametricTerms = c("estimate", "statistic", "p.value"),
                              var.model = c("dev.expl"), 
                              eff.size.term.index = NULL,
                              correct.p.value.smoothTerms = "none", correct.p.value.parametricTerms = "none",
                              verbose = TRUE, pbar = TRUE, n_cores = 1, ...){
  # data type assertions
  if(class(data) != "ModelArray") {
    stop("data's class is not ModelArray!")
  }
  
  ## element.subset:
  if (is.null(element.subset)) {  # request all elements
    num.element.total <- numElementsTotal(modelarray=data, scalar_name = scalar)
    element.subset <- 1:num.element.total
  }
  # checker for min and max of element.subset; and whether elements are integer
  if (min(element.subset) < 1) {
    stop("Minimal value in element.subset should >= 1")
  }
  if (max(element.subset) > nrow(scalars(data)[[scalar]])) {
    stop(paste0("Maximal value in element.subset should <= number of elements = "), as.character(nrow(scalars(data)[[scalar]])))
  }
  if (class(element.subset) != "integer") {
    stop("Please enter integers for element.subset!")
  }
  
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
  
  # print out the additional arguments in smooth terms:
  
  # # to check formula, we need to fit one element:
  # values <- scalars(data)[[scalar]][1,]
  # dat <- phenotypes
  # dat[[scalar]] <- values
  # onemodel <- mgcv::gam(formula = formula, data = dat)
    
  # checker_gam_formula(formula, gam.formula.breakdown, onemodel)

  checker_gam_formula(formula, gam.formula.breakdown)
    

  ### sanity check: whether they match: modelarray's source file list and phenotypes' source file list:
  sources.modelarray <- sources(data)[[scalar]]
  sources.phenotypes <- phenotypes[["source_file"]]
  if (is.null(sources.phenotypes)) {
    stop(paste0("Did not find column 'source_file' in argument 'phenotypes'. Please check!"))
  }
  
  ## length should be the same:
  if (length(sources.modelarray) != length(sources.phenotypes)) {
    stop(paste0("The length of source file list from phenotypes's column 'source_file' is not the same as that in ModelArray 'data'! Please check out! The latter one can be accessed by: sources(data)[[scalar]]"))
  }
  
  ## check if the list is unique:
  if (length(sources.modelarray) != length(unique(sources.modelarray))) {
    stop(paste0("The source files in ModelArray 'data' are not unique! Please check out! It can be accessed by: sources(data)[[scalar]]"))
  }
  if (length(sources.phenotypes) != length(unique(sources.phenotypes)) ) {
    stop(paste0("The source files from phenotypes's column 'source_file' are not unique! Please check out and remove the duplicated one!"))
  }
  
  if (identical(sources.modelarray, sources.phenotypes)) {
    # identical, pass
  } else {  # not identical (but length is the same):
    # check if two lists can be matched (i.e. no unmatched source filename)
    if ((all(sources.modelarray %in% sources.phenotypes)) & ((all(sources.phenotypes %in% sources.modelarray)))) {
      # can be matched, just the order is different. Use match() function:
      reorder_idx <- match(sources.modelarray,  # vector of values in the order we want
                           sources.phenotypes)   # vector to be reordered
      # apply to phenotypes:
      phenotypes <- phenotypes[reorder_idx, ]
      row.names(phenotypes) <- NULL # reset the row name, just to be safe for later adding scalar values... see ModelArray_paper/notebooks/test_match_sourceFiles.Rmd
      if (!identical(phenotypes[["source_file"]], sources.modelarray)) {
        stop("matching source file names were not successful...")
      }
    } else {
      stop(paste0("phenotypes's column 'source_file' have different element(s) from source file list in ModelArray 'data'! Please check out! The latter one can be accessed by: sources(data)[[scalar]]"))
    }
    
    # stop(paste0("The source file list from phenotypes's column 'source_file' is not identical to that in ModelArray 'data'! Please check out! The latter one can be accessed by: sources(data)[[scalar]] "))
  }



  ### display additional arguments: [only important one]
  dots <- list(...)    
  dots_names <- names(dots)
  
  FUN <- mgcv::gam
  
  # # family:  # it works; but family may not be important
  # m <- invisible(eval(formals(FUN)$family))    # should not use message(), but print() --> but will print out or invisible()
  # m1 <- paste0("Family: ", m$family, "; Link function: ", m$link)
  # printAdditionalArgu(FUN, "family", dots, m1)
  
  # method: (default: "GCV.Cp")
  printAdditionalArgu(FUN, "method", dots)  # default: "GCV.Cp"

  # TODO: optional: check if fx=FALSE; if so, add edf to the list of var + warning: fx=TRUE is recommended
  
  
  
  # when full.outputs = TRUE:
  var.smoothTerms.full <- c("edf","ref.df","statistic","p.value")
  var.parametricTerms.full <- c("estimate", "std.error","statistic","p.value")
  var.model.full <- c("adj.r.squared","dev.expl", "sp.criterion", "scale",
                      "df", "logLik","AIC", "BIC", "deviance", "df.residual", "nobs")

  if (full.outputs == TRUE) {   # full set of outputs
    var.smoothTerms <- var.smoothTerms.full
    var.parametricTerms <- var.parametricTerms.full
    var.model <- var.model.full
    
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

    terms.full.formula <- stats::terms(formula, keep.order = TRUE)   # not the re-order the terms | see: https://rdrr.io/r/stats/terms.formula.html
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
    message(glue::glue("Fitting element-wise GAMs for {scalar}", ))
    message(glue::glue("initiating....", ))
  }

  # initiate: get the example of one element and get the column names
  outputs_initiator <- analyseOneElement.gam(i_element=1, formula, data, phenotypes, scalar,
                                          var.smoothTerms, var.parametricTerms, var.model,
                                          flag_initiate = TRUE,
                                          ...)
  column_names <- outputs_initiator$column_names
  list.smoothTerms = outputs_initiator$list.smoothTerms
  list.parametricTerms = outputs_initiator$list.parametricTerms

  # loop (by condition of pbar and n_cores)
  if(verbose){
    message(glue::glue("looping across elements....", ))
  }

  # is it a multicore process?
  flag_initiate <- FALSE
  if(n_cores > 1){
    
    if (pbar) {
      
      fits <- pbmcapply::pbmclapply(element.subset,   # a list of i_element
                                    analyseOneElement.gam,  # the function
                                    mc.cores = n_cores,
                                    formula, data, phenotypes, scalar,
                                    var.smoothTerms, var.parametricTerms, var.model,
                                    flag_initiate = FALSE,
                                    ...)
      
    } else {
      
      # foreach::foreach
      
      fits <- parallel::mclapply(element.subset,   # a list of i_element 
                                 analyseOneElement.gam,  # the function
                                 mc.cores = n_cores,
                                 formula, data, phenotypes, scalar,
                                 var.smoothTerms, var.parametricTerms, var.model,
                                 flag_initiate = FALSE,
                                 ...)
      
    }
  } else  {  # n_cores ==1, not multi-core
    
    if (pbar) {
      
      fits <- pbapply::pblapply(element.subset,   # a list of i_element
                                analyseOneElement.gam,  # the function
                                formula, data, phenotypes, scalar,
                                var.smoothTerms, var.parametricTerms, var.model,
                                flag_initiate = FALSE,
                                ...)
      
    } else {
      
      fits <- lapply(element.subset,   # a list of i_element
                     analyseOneElement.gam,  # the function
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
      temp <- strsplit(eff.size.term.fullFormat, "[(]")[[1]]
      if (length(temp) == 1) {  # it's not a smooth term - as there is no ()
        str_valid <- eff.size.term.fullFormat
      } else {
        smooth.class <- temp[1]
        
        theEval <- eval(parse(text = eff.size.term.fullFormat))
        str_valid <- paste0(smooth.class, "_",
                            paste(theEval$term, collapse = "_"))  # ti(x,z) --> ti_x_z; s(x) --> s_x
        if (theEval$by != "NA") {
          str_valid <- paste0(str_valid, "_BY", theEval$by)   # s(age,by=oSex) --> s_age_BYoSex
        }
      }
      
      eff.size.term.shortFormat <- str_valid  
      eff.size.term.shortFormat.list <- append(eff.size.term.shortFormat.list, eff.size.term.shortFormat)
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
        reduced.formula <- stats::as.formula(paste0(temp[[1]][3], "~1"))
      } else {
        reduced.formula <- formula(stats::drop.terms(terms.full.formula, idx.eff.size.term, keep.response = TRUE))  # index on RHS of formula -> change to class of formula (stats::as.formula does not work)
      }
      
      message(paste0("* Getting effect size for term: ", eff.size.term.fullFormat, " via reduced model as below","; will show up as ", eff.size.term.shortFormat, " in final dataframe"))   # NOTES: it would be great to figure out how to paste and print formula without format changed. May try out stringf?
      print(reduced.formula)

      # var* for reduced model: only adjusted r sq is enough
      # initiate:
      reduced.model.outputs_initiator <- analyseOneElement.gam(i_element=1, reduced.formula, data, phenotypes, scalar,
                                            var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                                            flag_initiate = TRUE,
                                            ...)
      reduced.model.column_names <- reduced.model.outputs_initiator$column_names
  
      # run on reduced model, get the adj r sq of reduced model
      if(n_cores > 1){
      
        if (pbar) {
          
          reduced.model.fits <- pbmcapply::pbmclapply(element.subset,   # a list of i_element
                                        analyseOneElement.gam,  # the function
                                        mc.cores = n_cores,
                                        reduced.formula, data, phenotypes, scalar,
                                        var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                                        flag_initiate = FALSE,
                                        ...)
          
        } else {
          
          # foreach::foreach
          
          reduced.model.fits <- parallel::mclapply(element.subset,   # a list of i_element 
                                     analyseOneElement.gam,  # the function
                                     mc.cores = n_cores,
                                     reduced.formula, data, phenotypes, scalar,
                                     var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                                     flag_initiate = FALSE,
                                     ...)
          
        }
      } else  {  # n_cores ==1, not multi-core
        
        if (pbar) {
          
          reduced.model.fits <- pbapply::pblapply(element.subset,   # a list of i_element
                                    analyseOneElement.gam,  # the function
                                    reduced.formula, data, phenotypes, scalar,
                                    var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                                    flag_initiate = FALSE,
                                    ...)
          
        } else {
          
          reduced.model.fits <- lapply(element.subset,   # a list of i_element
                         analyseOneElement.gam,  # the function
                         reduced.formula, data, phenotypes, scalar,
                         var.smoothTerms=c(), var.parametricTerms=c(), var.model=c("adj.r.squared"),
                         flag_initiate = FALSE,
                         ...)
        }
      }  # end of loop for calculating reduced model across elements
  
      reduced.model.df_out <- do.call(rbind, reduced.model.fits)
      reduced.model.df_out <- as.data.frame(reduced.model.df_out)
      colnames(reduced.model.df_out) <- reduced.model.column_names
      
      # rename "adj.r.squared" as "redModel.adj.r.squared" before merging into df_out
      names(reduced.model.df_out)[names(reduced.model.df_out) == 'model.adj.r.squared'] <- "redModel.adj.r.squared"
        
      # combine new df_out to original one:
      df_out <- merge(df_out, reduced.model.df_out, by = "element_id")
      
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
          temp.corrected <- stats::p.adjust(df_out[[tempstr.raw]], method = methodstr)
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
          temp.corrected <- stats::p.adjust(df_out[[tempstr.raw]], method = methodstr)
          df_out <- df_out %>% tibble::add_column( "{tempstr.corrected}" := temp.corrected, .after = tempstr.raw)
        }
        
      }
      
    }
    
  }



  ### return
  df_out
}




## TODO: add an example function for developer: ModelArray.gam()

