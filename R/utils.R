#' check if an object in .h5 exists
#' @param fn_h5 filename of the .h5 file
#' @param group_name full directory of this object in .h5 name
#' @param object_name name of the object, should be a string without "/"
#' @noRd
flagObjectExistInh5 <- function(fn_h5, group_name = "/results", object_name = "myAnalysis") {
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  h5.nrow <- nrow(h5[h5$group == group_name & h5$name == object_name, ])
  h5.nrow > 0
}


#' check if h5 group "results" exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @noRd
flagResultsGroupExistInh5 <- function(fn_h5) {
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  h5.nrow <- nrow(h5[h5$group == "/" & h5$name == "results", ])
  h5.nrow > 0
}

#' check if a subfolder of results exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @param analysis_name The subfolder name in "results" in .h5 file
#' @noRd
flagAnalysisExistInh5 <- function(fn_h5, analysis_name) {
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  h5.nrow <- nrow(h5[h5$group == "/results" & h5$name == analysis_name, ])
  h5.nrow > 0
}


#' print the additional arguments settings
#' @param FUN The function, e.g. mgcv::gam, without "()"
#' @param argu_name The argument name of the function
#' @param dots list of additional arguments
#' @param message_default The message for default
#' @param message_usr_input The message describing user's input
#' @importFrom dplyr %>%
#' @noRd
printAdditionalArgu <- function(FUN, argu_name, dots, message_default = NULL, message_usr_input = NULL) {
  dots_names <- names(dots)
  if (argu_name %in% dots_names) {
    if (is.null(message_default)) {
      message_default <- invisible(eval(formals(FUN)[[argu_name]]))
    }
    if (is.null(message_default)) { # if the default is NULL:
      message_default <- "NULL"
    }

    if (is.null(message_usr_input)) {
      m1 <- paste0(argu_name, " = ", dots[[argu_name]], " (default: ", message_default, ")") %>%
        crayon::black() %>%
        cat() # or, %>% message()
    } else { # specified the message:
      m1 <- paste0(argu_name, " = ", message_usr_input, " (default: ", message_default, ")") %>%
        crayon::black() %>%
        cat()
    }
  } else {
    m1 <- paste0(argu_name, ": default") %>%
      crayon::black() %>%
      cat()
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
#' @noRd
check_validity_correctPValue <- function(correct.list, name.correct.list,
                                         var.list, name.var.list) {
  p.adjust.methods.full <- stats::p.adjust.methods[stats::p.adjust.methods != "none"]
  if (all(correct.list == "none") == FALSE) { # any element is not "none"
    checker.method.in <- correct.list %in% p.adjust.methods.full
    if (all(checker.method.in) == FALSE) { # not all "TRUE"
      stop(
        paste0(
          "Some of elements in ",
          name.correct.list,
          " are not valid. Valid inputs are: ",
          paste(p.adjust.methods.full, collapse = ", ")
        )
      )
    }

    if ("p.value" %in% var.list == FALSE) { # not in the list | # check whether there is "p.value" in var.list
      warning(
        paste0(
          "p.value was not included in ",
          name.var.list,
          ", so not to perform its p.value corrections"
        )
      ) # TODO: why this warning comes out after ModelArray.aModel is done?
    }
  }
}


#' Print out important arguments in smooth terms s() in mgcv::gam() formula
#'
#' @details
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/s
#'
#' @param ofInterest got via: `gam.formula.breakdown <- mgcv::interpret.gam(formula)`;
#' `ofInterest <- gam.formula.breakdown$smooth.spec[[i]]`
#' @importFrom dplyr %>%
#' @noRd
#'
checker_gam_s <- function(ofInterest) {
  FUN <- mgcv::s

  # add "by=?" back to term name:
  if (ofInterest$by == "NA") {
    term_name <- ofInterest$label
  } else {
    term_name <- paste0(substr(ofInterest$label, 1, nchar(ofInterest$label) - 1), ", by=", ofInterest$by, ")")
  }

  paste0(term_name, ": ") %>%
    crayon::black() %>%
    cat()

  ### k (or bs.dim):   # could be multiple values
  m1 <- invisible(eval(formals(FUN)[["k"]])) # default
  m2 <- ofInterest$bs.dim # could be a list of multiple values
  if ((length(unique(m2)) == 1) && (m1 %in% unique(m2))) { # default
    msg_k <- " (default)"
  } else {
    msg_k <- ""
  }

  paste0(
    "  k = ",
    paste(as.character(m2), collapse = ", "), msg_k, "; "
  ) %>%
    crayon::black() %>%
    cat()

  ### fx:
  m1 <- invisible(eval(formals(FUN)[["fx"]])) %>% as.character() # default
  m2 <- ofInterest$fixed # actual

  if (as.character(!as.logical(m1)) %in% m2) { # there is an opposite logical value in m2 (actual)
    msg_fx <- ""
  } else {
    msg_fx <- " (default)"
  }

  paste0("  fx = ", toString(ofInterest$fixed), msg_fx, "; ") %>%
    crayon::black() %>%
    cat()

  ### bs:
  m1 <- invisible(eval(formals(FUN)[["bs"]])) # default
  # actual:
  mybs <- gsub(".smooth.spec", "", class(ofInterest)) # if there are multiple elements in bs, mybs will be a list
  if ((length(unique(mybs)) == 1) && (m1 %in% unique(mybs))) { # default
    msg_bs <- " (default)"
  } else {
    msg_bs <- ""
  }

  paste0(
    "  bs = ",
    paste(as.character(mybs), collapse = ", "), msg_bs
  ) %>%
    crayon::black() %>%
    cat()

  cat("\n")
}

#' Print out important arguments in smooth term te() or ti() or t2() in mgcv::gam() formula
#'
#' @details
#' Why a separate function is needed for t(), cannot using s(): in ofInterest,
#' "fx" is "fx" for t(), but "fixed" for s() - so they are different.
#' ref: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/te or /t2()
#'
#' @param FUN could be mgcv::te(), ti() or t2()
#' @param ofInterest got via: `gam.formula.breakdown <- mgcv::interpret.gam(formula)`;
#' `ofInterest <- gam.formula.breakdown$smooth.spec[[i]]`
#' @importFrom dplyr %>%
#' @noRd
#'
checker_gam_t <- function(FUN, ofInterest) {
  paste0(ofInterest$label, ": ") %>%
    crayon::black() %>%
    cat()

  # t()'s interpret.gam()'s smooth.spec does not have k or bs.dim as s() does;
  # may provided in ofInterest$margin[[i-xx]]$bs.dim but not fully clear

  ### fx:
  m1 <- invisible(eval(formals(FUN)[["fx"]])) %>% as.character() # default
  m2 <- ofInterest$fx # actual

  if (as.character(!as.logical(m1)) %in% m2) { # there is an opposite logical value in m2 (actual)
    msg_fx <- ""
  } else {
    msg_fx <- " (default)"
  }

  paste0("  fx = ", toString(ofInterest$fx), msg_fx, "; ") %>%
    crayon::black() %>%
    cat()

  ### bs:  # also different way of extracting from s()'s
  m1 <- invisible(eval(formals(FUN)[["bs"]])) %>% as.character() # default

  mybs <- list() # actual, as a list
  for (i in seq_along(ofInterest$margin)) {
    temp <- gsub(".smooth.spec", "", ofInterest$margin[[i]] %>% class())
    mybs[i] <- temp
  }

  # then check if all elements are default value of bs
  if ((length(unique(mybs)) == 1) && (m1 %in% unique(mybs))) { # all elements are the same, and = default
    msg_bs <- " (default)"
  } else {
    msg_bs <- ""
  }
  # print out
  paste0(
    "  bs = ",
    paste(as.character(mybs), collapse = ", "), msg_bs
  ) %>%
    crayon::black() %>%
    cat()

  cat("\n")
}

#' A checker for formula in gam for ModelArray.gam()
#' @param formula The formula
#' @param gam.formula.breakdown Got from mgcv::interpret.gam(formula)
#' @param onemodel The model of one element got from mgcv::gam()
#' @importFrom dplyr %>%
#' @noRd
#'
checker_gam_formula <- function(formula, gam.formula.breakdown, onemodel = NULL) {
  # print out the formula:
  temp <- formula %>% as.character()
  str_formula <- paste0(temp[2], " ~ ", temp[3])

  m1 <- paste0("The formula requested: ", str_formula) %>%
    crayon::black() %>%
    cat()
  cat(m1, "\n")


  if (length(gam.formula.breakdown$smooth.spec) != 0) { # if there is smooth term
    list_smooth_terms <- character(length(gam.formula.breakdown$smooth.spec))
    for (i_smoothTerm in seq_along(gam.formula.breakdown$smooth.spec)) {
      ofInterest <- gam.formula.breakdown$smooth.spec[[i_smoothTerm]]
      list_smooth_terms[i_smoothTerm] <- ofInterest$label

      # check what class of smooth term (s or ti or ?); then call the corresponding function
      # - throw out the message for important argument for this class
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
  } else { # no smooth term
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
#' Generates a formula in the format
#' \code{y ~ orderedFactor + s(x) + s(x, by = orderedFactor)},
#' where \code{y} is \code{response.var}, \code{x} is \code{smooth.var},
#' and \code{orderedFactor} is \code{factor.var}.
#' The formula generated can be further modified, e.g. by adding covariates.
#'
#' @details
#' This helper exists because setting up factor-smooth interactions in
#' \code{\link[mgcv]{gam}} requires an ordered factor and a specific formula
#' structure. If \code{factor.var} in \code{phenotypes} is not already an
#' ordered factor, this function creates one using \code{reference.group} as
#' the baseline level and adds it as a new column (named with
#' \code{prefix.ordered.factor} prepended to \code{factor.var}).
#'
#' The returned \code{phenotypes} data.frame must be used in the subsequent
#' \code{\link{ModelArray.gam}} call so that the ordered factor column is
#' available to the model.
#'
#' @param response.var Character. The variable name for the response
#'   (dependent variable), typically a scalar name like \code{"FD"}.
#' @param factor.var Character. The variable name for the factor. It should
#'   be an ordered factor in \code{phenotypes}. If not, an ordered factor
#'   will be generated as a new column, which requires \code{reference.group}.
#' @param smooth.var Character. The variable name for the smooth term main
#'   effect (e.g. \code{"age"}).
#' @param phenotypes A data.frame of the cohort with columns of independent
#'   variables, including \code{factor.var} and \code{smooth.var}.
#' @param reference.group Character. The reference (baseline) group for the
#'   ordered factor of \code{factor.var}. Required when \code{factor.var}
#'   in \code{phenotypes} is not already an ordered factor.
#' @param prefix.ordered.factor Character. Prefix for the ordered factor
#'   column name. Required when \code{factor.var} in \code{phenotypes} is
#'   not already an ordered factor. Default is \code{"o"}.
#' @param fx Logical. Passed to \code{\link[mgcv]{s}}. If \code{TRUE}
#'   (recommended), the smooth is treated as fixed degrees of freedom.
#'   Default is \code{TRUE}.
#' @param k Integer or \code{NULL}. Basis dimension passed to
#'   \code{\link[mgcv]{s}} for both the main smooth and interaction terms.
#'   If \code{NULL} (default), uses the default from \code{mgcv::s()}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{formula}{The generated \code{\link[stats]{formula}} object.}
#'     \item{phenotypes}{The (possibly updated) data.frame. If
#'       \code{factor.var} was not already an ordered factor, a new column
#'       named \code{paste0(prefix.ordered.factor, factor.var)} is added.
#'       Otherwise identical to the input.}
#'   }
#'
#' @seealso \code{\link{gen_gamFormula_contIx}} for continuous-by-continuous
#'   interactions, \code{\link{ModelArray.gam}} which accepts the generated
#'   formula.
#'
#' @examples
#' \dontrun{
#' phenotypes <- read.csv("cohort.csv")
#'
#' # factor.var is not yet ordered - function creates it
#' result <- gen_gamFormula_fxSmooth(
#'   response.var = "FD",
#'   factor.var = "sex",
#'   smooth.var = "age",
#'   phenotypes = phenotypes,
#'   reference.group = "female"
#' )
#' result$formula
#'
#' # Use the updated phenotypes (contains the ordered factor column)
#' results <- ModelArray.gam(
#'   result$formula,
#'   data = ma,
#'   phenotypes = result$phenotypes,
#'   scalar = "FD"
#' )
#' }
#'
#' @rdname gen_gamFormula_fxSmooth
#' @export
#'
gen_gamFormula_fxSmooth <- function(response.var, factor.var, smooth.var, phenotypes,
                                    reference.group = NULL, prefix.ordered.factor = "o",
                                    fx = TRUE, k = NULL) {
  class.factor.var <- class(phenotypes[[factor.var]])
  if (
    !(
      (length(class.factor.var) == 2) &&
        (class.factor.var[1] == "ordered") &&
        (class.factor.var[2] == "factor")
    )
  ) { # class is not c("ordered", "factor")

    message(
      "input `factor.var` is not an ordered factor; will generate",
      " one in data.frame `phenotypes` which will be returned"
    )
    if (is.null(reference.group)) {
      stop("requires a reference.group to generate the ordered factor")
    }

    # name of the ordered factor:
    unordered.factor.var <- factor.var
    factor.var <- paste0(prefix.ordered.factor, unordered.factor.var)

    message(paste0("the ordered factor will be named as: ", factor.var))

    # check if factor.var already exists:
    if (factor.var %in% colnames(phenotypes)) {
      stop(
        paste0(
          "a column with the same name '",
          factor.var,
          "' already exists in `phenotypes` data frame! Please change the `prefix.ordered.factor`!"
        )
      )
    }

    # add the column to phenotypes:
    list.groups <- unique(phenotypes[[unordered.factor.var]])
    list.groups <- list.groups[list.groups != reference.group] # temporarily drop the reference group
    list.groups <- c(reference.group, list.groups) # add as first
    phenotypes[[factor.var]] <- ordered(phenotypes[[unordered.factor.var]], levels = list.groups)
    # the first element in the list.groups would be the reference group
  }

  if (is.null(k)) {
    k <- invisible(eval(formals(mgcv::s)[["k"]]))
  }

  # generate the formula:
  formula <- paste0(response.var, "~", factor.var, "+")
  formula <- paste0(formula, "s(", smooth.var, ",k=", toString(k), ",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "s(", smooth.var, ",k=", toString(k), ",by=", factor.var, ",fx=", toString(fx), ")")

  # return
  formula <- stats::as.formula(formula) # when printing formula, it could be shown in more than one lines...
  toReturn <- list(
    formula = formula,
    phenotypes = phenotypes
  )
  return(toReturn)
}


#' Generate GAM formula with continuous-by-continuous interaction
#'
#' @description
#' Generates a formula in the format
#' \code{y ~ ti(x) + ti(z) + ti(x, z)}, where \code{y} is
#' \code{response.var}, \code{x} is \code{cont1.var}, and \code{z} is
#' \code{cont2.var}. The formula generated can be further modified, e.g.
#' by adding covariates.
#'
#' @details
#' This helper uses \code{\link[mgcv]{ti}} (tensor product interaction)
#' terms so that the interaction \code{ti(x, z)} captures only the
#' interaction effect, separate from the main effects \code{ti(x)} and
#' \code{ti(z)}. This decomposition is important for interpretability and
#' for requesting \code{changed.rsq.term.index} in
#' \code{\link{ModelArray.gam}}.
#'
#' @param response.var Character. The variable name for the response
#'   (dependent variable), typically a scalar name like \code{"FD"}.
#' @param cont1.var Character. The name of the first continuous variable.
#' @param cont2.var Character. The name of the second continuous variable.
#' @param fx Logical. Passed to \code{\link[mgcv]{ti}}. If \code{TRUE}
#'   (recommended), the smooth is treated as fixed degrees of freedom.
#'   Default is \code{TRUE}.
#' @param k Integer or \code{NULL}. Basis dimension passed to
#'   \code{\link[mgcv]{ti}} for all three terms (both main effects and the
#'   interaction). If \code{NULL} (default), uses the default from
#'   \code{mgcv::ti()}.
#'
#' @return A \code{\link[stats]{formula}} object.
#'
#' @seealso \code{\link{gen_gamFormula_fxSmooth}} for factor-smooth
#'   interactions, \code{\link{ModelArray.gam}} which accepts the generated
#'   formula.
#'
#' @examples
#' \dontrun{
#' formula <- gen_gamFormula_contIx(
#'   response.var = "FD",
#'   cont1.var = "age",
#'   cont2.var = "cognition"
#' )
#' formula
#'
#' # Use in ModelArray.gam with changed R-squared for the interaction
#' results <- ModelArray.gam(
#'   formula,
#'   data = ma,
#'   phenotypes = phenotypes,
#'   scalar = "FD",
#'   changed.rsq.term.index = list(3)
#' )
#' }
#'
#' @rdname gen_gamFormula_contIx
#' @export
gen_gamFormula_contIx <- function(response.var, cont1.var, cont2.var,
                                  fx = TRUE, k = NULL) {
  if (is.null(k)) {
    k <- invisible(eval(formals(mgcv::ti)[["k"]]))
  }

  formula <- paste0(response.var, "~")
  formula <- paste0(formula, "ti(", cont1.var, ",k=", toString(k), ",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "ti(", cont2.var, ",k=", toString(k), ",fx=", toString(fx), ")", "+")
  formula <- paste0(formula, "ti(", cont1.var, ",", cont2.var, ",k=", toString(k), ",fx=", toString(fx), ")")

  formula <- stats::as.formula(formula)
  return(formula)
}


#' Bind two tibbles together, considering one or both of them is empty tibble()
#' @details Without this function, bind_cols(A, tibble()) will becomes tibble()...
#' @param a A tibble, can be empty tibble()
#' @param b A tibble, can be empty tibble()
#' @return c, A tibble after binding a and b together
#' @noRd
bind_cols_check_emptyTibble <- function(a, b) {
  flag_a_empty <- !all(dim(a)) # if TRUE, a is empty
  flag_b_empty <- !all(dim(b)) # if TRUE, b is empty

  if (flag_a_empty && flag_b_empty) c <- tibble::tibble() # both are empty
  if (flag_a_empty && (!flag_b_empty)) c <- b # b is not empty ==> taking b
  if ((!flag_a_empty) && flag_b_empty) c <- a # a is not empty ==> taking a
  if ((!flag_a_empty) && (!flag_b_empty)) c <- dplyr::bind_cols(a, b) # neither of them is empty ==> simply combine

  c
}


#' Summarize an HDF5 file without loading a full ModelArray
#'
#' @description
#' Reads the HDF5 file structure and returns a summary of available scalars,
#' their dimensions, and any saved analyses. Useful for inspecting large files
#' without constructing a full \linkS4class{ModelArray} object.
#'
#' @details
#' This function opens the HDF5 file read-only via \code{\link[rhdf5]{h5ls}},
#' inspects the group structure under \code{/scalars/} and \code{/results/},
#' and closes the file. It does not load any data into memory. The returned
#' object has a \code{print} method that displays a formatted summary.
#'
#' @param filepath Character. Path to an HDF5 (\code{.h5}) file.
#'
#' @return An object of class \code{"h5summary"}, which is a list with
#'   components:
#'   \describe{
#'     \item{scalars}{A data.frame with columns \code{name},
#'       \code{nElements}, and \code{nInputFiles}.}
#'     \item{analyses}{Character vector of analysis names found under
#'       \code{/results/}.}
#'     \item{filepath}{The input filepath.}
#'   }
#'
#' @seealso \code{\link{ModelArray}} for loading the full object,
#'   \linkS4class{ModelArray} for the class definition.
#'
#' @examples
#' \dontrun{
#' h5summary("path/to/data.h5")
#'
#' # Inspect before deciding which scalars to load
#' info <- h5summary("path/to/data.h5")
#' info$scalars$name
#' ma <- ModelArray("path/to/data.h5", scalar_types = info$scalars$name)
#' }
#'
#' @rdname h5summary
#' @export
h5summary <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  listing <- rhdf5::h5ls(filepath)
  rhdf5::h5closeAll()

  # Find scalar datasets: /scalars/<name>/values
  scalar_rows <- listing[
    listing$group != "/scalars" &
      grepl("^/scalars/", listing$group) &
      !grepl("^/scalars/scalars", listing$group) &
      listing$name == "values",
  ]
  scalar_info <- data.frame(
    name = sub("^/scalars/", "", scalar_rows$group),
    dim = scalar_rows$dim,
    stringsAsFactors = FALSE
  )

  # Parse dimensions
  if (nrow(scalar_info) > 0) {
    dims <- strsplit(scalar_info$dim, " x ")
    scalar_info$nElements <- as.integer(sapply(dims, `[`, 1))
    scalar_info$nInputFiles <- as.integer(sapply(dims, `[`, 2))
    scalar_info$dim <- NULL
  }

  # Find analyses: /results/<name>
  result_groups <- listing[
    listing$group == "/results" & listing$otype == "H5I_GROUP",
  ]
  analyses <- result_groups$name

  structure(
    list(
      scalars = scalar_info,
      analyses = analyses,
      filepath = filepath
    ),
    class = "h5summary"
  )
}

#' @rdname h5summary
#'
#' @param x An \code{h5summary} object as returned by
#'   \code{\link{h5summary}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisible \code{x}. Called for its side effect of printing a
#'   human-readable summary to the console.
#'
#' @method print h5summary
#' @export
print.h5summary <- function(x, ...) {
  cat("H5 file:", x$filepath, "\n\n")
  if (nrow(x$scalars) > 0) {
    cat("Scalars:\n")
    for (i in seq_len(nrow(x$scalars))) {
      cat("  ", x$scalars$name[i], ": ",
        x$scalars$nElements[i], " elements x ",
        x$scalars$nInputFiles[i], " input files\n",
        sep = ""
      )
    }
  } else {
    cat("  (no scalars found)\n")
  }
  cat("\nAnalyses:", if (length(x$analyses) > 0) {
    paste(x$analyses, collapse = ", ")
  } else {
    "(none)"
  }, "\n")
  invisible(x)
}
