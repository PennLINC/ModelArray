% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyse.R
\name{ModelArray.lm}
\alias{ModelArray.lm}
\title{Fit linear model for element-wise data}
\usage{
ModelArray.lm(
  formula,
  data,
  phenotypes,
  scalar,
  element.subset = NULL,
  full.outputs = FALSE,
  var.terms = c("estimate", "statistic", "p.value"),
  var.model = c("adj.r.squared", "p.value"),
  correct.p.value.terms = c("fdr"),
  correct.p.value.model = c("fdr"),
  num.subj.lthr.abs = 10,
  num.subj.lthr.rel = 0.2,
  verbose = TRUE,
  pbar = TRUE,
  n_cores = 1,
  ...
)
}
\arguments{
\item{formula}{Formula (passed to `stats::lm()`)}

\item{data}{ModelArray class}

\item{phenotypes}{A data.frame of the cohort with columns of independent variables
and covariates to be added to the model. It should contains a column called "source_file",
and this column should match to that in \code{data}.}

\item{scalar}{A character. The name of the element-wise scalar to be analysed}

\item{element.subset}{A list of positive integers (min = 1, max = number of elements).
The subset of elements you want to run. Default is `NULL`, i.e. requesting all elements in `data`.}

\item{full.outputs}{TRUE or FALSE, Whether to return full set of outputs.
If FALSE, it will only return those requested in arguments \code{var.*} and \code{correct.p.value.*};
if TRUE, arguments \code{var.*} will be ignored, and will return all possible statistics for
\code{var.*} and any options requested in arguments \code{correct.p.value.*}.}

\item{var.terms}{A list of characters.
The list of variables to save for terms (got from `broom::tidy()`).
See "Details" section for more.}

\item{var.model}{A list of characters.
The list of variables to save for the model (got from `broom::glance()`).
See "Details" section for more.}

\item{correct.p.value.terms}{A list of characters.
To perform and add a column for p.value correction for each term. Default: "fdr".
See "Details" section for more.}

\item{correct.p.value.model}{A list of characters.
To perform and add a column for p.value correction for the model.
Default: "fdr". See "Details" section for more.}

\item{num.subj.lthr.abs}{An integer, lower threshold of absolute number of subjects.
For an element, if number of subjects who have finite values (defined by `is.finite()`,
i.e. not NaN or NA or Inf) in h5 file > \code{num.subj.lthr.abs},
then this element will be run normally; otherwise,
this element will be skipped and statistical outputs will be set as NaN.
Default is 10.}

\item{num.subj.lthr.rel}{A value between 0-1, lower threshold of relative number of subjects.
Similar to \code{num.subj.lthr.abs},
if proportion of subjects who have valid value > \code{num.subj.lthr.rel},
then this element will be run normally; otherwise,
this element will be skipped and statistical outputs will be set as NaN.
Default is 0.2.}

\item{verbose}{TRUE or FALSE, to print verbose message or not}

\item{pbar}{TRUE or FALSE, to print progress bar or not}

\item{n_cores}{Positive integer, The number of CPU cores to run with}

\item{...}{Additional arguments for `stats::lm()`}
}
\value{
Tibble with the summarized model statistics for all elements requested
}
\description{
`ModelArray.lm` fits linear model (`stats::lm()`) for each of elements requested, and returns a tibble
dataframe of requested model statistics.
}
\details{
You may request returning specific statistical variables by setting \code{var.*},
or you can get all by setting \code{full.outputs=TRUE}.
Note that statistics covered by \code{full.outputs} or \code{var.*} are the ones from broom::tidy()
and broom::glance() only, and do not include corrected p-values.
However FDR-corrected p-values ("fdr") are generated by default.
List of acceptable statistic names for each of \code{var.*}:
\itemize{
  \item \code{var.terms}: c("estimate","std.error","statistic","p.value");
     For interpretation please see \link[broom]{tidy.lm}.
  \item \code{var.model}: c("r.squared", "adj.r.squared", "sigma", "statistic",
     "p.value", "df", "logLik", "AIC", "BIC", "deviance", "df.residual", "nobs");
     For interpretation please see \link[broom]{glance.lm}.
}
For p-value corrections (arguments \code{correct.p.value.*}),
supported methods include all methods in `p.adjust.methods` except "none".
Can be more than one method. FDR-corrected p-values ("fdr") are calculated by default.
Turn it off by setting to "none". \cr
Arguments \code{num.subj.lthr.abs} and \code{num.subj.lthr.rel}
are mainly for input data with subject-specific masks,
 i.e. currently only for volume data. For fixel-wise data, you may ignore these arguments.
}
