rm(list=ls())

library(dplyr)
library(mgcv)

source("R/utils-pipe.R")
source("R/analyse.R")
source("R/utils.R")
source("R/ModelArray_Constructor.R")
source("R/ModelArray_S4Methods.R")

h5_path <- system.file("extdata", "n50_fixels.h5", package = "ModelArray")

scalar_name <- c("FD")
modelarray <- ModelArray(h5_path,
                         scalar_types = scalar_name,
                         analysis_names = c("my_analysis"))

csv_path <- system.file("extdata", "n50_cohort.csv", package = "ModelArray")
phenotypes <- read.csv(csv_path)
element.subset = 1:10

var.smoothTerms = c("statistic","p.value")
var.parametricTerms = c("estimate", "statistic", "p.value")
var.model = c("dev.expl", "adj.r.squared")

full.formula <- FD ~ s(age) + factorA
reduced.formula <- FD ~ factorA

results <- analyseOneElement.gam.fullNred(i_element=1, full.formula, reduced.formula, # here is the only diff with analyseOneElement.gam()'s arguments is: there are two input formula for this function
                               modelarray, phenotypes, scalar = scalar_name, 
                               var.smoothTerms, var.parametricTerms, var.model,
                               flag_initiate = TRUE)



# mygam <- ModelArray.gam(FD ~  s(age) + factorA, data = modelarray, phenotypes = phenotypes, scalar = scalar_name, element.subset = element.subset,
#                var.model = c("dev.expl", "adj.r.squared"),
#                changed.rsq.term.index = c(1),
#                n_cores = 2, pbar = FALSE)
