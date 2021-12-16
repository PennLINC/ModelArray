
# HOW TO RUN:
# in bash, same folder as this current file:
# $ Rscript ./memoryProfiling_ModelArray.gam.R  > xxx.txt 2>&1 &


# set ups
rm(list = ls())

# the working directory of this Demo is where it locates

flag_library_what <- "manually"   # "ModelArray" or "manually"

if (flag_library_what == "ModelArray") {
  message("run: devtools::install() to install ModelArray package")
  library(devtools)
  devtools::install()
  library(ModelArray)
  
  library(tictoc)
  
} else if (flag_library_what == "manually") {
  
  message("run: source several R scripts and library some R packages...")
  
  source("../R/ModelArray_Constructor.R")
  source("../R/ModelArray_S4Methods.R")
  source("../R/utils.R")
  source("../R/analyse.R")
  # library(ModelArray)
  suppressMessages(library(dplyr))
  library(broom)
  library(hdf5r)
  library(tictoc)
  library(mgcv)
  # library(lineprof)
  # library(profvis)
  # library(peakRAM)
  suppressMessages(library(doParallel))
}


flag_whichdataset <- "josiane"   # "val" or "test_n50" or "josiane"
num.subj <- 938  # [integer]   
num.fixels <- 0  # 0 = full 
flag_which_subset <- ""
flag_where <- "vmware"   # "CUBIC" or "vmware"


#####
now_str <- format(Sys.time(), "%Y%m%d-%H%M%S")

if (flag_whichdataset == "test_n50") {
  fn <- "../inst/extdata/n50_fixels.h5"
  
  if (flag_where == "CUBIC") {
    fn.output <- "../../dropbox/data_forCircleCI_n50/n50_fixels_output.h5"
  } else if (flag_where == "vmware") {
    fn.output <- "/home/chenying/Desktop/fixel_project/data/data_forCircleCI_n50/n50_fixels_output.h5"  
    # absoluate path: "/home/chenying/Desktop/fixel_project/data/data_forCircleCI_n50/n50_fixels_output.h5";  
    # relative path: "../../data/data_forCircleCI_n50/n50_fixels_output.h5"
  }
  
  fn_csv <- "../inst/extdata/n50_cohort.csv"
  
  scalar = c("FD")
  
} else if (flag_whichdataset == "josiane") {
  if (flag_where == "vmware") {
    fn <- paste0("../../data/data_from_josiane/ltn_FDC_n", toString(num.subj), ".h5")
    fn.output <-  paste0("../../data/data_from_josiane/results/ltn_FDC_n", toString(num.subj), "_wResults_nfixel-",toString(num.fixels), "_",now_str, ".h5")
    fn_csv <- paste0("../../data/data_from_josiane/df_example_n", toString(num.subj), ".csv")
    
  }
  
  scalar <- c("FDC")
}


# generate fn.output:
if (fn != fn.output) {
  file.copy(from=fn, to=fn.output, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)   # , recursive = TRUE
}

# h5closeAll()
fixelarray <- ModelArray(fn.output, scalar_types = scalar)


#fixelarray
#scalars(fixelarray)[[scalar]]

#####
# check # subjects matches:
if (dim(scalars(fixelarray)[[scalar]])[2] != num.subj) {
  stop(paste0("number of subjects in .h5 = ", dim(scalars(fixelarray)[[scalar]])[2], ", is not equal to entered number = ", toString(num.subj)))
}  
  
phenotypes <- read.csv(fn_csv)
# check # subjects matches:
if (nrow(phenotypes) != num.subj) {
  stop(paste0("number of subjects in .csv = ", toString(nrow(phenotypes)), ", is not equal to entered number = ", toString(num.subj)))
}

if (flag_whichdataset == "test_n50") {
  formula <- FD ~ s(age, k=4, fx=TRUE) + s(factorA)         # ++++ to do: add motion quantification   # FD ~ s(age, k=4) + sex  # FD ~ s(age) + sex
} else if (flag_whichdataset == "val") {
  formula <- FD ~ Age
} else if (flag_whichdataset == "josiane") {
  formula <- FDC ~ s(Age, k=4, fx = TRUE) + sex
}

gam.method = "REML"   # "GCV.Cp", "REML"  # any other methods usually used?
gam.optimizer = c("outer","newton")   # default: c("outer","newton") | # number of optimizers here will not change number of rows or columns in output gam model

if (num.fixels == 0) {
  num.fixels <- dim(scalars(fixelarray)[[scalar]])[1]
}
element.subset <- 1:num.fixels  

### running on real data #####
print(formula)
tic("Running ModelArray.gam()")
# ++++++++++++++= NEXT TIME: include method = gam.method!!! ++++++++++++++++++++=
# +++++++++++++++ NEXT TIME: sex --> ordered factor, and use oSex in formula! (this may make the plots - e.g. Bart's function more making sense? as there will be a reference level of female or male)++++++++++++++++++++++++++
gam_real <- ModelArray.gam(formula = formula, data = fixelarray, phenotypes = phenotypes, scalar = scalar, 
                           element.subset = element.subset, full.outputs = TRUE,
                           eff.size.term.index = c(1),
                           correct.p.value.smoothTerms = c("fdr", "bonferroni"),
                           correct.p.value.parametricTerms = c("fdr", "bonferroni"),
                           n_cores=4, pbar = TRUE)
toc(log = TRUE)
message("")


head(gam_real)
# write:
analysis_name <- "gam_allOutputs"
writeResults(fn.output, df.output = gam_real, analysis_name=analysis_name, overwrite=TRUE)

# read and see
fixelarray_new <- ModelArray(fn.output, scalar_types = scalar, analysis_names = analysis_name)
message("after saving to .h5:")
fixelarray_new@results$gam_allOutputs

message("after saving to .h5, head of gam_real:")
head(gam_real)
# 
# message("dimension of gam_real:")
# dim(gam_real)
