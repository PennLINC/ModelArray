" This script is to analyze the fixel cluster of significance.
Steps:
  // python file = ConFixel/notebooks/analyze_signi_cluster.py
1. [python]: function: convert_voxelMask_to_fixelIndex()
  This will generate the list of fixel's ids that included in manually drawn mask (ROI.mif)
  
2. [R] All significant fixels: the list 
    - load in the .h5 file (with output such as p.value after bonferroni correction)
    - threshold, get the list of all significant fixels
3. [R] take the intersect of 1 and 2; save it.

4. [python]: function: save_selectedFixel_mask()
  This is to verify that the intersect list we saved is visually correct. 
  The function will save a 'mask' of intersected fixels we just saved.
  You should visually verify it in mrview. See instructions at the end of the function.
  
5. [R] Analyze the intersect


"



rm(list=ls())
# set up
library("dplyr")  # for %>%
library("mgcv")
library("broom")
library("testthat")
library("gratia")
library("ggplot2")

source("notebooks/GAMM_plotting.R")   # Bart Larsen's function for visualizing gam results
source("R/FixelArray_Constructor.R")
source("R/FixelArray_S4Methods.R")
source("R/utils.R")
source("R/analyse.R")

step2_thresholding <- function(fn.h5.results, scalar_name, analysis_name, stat_name_thr, thr, flag_compare, folder.h5.results,
                               results_matrix,
                               flag_run_step2) {
  fn.fixel_id_list_thr <- paste0(folder.h5.results, "/", analysis_name, "_", stat_name_thr, "_",
                                 flag_compare, "_", toString(thr),
                                 "_fixelIdList.txt")
  
  if (flag_run_step2 == TRUE) {
    # after thresholding, the fixel_id list:
    if (flag_compare == "gt") {
      fixel_id_list_thr <- results_matrix[(results_matrix[,stat_name_thr] > thr), "fixel_id"]  
    } else if (flag_compare == "lt") {
      fixel_id_list_thr <- results_matrix[(results_matrix[,stat_name_thr] < thr), "fixel_id"]  
    } else {
      stop("invalid flag_compare!")
    }
    
    message(paste0("found ", toString(length(fixel_id_list_thr)), " fixels after thresholding"))   # number of fixels  
    
    # save the list of fixel_id that reaches the threshold
    write.table(fixel_id_list_thr, fn.fixel_id_list_thr, row.names=FALSE, col.names=FALSE, quote = FALSE)
  } 
  
  return(fn.fixel_id_list_thr)
}

step3_intersect <- function(folder.h5.results, 
                           filename.fixelIdListMask, fn.fixel_id_list_thr,
                           flag_flipFixel_roi, flag_flipFixel_signi, num_fixel_total,
                           flag_run_step3){
  # filename for saving the list of fixel ids of the intersection:
  fn.fixel_id_list_intersect <- gsub("_fixelIdList.txt",
                                     paste0("__Intersect__", filename.fixelIdListMask),
                                     fn.fixel_id_list_thr)
  
  if (flag_run_step3 == TRUE) {
    # load: list of fixels' ids after thresholding (without manually drawn mask)
    fixel_id_list_thr <- scan(fn.fixel_id_list_thr, what="", sep="\n") %>% as.integer()
    if (flag_flipFixel_signi == TRUE) {
      fixel_id_list_thr <- num_fixel_total - 1 - fixel_id_list_thr  # "-1": because the fixel id starts from 0
    }
    
    # load: list of fixels' ids within manually drawn mask:
    fn.fixelIdListMask <- file.path(folder.h5.results, filename.fixelIdListMask)
    fixel_id_list_mask <- scan(fn.fixelIdListMask, what = "", sep = "\n") %>% as.integer()
    if (flag_flipFixel_roi == TRUE) {
      fixel_id_list_mask <- num_fixel_total - 1 - fixel_id_list_mask  # "-1": because the fixel id starts from 0
    }
    
    # get the intersection:
    
    fixel_id_list_intersect <- intersect(fixel_id_list_thr, fixel_id_list_mask)
    length(fixel_id_list_intersect)
    
    # save the intersection:
    
    if (fn.fixel_id_list_intersect == fn.fixel_id_list_thr) {
      stop("The filename for fixel intersection is not correct and will overwrite another file!")
    }
    write.table(fixel_id_list_intersect, fn.fixel_id_list_intersect, row.names=FALSE, col.names=FALSE, quote = FALSE)
  } 
  
  return(fn.fixel_id_list_intersect)
}
  
  
### inputs: #####
num.subj <- 938
fn.h5.results <- paste0("/home/chenying/Desktop/fixel_project/data/data_from_josiane/results/ltn_FDC_n",toString(num.subj),"_wResults_nfixel-0_20211126-182543.h5")
fn_csv <- paste0("../data/data_from_josiane/df_example_n", toString(num.subj), ".csv")
scalar_name <- c("FDC")

analysis_name <- "gam_allOutputs"

stat_name_thr <- "s_Age.p.value.bonferroni"   # the stat name for thresholding
flag_compare <- "lt"
thr <- 1e-20

# stat_name_thr <- "s_Age.p.value"  # the stat name for thresholding
# flag_compare <- "lt"
# thr <- 1e-15

# stat_name_thr <- "s_Age.eff.size"
# flag_compare <- "gt"
# thr <- 0.2

## step 2: 
flag_run_step2 <- FALSE   # run once is enough; independent from python's output

flag_flipFixel_roi <- TRUE
flag_flipFixel_signi <- FALSE

## step 3:
flag_run_step3 <- FALSE
filename.fixelIdListMask <- "ROI_x65_sage_p_bonfer_lt_1e-20_fixelIdList.txt"  # for step 3
#filename.fixelIdListMask <- "ROI_x69_sage_p_lt_1e-15_fixelIdList.txt"  # for step 3

## step 5:
stat_toPlot <- "s_Age.eff.size"
formula <- FDC ~ s(Age, k=4, fx=TRUE) + sex
method.gam.refit <- "GCV.Cp"   # +++++++++++++++

### load data #####
folder.h5.results <- gsub(".h5", "", fn.h5.results, fixed=TRUE)
fixelarray <- FixelArray(fn.h5.results, scalar_types = scalar_name, analysis_names = analysis_name)
num_fixel_total <- nrow(fixelarray@fixels)
if (num.subj != fixelarray@subjects[[scalar_name]] %>% length()) {
  stop("number of subjects in fixelarray is not equal to requested one!")
}
results_matrix <- fixelarray@results[[analysis_name]]$results_matrix 
# colnames(fixelarray@results$gam_allOutputs$results_matrix )

phenotypes <- read.csv(fn_csv)
# check # subjects matches:
if (nrow(phenotypes) != num.subj) {
  stop(paste0("number of subjects in .csv = ", toString(nrow(phenotypes)), ", is not equal to entered number = ", toString(num.subj)))
}

### Step 2: Thresholding #####

fn.fixel_id_list_thr <- step2_thresholding(fn.h5.results, scalar_name, analysis_name, stat_name_thr, thr, flag_compare, folder.h5.results, 
                     results_matrix, flag_run = flag_run_step2)



### Step 3: get intersection #####
fn.fixel_id_list_intersect <- step3_intersect(folder.h5.results, 
                           filename.fixelIdListMask, fn.fixel_id_list_thr,
                           flag_flipFixel_roi, flag_flipFixel_signi, num_fixel_total,
                           flag_run_step3)

### Step 4: Please verify selected fixel ids! See python file. #####

### Step 5: Average and plot #####
# load the final list:
fixel_id_list_intersect <- scan(fn.fixel_id_list_intersect, what="", sep="\n") %>% as.integer()

# avg
scalar_matrix <- scalars(fixelarray)[[scalar_name]]
if (nrow(scalar_matrix) != num_fixel_total) {
  stop("scalar_matrix does not contain full list of fixels!")
}
matrix_selected <- scalar_matrix[fixel_id_list_intersect+1, ]    # # of selected fixels x # of subjects | !! fixel_id starts from 0 so need to +1 !!!


# double check they are the "selected" fixels: meeting the criteria when selecting
dat_selectedFixels_metric <- data.frame(fixel_id = fixel_id_list_intersect,
                                        selecting_metric = numeric(length(fixel_id_list_intersect)),
                                        s_Age_p.value = numeric(length(fixel_id_list_intersect)))   # all zeros
for (i_fixel_selected in 1:length(fixel_id_list_intersect)) {
  # re-fit:
  fixel_id <- fixel_id_list_intersect[i_fixel_selected]
  
  values <- scalars(fixelarray)[[scalar_name]][(fixel_id + 1),]    # fixel_id starts from 0
  
  dat <- phenotypes
  dat[[scalar_name]] <- values
  
  onemodel <- mgcv::gam(formula = formula, data = dat,
                        method = method.gam.refit)
  onemodel.tidy.smoothTerms <- onemodel %>% broom::tidy(parametric = FALSE)
  onemodel.tidy.parametricTerms <- onemodel %>% broom::tidy(parametric = TRUE)
  onemodel.glance <- onemodel %>% broom::glance()
  onemodel.summary <- onemodel %>% summary()
  
  temp <- results_matrix[fixel_id + 1, stat_name_thr]   # fixel_id starts from 0
  dat_selectedFixels_metric[i_fixel_selected, "selecting_metric"] <- temp   # from results_matrix
  
  dat_selectedFixels_metric[i_fixel_selected, "s_Age_p.value"] <- onemodel.tidy.smoothTerms$p.value
}


print("max p.value after bonferroni:")
dat_selectedFixels_metric$s_Age_p.value %>% max() * num_fixel_total
if (flag_compare == "lt") {
  expect_true(max(dat_selectedFixels_metric$selecting_metric) < thr)
} else if (flag_compare == "gt") {
  expect_true(min(dat_selectedFixels_metric$selecting_metric) > thr)
}
  
# # if selecting_metric == s_Age.p.value.bonferroni
# testthat::expect_equal(dat_selectedFixels_metric$s_Age_p.value * num_fixel_total,
#                        dat_selectedFixels_metric$selecting_metric)   


### take the average within the cluster (selected fixels) ########
# # averaged across fixels x # of subjects:
#avgFixel_subj = list of number of subjects
# then loop across subjects (columns), get the avg 

avgFixel_subj <- numeric(num.subj)
for (i_subj in 1:num.subj) {
  # take the average:
  avgFixel_subj[i_subj] <- mean(matrix_selected[,i_subj])
}
df_avgFixel <- phenotypes
df_avgFixel[[scalar_name]] <- avgFixel_subj

### plot ######

#' @param fixel_id starting from 0!
plot_oneFixel <- function(fixelarray, fixel_id, scalar_name,
                          phenotypes,
                          dat = NULL, return_else = FALSE) {
  if (is.null(dat)) {
    values <- scalars(fixelarray)[[scalar_name]][(fixel_id + 1),]    # fixel_id starts from 0
    
    dat <- phenotypes
    dat[[scalar_name]] <- values
    
  } else {
    # directly using dat
    
  }
  
  onemodel <- mgcv::gam(formula = formula, data = dat,
                        method = method.gam.refit)
  
  
  #f <- vis.gam(onemodel)
    
  f <- visualize_model(onemodel, smooth_var = 'Age')
  
  ## below: tried but there is an offset in y (Bart: y is centered at 0)
  #f <- gratia::draw(onemodel, select = "s(Age)", residuals=TRUE, rug = FALSE)   # rug is at the bottom - for displaying the data density on x-axis
  #f + theme_classic()
  
  if (return_else == FALSE) {
    return(f)
  } else {
    results = list(f = f,
                   onemodel = onemodel)
    return(results)
  }
  
  
}

# plot one fixel:
f_1 <- plot_oneFixel(fixelarray, fixel_id_list_intersect[1], scalar_name, phenotypes)
f_last <- plot_oneFixel(fixelarray, fixel_id_list_intersect[length(fixel_id_list_intersect)], scalar_name, phenotypes)

results <- plot_oneFixel(fixelarray=NULL, NULL, scalar_name, phenotypes, dat=df_avgFixel, return_else = TRUE)
f_avgFixel <- results$f
onemodel_avgFixel <- results$onemodel
f_avgFixel

onemodel_avgFixel.summary <- summary(onemodel_avgFixel)
onemodel_avgFixel.smoothTerm <- broom::tidy(onemodel_avgFixel, parametric=FALSE)
onemodel_avgFixel.parametricTerm <- broom::tidy(onemodel_avgFixel, parametric=TRUE) 
onemodel_avgFixel.model <- broom::glance(onemodel_avgFixel)

red.formula <- formula(drop.terms(terms(formula, keep.order = TRUE), 
                                  c(1), keep.response = TRUE))   # drop the first term, i.e. smooth of age
redmodel_avgFixel <- mgcv::gam(formula=red.formula, data = dat,
                               method = method.gam.refit)
redmodel_avgFixel.summary <- summary(redmodel_avgFixel)
eff.size.avgFixel <- onemodel_avgFixel.summary$r.sq - redmodel_avgFixel.summary$r.sq  

s_Age.p.value_avgFixel <- onemodel_avgFixel.smoothTerm[onemodel_avgFixel.smoothTerm$term=="s(Age)","p.value"]

onemodel_avgFixel.summary
print(paste0("s(Age)'s p.value of re-fit after avg in this cluster = ", toString(s_Age.p.value_avgFixel)))  # if =0, it means <1e-16
print(paste0("s(Age)'s effect size of re-fit after avg in this cluster = ", sprintf("%.3f",eff.size.avgFixel)))  

# and add to the plot!
x = 12; y = 1.65
# x = 20; y = 0.55
if (s_Age.p.value_avgFixel < 0.001) {
  txt.s_Age.p.value_avgFixel = "s(Age)'s p.value < 0.001"
} else {
  txt.s_Age.p.value_avgFixel = paste0("s(Age)'s p.value = ", sprintf("%.3f",s_Age.p.value_avgFixel))
}

 
f_avgFixel + geom_text(x=x, y=y, size = 6,
                       label=paste0(txt.s_Age.p.value_avgFixel, "\n",
                                    "s(Age)'s effect size = ", sprintf("%.3f",eff.size.avgFixel)))

# NOTE: as there is more sex=2 than sex=1, and sex is numeric, so the median(df[,sex]) = 2, the gam curve is fitted upon sex=2
# TODO: check how many fixels' model: sex is significant

