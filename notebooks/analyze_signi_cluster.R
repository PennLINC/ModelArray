" This script is to analyze the fixel cluster of significance.
Steps:
1. Run python file: ConFixel/notebooks/analyze_signi_cluster.py
  This will generate the list of fixel's ids that included in manually drawn mask (ROI.mif)
  
Then, run this R script:
2. All significant fixels: the list 
    - load in the .h5 file (with output such as p.value after bonferroni correction)
    - threshold, get the list of all significant fixels
3. take the intersect of 1 and 2; save it.
4. Analyze the intersect
"



rm(list=ls())
# set up
library("dplyr")  # for %>%

source("R/FixelArray_Constructor.R")
source("R/FixelArray_S4Methods.R")
source("R/utils.R")
source("R/analyse.R")

step2_thresholding <- function(fn.h5.results, scalar_name, analysis_name, stat_name_thr, thr, flag_compare, folder.h5.results,
                               results_matrix,
                               flag_run) {
  fn.fixel_id_list_thr <- paste0(folder.h5.results, "/", analysis_name, "_", stat_name_thr, "_",
                                 flag_compare, "_", toString(thr),
                                 "_fixelIdList.txt")
  
  if (flag_run == TRUE) {
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


### inputs: #####
fn.h5.results <- "/home/chenying/Desktop/fixel_project/data/data_from_josiane/results/ltn_FDC_n938_wResults_nfixel-0_20211126-182543.h5"
scalar_name <- c("FDC")

analysis_name <- "gam_allOutputs"
# stat_name_thr <- "s_Age.p.value.bonferroni"  # the stat name for thresholding
# flag_compare <- "lt"
# thr <- 1e-20

stat_name_thr <- "s_Age.eff.size"
flag_compare <- "gt"
thr <- 0.2

flag_run_step2 <- TRUE   # run once is enough; independent from python's output

filename.fixelIdListMask <- "ROI_x65_sage_p_bonfer_lt_1e-20_fixelIdList.txt"  # for step 3

folder.h5.results <- gsub(".h5", "", fn.h5.results, fixed=TRUE)

### load the output .h5 data #####

fixelarray <- FixelArray(fn.h5.results, scalar_types = scalar_name, analysis_names = analysis_name)
results_matrix <- fixelarray@results[[analysis_name]]$results_matrix 
# colnames(fixelarray@results$gam_allOutputs$results_matrix )

### Step 2: Thresholding #####

fn.fixel_id_list_thr <- step2_thresholding(fn.h5.results, scalar_name, analysis_name, stat_name_thr, thr, flag_compare, folder.h5.results, 
                     results_matrix, flag_run = flag_run_step2)



### Step 3: get intersection #####
# load: list of fixels' ids after thresholding (without manually drawn mask)
fixel_id_list_thr <- scan(fn.fixel_id_list_thr, what="", sep="\n") %>% as.integer()

# load: list of fixels' ids within manually drawn mask:

fn.fixelIdListMask <- file.path(folder.h5.results, filename.fixelIdListMask)
fixel_id_list_mask <- scan(fn.fixelIdListMask, what = "", sep = "\n") %>% as.integer()
 
# get the intersection:

fixel_id_list_intersect <- intersect(fixel_id_list_thr, fixel_id_list_mask)
length(fixel_id_list_intersect)
