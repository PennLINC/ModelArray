# preparation for Josiane's dataset

source("../R/ModelArray_Constructor.R")
source("../R/ModelArray_S4Methods.R")
source("../R/utils.R")
source("../R/analyse.R")

suppressMessages(library(dplyr))
library(broom)
library(rhdf5)
library(hdf5r)
library(tictoc)
library(testthat)

h5closeAll()


folder_josiane <- "/cbica/projects/fixel_db/dropbox/data_from_josiane"

### check .csv ##### 
df_example <- read.csv(paste0(folder_josiane, "/", "df_example_n938.csv"))
ltn_FDC_n941 <- read.csv(paste0(folder_josiane, "/", "ltn_FDC_n941.csv")) 

a <- df_example$bblid
b <- ltn_FDC_n941$subject

all(diff(a) > 0)   # all increasing
all(diff(b) > 0)   # all increasing

a_in_b <- a %in% b
setdiff_a_b <- setdiff(a,b)
setdiff_b_a <- setdiff(b,a)

setdiff_b_a_ids <- match(setdiff_b_a, b)  # the ids in b, whose values a does not have

a_prime <- append(a,       setdiff_b_a[1], setdiff_b_a_ids[1]-1)
a_prime <- append(a_prime, setdiff_b_a[2], setdiff_b_a_ids[2]-1)
a_prime <- append(a_prime, setdiff_b_a[3], setdiff_b_a_ids[3]-1)

all(diff(a_prime) > 0) 
expect_equal(a_prime,b)

ids_tokeep <- setdiff(1:941, setdiff_b_a_ids)


### prepare full .h5 file n=941-3 = 938 #####
# by removing 3 subjects from ltn_FDC.h5
fn_fa_n941 <- paste0(folder_josiane, "/", "ltn_FDC_n941.h5")
h5ls(fn_fa_n941)
fa_n941 <- ModelArray(fn_fa_n941, scalar_types=c("FDC"))   # number of subjects: 941
# n941.h5 <- H5File$new(fn_fa_n941, mode="a")
# n941.fixels_ds <- n941.h5[["fixels"]]

# create a new .h5 file:
fn_fa_n938 <- paste0(folder_josiane, "/", "ltn_FDC_n938.h5")
n938.h5 <- H5File$new(fn_fa_n938, mode="a")

# fixels dataset: (copy)
# n938.h5[["fixels"]] <- n941.fixels_ds  # not working.....
# n938.h5[["fixels"]] <- matrix(1:10, ncol=2)
n938.h5[["fixels"]] <- h5read(fn_fa_n941, "fixels")

# voxels dataset: (copy)
n938.h5[["voxels"]] <- h5read(fn_fa_n941, "voxels")

# scalars:
scalars.grp <- n938.h5$create_group("scalars")
FDC.grp <- scalars.grp$create_group("FDC")
FDC.grp[["ids"]] <- matrix(0:(938-1), ncol=938)   # 0:(938-1)
FDC.grp[["values"]] <- h5read(fn_fa_n941, "scalars/FDC/values")[,ids_tokeep]   # will take ~5min

n938.h5$close_all()
h5closeAll()

h5ls(fn_fa_n941)
h5ls(fn_fa_n938)

# n941.h5$close_all()
# n938.h5$close_all()



### save .h5 file with different number of subjects #####
# from 1st subject to n_th subject in ltn_FDC_n938.h5

list_num_subj <- c(50,100,200,300,500,750)  # max=938  # c(30,50,100,200,300,500,750) 
for (num_subj in list_num_subj) {
  message(paste0("number of subjects = ",toString(num_subj)))
  
  tic()
  
  fn.new <- paste0(folder_josiane, "/", "ltn_FDC_n",toString(num_subj), ".h5")
  new.h5 <- H5File$new(fn.new, mode="a")
  
  # fixels dataset: copy
  new.h5[["fixels"]] <- h5read(fn_fa_n938, "fixels")
  
  # voxels dataset: copy
  new.h5[["voxels"]] <- h5read(fn_fa_n938, "voxels")
  
  # scalars:
  scalars.grp <- new.h5$create_group("scalars")
  FDC.grp <- scalars.grp$create_group("FDC")
  FDC.grp[["ids"]] <- matrix(0:(num_subj-1), ncol=num_subj)   
  FDC.grp[["values"]] <- h5read(fn_fa_n938, "scalars/FDC/values")[,1:num_subj] 
  
  new.h5$close_all()
  
  h5ls(fn.new)
  
  toc(log=TRUE)    # <=1min for num_subj<=100; when num_subj=750: almost 3min
}



### save .csv file with different number of subjects #####
list_num_subj <- c(30,50,100,200,300,500,750)
for (num_subj in list_num_subj) {
  message(paste0("number of subjects = ",toString(num_subj)))
  
  fn.csv.new <- paste0(folder_josiane,"/","df_example_n",toString(num_subj))
  df.new <- df_example[1:num_subj,]
  
  write.csv(df.new, file = fn.csv.new, row.names = FALSE)
  
}








