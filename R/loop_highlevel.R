# test out our strategy of looping and computing 

### test: with one fixel, convert model to 3 tables via broom::tidy(), glance(), and augment()  #####

source("R/FixelArray_Constructor.R")
source("R/FixelArray_S4Methods.R")
source("R/utils.R")
source("R/analyse.R")
# library(FixelArray)
library(dplyr)
library(broom)
library(hdf5r)

# getwd()    # check out the current working directory

# INPUTS: #
fn <- "inst/extdata/n50_fixels.h5"
fn.output <- "/home/chenying/Desktop/fixel_project/data/data_forCircleCI_n50/n50_fixels_output.h5"    # TODO: now have to use the absolute path instead of relative... not sure why | relative:  "../data/data_forCircleCI_n50/n50_fixels_output.h5" 
file.copy(from=fn, to=fn.output, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE, recursive = TRUE) 

# h5closeAll()
fa <- FixelArray(fn)  # TODO: error with fa <- FixelArray(fn.output)  ???

fn.output.h5 <- H5File$new(fn.output, mode="a")    # open; "a": creates a new file or opens an existing one for read/write
fn.output.h5

fn_csv <- "inst/extdata/n50_cohort.csv"
phenotypes <- read.csv(fn_csv)

scalar <- "FD"
formula <- FD~age

var.terms <- c("estimate", "p.value")   # list of columns to keep
var.model <- c("r.squared", "p.value", "AIC")

analysis_name <- "lm"

fixel.subset <- 1:100

# loop starts, OR, for the first case #
i <- 1
values <- scalars(fa)[[scalar]][i,]
dat <- phenotypes
dat[[scalar]] <- values
onemodel <- stats::lm(formula, data = dat)   # TODO: add ... (additional arguments) here into lm!
onemodel.tidy <- onemodel %>% tidy()
onemodel.glance <- onemodel %>% glance()
# Augment accepts a model object and a dataset and adds information about each observation in the dataset. 
#   also accepts new data: Users may pass data to augment via either the data argument or the newdata argument. 
# onemodel.augment <- onemodel %>% augment()  

# delete columns you don't want:
var.terms.full <-names(onemodel.tidy)

var.model.full <- names(onemodel.glance)

# list to remove:
var.terms.orig <- var.terms
var.terms <- list("term", var.terms) %>% unlist()    # we will always keep "term" column
var.terms.remove <- list()   
for (l in var.terms.full) {
  if (!(l %in% var.terms)) {
    var.terms.remove <- var.terms.remove %>% append(., l) %>% unlist()  # the order will still be kept
  }
}

var.model.remove <- list()
for (l in var.model.full) {
  if (!(l %in% var.model)) {
    var.model.remove <- var.model.remove %>% append(., l) %>% unlist()  # the order will still be kept
  }
}

# remove those columns:
onemodel.tidy <- select(onemodel.tidy, -all_of(var.terms.remove))
onemodel.glance <- select(onemodel.glance, -all_of(var.model.remove))

# adjust:
onemodel.tidy$term[onemodel.tidy$term == "(Intercept)"] <- "Intercept"  # change the term name from "(Intercept)" to "Intercept"
onemodel.glance <- onemodel.glance %>% mutate(term="model")   # add a column 

# flatten .tidy results into one row:
onemodel.tidy.onerow <- onemodel.tidy %>% tidyr::pivot_wider(names_from = term,
                                                             values_from = all_of(var.terms.orig),
                                                             names_glue = "{term}.{.value}")
onemodel.glance.onerow <- onemodel.glance %>%  tidyr::pivot_wider(names_from = term, 
                                                                  values_from = all_of(var.model),
                                                                  names_glue = "{term}.{.value}")
# TODO: change the potential strings in the table into numerics + lut


# combine the tables:
onemodel.onerow <- bind_cols(onemodel.tidy.onerow, onemodel.glance.onerow)

# now you can get the headers, # of columnes, etc of the output results


# initiate the saving:
# TODO: check if group "results" already exists!
results.grp <- fn.output.h5$create_group("results")
results.analysis.grp <- results.grp$create_group(analysis_name)  # create a subgroup called analysis_name under results.grp
results.analysis.grp[["results_matrix"]] <- matrix(0, nrow=length(fixel.subset), ncol = ncol(onemodel.onerow))   # all 0
results_matrix_ds <- results.analysis.grp[["results_matrix"]]   # name it
# attach column names:
h5attr(results.analysis.grp[["results_matrix"]], "colnames") <- colnames(onemodel.onerow)   # TODO: confirm with Matt that this is fine for ConFixel

# then flush into .h5 file: 
results_matrix_ds[i,] <- as.numeric(onemodel.onerow)


# close the file
fn.output.h5$close_all()


### test using multicore to save the results #####
fn.test <- "/home/chenying/Desktop/fixel_project/data/data_forCircleCI_n50/testsaving.h5"
fn.test.h5 <- H5File$new(fn.test, mode="a")
fn.test.h5

mat <- 1:100

if (fn.test.h5$exists("results") == TRUE) {
  fn.test.h5$link_delete("results")
}

results.grp <- fn.test.h5$create_group("results")
# results.grp <- fn.test.h5$open("results")

results.grp[["mat"]] <- matrix(0, nrow=length(mat), ncol = 1) 
mat_ds <- results.grp[["mat"]]

lapply(1:100, function(i,...){# works
  mat_ds[i,] <- mat[i]
})


parallel::mclapply(1:100, function(i,...){    # does not work......
  mat_ds[i,] <- mat[i]
}, mc.cores = 2)


library(foreach)
library(doParallel)

# registerDoParallel(cores = 2)

cl <- makeCluster(2)
registerDoParallel(cl)


# getDoParWorkers()   # get how many workers foreach is going to use

tic()
foreach(i = 1:100, .combine=combine, .packages='hdf5r') %dopar%  {    # Seems foreach always throw out outputs... remove it with rm(temp) and call garbage collection gc()
  mat_ds[i,] <- mat[i]
  # return(TRUE)
}
toc()

fn.test.h5$close_all()

unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

stopCluster(cl)
registerDoSEQ()

### size of different object #####
onemodel <- stats::lm(formula, data = dat)   # raw model results is very big
object.size(onemodel)    # 26720 bytes

# if we only extract what we want (and change it into one row), it will be much smaller:
onemodel.tidy.onerow <- onemodel %>% tidy() %>% tidyr::pivot_wider(names_from = term,values_from = c(estimate, std.error,statistic, p.value), names_glue="{term}.{.value}")
onemodel.tidy.onerow
object.size(onemodel.tidy.onerow)     # 1896 bytes

# if we discard the column names and other info from tibble, and change the tibble to numeric list, the object size will be even smaller! (as the header size ~ 1KB)
onemodel.tidy.onerow.numeric <- as.numeric(onemodel.tidy.onerow)
object.size(onemodel.tidy.onerow.numeric)   # 112 bytes

# if we convert all the fixels' results into a numeric matrix, how large it will be?
a <- rnorm(1000000*10)    # 1M fixels * 10 columns to save
object.size(a)     # 80,000,048 bytes
object.size(a)/1024/1024    # 76.3MB

# turning the result matrix into data.frame (even a tibble) will only increase a little bit of size (of header, ~several KB):
# considering it's a matrix with 1M of fixels and 10 columns:
a <- matrix(rnorm(1000000*10, mean=0, sd=1), 1000000,10)
object.size(a)     # 80,000,216 bytes, around 80MB
# convert into data.frame:
a.df <- as.data.frame(a)
head(a.df)  # first several rows
colnames(a.df)  # there are already column names
class(a.df$V10)   # class is still numeric (i.e. double)
object.size(a.df)    # 80,001,904 bytes     # almost no increase
object.size(a.df) - object.size(a)     # 1688 bytes
# even tibble won't increase much size:
a.tibble <- as_tibble(a)
head(a.tibble)
dim(a.tibble)
library(pryr)
object_size(a.tibble)   # 80,002,040 B


### will a function use variable in main? #####
mya <- 1
myadd <- function(mya,b) {
  mya+b
}
myadd(4,5)  # returns 9, using the input value of mya (=4) instead of mya in the outside of function
myadd(mya=4, b=5)

myadd_wrong <- function(a,b) {  # did not define 
  mya+b
}
myadd_wrong(4,5)   # returns 1+5 = 6, as mya was not defined in the input
