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
