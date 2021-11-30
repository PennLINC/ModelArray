# Exported Functions

### setClass of "FixelArray" #####
#' An S4 class to represent a bank account.
#'
#' @slot fixels A DelayedArray object of fixel data
#' @slot voxels A DelayedArray object of voxel indeces
#' @slot results An h5 group of FixelArray analysis outputs
#' @slot subjects A list of subject labels
#' @slot scalars A list of scalars measured by the fixels
#' @slot path Path to the h5 file on disk
#' @importClassesFrom DelayedArray DelayedArray
FixelArray <- setClass(
  "FixelArray",
  #contains="DelayedArray",
  slots = c(
    fixels="DelayedArray",
    voxels="DelayedArray",
    results="list",
    subjects="list",
    scalars="list",
    path="character"
  )
)


#' FixelArraySeed
#'
#' Generates a "seed" for the h5 file format. A wrapper around HDF5ArraySeed
#' used to instantiate a delayed array
#' 
#' @param filepath Path to an existing h5 file.
#' @param name Name of the group/field in the h5 file
#'
#' @noRd
FixelArraySeed <- function(
  # TODO write a test for this: checks that the h5 file has the right fields
  
  filepath,
  name = "fixels",
  type = NA) {
  
  if(all(
    c("fixels", "voxels", "scalars")
    %in%
    rhdf5::h5ls(filepath)$name
  )
  ) {
    
    seed = HDF5Array::HDF5ArraySeed(    
      filepath, name = name, type = type)   # HDF5Array is also from BioConductor...
    
    seed
    
  } else {
    
    stop("Improperly formatted Fixel data")
    
  }
  
}


#' Load fixel data output from mrtrix as an h5 file into R as a FixelArray object
#' Tips for debugging: 
#' if you run into this error: "Error in h(simpleError(msg, call)) : error in evaluating the argument 'seed' in selecting a method for function 'DelayedArray': HDF5. Symbol table. Can't open object." Then please check if you give correct "scalar_types" - check via h5ls(filename_for_h5)
#' TODO: IN THE FUTURE, THE SCALAR_TYPES AND ANALYSIS_NAMES SHOULD BE AUTOMATICALLY DETECTED!
#' @param filepath file
#' @param scalar_types expected scalars
#' @param analysis_names the subfolder names for results in .h5 file
#' @return FixelArray object
#' @export
#' @import methods
#' @import dplyr  # for %>%

FixelArray <- function(filepath, scalar_types = c("FD"), analysis_names = c("myAnalysis")) {
  ## fixel_data: 
  
  # TODO: try and use hdf5r instead of rhdf5 and delayedarray here
  # fn.h5 <- H5File$new(filepath, mode="a")    # open; "a": creates a new file or opens an existing one for read/write
  # fixel_data <- fn.h5[["fixels"]]
  # NOTE: without DelayedArray (Bioconductor), the fixel_data won't look like a regular matrix in R or get transposed; 
  # NOTE: I also need to test if only using hdf5r can still extract scalars(fa)[["FD"]]
  
  fixel_data <- FixelArraySeed(filepath, name = "fixels", type = NA) %>%
    DelayedArray::DelayedArray()    # NOTE: without DelayedArray (BioConductor), the fixel_data won't look like a regular matrix in R or get transposed

  if(dim(fixel_data)[2] != 5) {

    fixel_data <- t(fixel_data)

  }

  colnames(fixel_data) <- c("Fixel_id", "Voxel_id", "x", "y", "z")
  
  ## voxel_data:
  voxel_data <- FixelArraySeed(filepath, name = "voxels", type = NA) %>%
    DelayedArray::DelayedArray()

  if(dim(voxel_data)[2] != 4) {

    fixel_data <- t(fixel_data)   # TODO: transpose voxel_data? leave it here for now.

  }

  colnames(voxel_data) <- c("Voxel_id", "x", "y", "z")

  ids <- vector("list", length(scalar_types))

  ## scalar_data:
  scalar_data <- vector("list", length(scalar_types))

  for(x in 1:length(scalar_types)){
    
    # TODO: IT'S BETTER TO CHECK IF THIS SCALAR_TYPE EXISTS OR NOT..... - Chenying
    
    # /scalars/<scalar_type>/values:
    scalar_data[[x]] <- FixelArraySeed(filepath, name = sprintf("scalars/%s/values", scalar_types[x]), type = NA) %>%
      DelayedArray::DelayedArray()

    if(dim(scalar_data[[x]])[1] < dim(scalar_data[[x]])[2]){
      scalar_data[[x]] <- t(scalar_data[[x]])
    }
    
    # /scalars/<scalar_type>/ids:
    ids[[x]] <- FixelArraySeed(filepath, name = sprintf("scalars/%s/ids", scalar_types[x]), type = NA) %>%
      DelayedArray::DelayedArray()

    if(dim(ids[[x]])[1] < dim(ids[[x]])[2]){
      ids[[x]] <- t(ids[[x]])
    }
  }

  names(scalar_data) <- scalar_types
  names(ids) <- scalar_types

  
  ## results:
  # first, we need to check if results group exists in this .h5 file
  flag_results_exist <- flagResultsGroupExistInh5(filepath)
  # message(flag_results_exist)
  if (flag_results_exist==FALSE) {
    results_data <- list()

  } else {     # results group exist --> to load subfolders
    results_data <- vector("list", length(analysis_names))
    
    for (x in 1:length(analysis_names)) {
      analysis_name <- analysis_names[x]

      # we need to check if this subfolder exists in this .h5 file:
      flag_analysis_exist <- flagAnalysisExistInh5(filepath, analysis_name =analysis_name)
      if (flag_analysis_exist==FALSE) {
        stop(paste0("This analysis: ",analysis_name, " does not exist..."))
      } else {    # exists
        # /results/<analysis_name>/has_names:
        names_results_matrix <- rhdf5::h5readAttributes(filepath, name = sprintf("results/%s/results_matrix", analysis_name))$colnames # after updating writeResults()
        
        # names_results_matrix <- FixelArraySeed(filepath, name = sprintf("results/%s/has_names", analysis_name), type = NA) %>%
        #   DelayedArray::DelayedArray()
        # if (dim(names_results_matrix)[1]<dim(names_results_matrix[2]){
        #   names_results_matrix <- t(names_results_matrix)
        # }
        
        # /results/<analysis_name>/results_matrix:
        results_data[[x]]$results_matrix <- FixelArraySeed(filepath, name = sprintf("results/%s/results_matrix", analysis_name), type = NA) %>%
          DelayedArray::DelayedArray()
        
        if (dim(results_data[[x]]$results_matrix)[2] != length(names_results_matrix)) {  # transpose if needed
          results_data[[x]]$results_matrix <- t(results_data[[x]]$results_matrix)
        }
        
        colnames(results_data[[x]]$results_matrix) <- as.character(DelayedArray::realize(names_results_matrix))    # designate the column names
        

        # /results/<analysis_name>/lut_col?:   # LOOP OVER # OF COL OF $RESULTS_MATRIX, AND SEE IF THERE IS LUT_COL
        for (i_col in 1:length(names_results_matrix)) {
          object_name <- paste0("lut_forcol",as.character(i_col))
          flag_lut_exist <- flagObjectExistInh5(filepath, group_name=paste0("/results/",analysis_name),object_name=object_name)
          if (flag_lut_exist == TRUE) {
            
            lut <- FixelArraySeed(filepath, name = paste0("results/", analysis_name,"/",object_name), type = NA) %>%
              DelayedArray::DelayedArray() 
            
            # results_data[[x]]$lut[[i_col]] <- lut
            
            # turn values in results_matrix into factors | HOWEVER, this also makes the entire $results_matrix into type "character"....
            lut %>% as.character() -> lut
            for (j_lut in 1:length(lut)) {
              str_lut <- lut[j_lut]
              idx_list <- results_data[[x]]$results_matrix[,i_col] %in% c(j_lut)
              results_data[[x]]$results_matrix[idx_list,i_col] <- lut[j_lut]
            }
            
          # } else {  # the lut for this column does not exist
          #   results_data[[x]]$lut[[i_col]] <- NULL
          }
          
        }

        # name the analysis:
        names(results_data)[[x]] <- analysis_name
        
        
        # NOTES:
        # if there is no "$lut", we can remove "$results_matrix", so that results(FixelArray) would look like: $<myAnalysis>, instead of $<myAnalysis>$results_matrix
        
      }
    }
  }
  
  
  
    
  
  new(
    "FixelArray",
    fixels = fixel_data,
    voxels = voxel_data,
    subjects = ids,
    scalars = scalar_data,
    results = results_data,   # ISSUE: LHS SHOULD BE THE SAME AS THE NAME IN THE H5 FILE, NOT NECESSARY CALLED "results"
    path = filepath
  )

}

#' Analyse (fit statistical model) and write the outputs for 1 fixel
#'
#' @param formula Formula (passed to `lm()`)
#' @param fa FixelArray class
#' @param phenotypes The cohort matrix with covariates to be added to the model  
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param fn.output.h5 Opened h5 file (via H5File$new(filename, mode="a"))
#' @param analysis_name The subgroup name in results, holding the analysis results 
#' @param i_fixel The i_th fixel, starting from 1, integer. For initiating (flag_initiate = TRUE), use i_fixel=1
#' @param var.terms The list of variables to save for terms (got from lm %>% tidy())
#' @param var.model The list of variables to save for the model (got from lm %>% glance())
#' @param flag_initiate Whether this is to initiate the new group (TRUE or FALSE) - if this is the first i_fixel, then TRUE.
#' @param results.grp # TODO: to explain
#' @param results.analysis.grp
#' 
#' @param overwrite Whether to overwrite or not    # TODO: to make this description more clear
#' @param verbose Print progress messages
#' @import hdf5r
#' @import broom
#' @import dplyr
 
analyseNwriteOneFixel.lm <- function(i_fixel, 
                                     formula, fa, phenotypes, scalar, fn.output.h5, analysis_name = "lm", 
                                     var.terms, var.model, 
                                     flag_initiate = FALSE, overwrite = TRUE,
                                     results.grp = NULL, results.analysis.grp = NULL, results_matrix_ds = NULL,
                                     verbose = TRUE, ...) {
  values <- scalars(fa)[[scalar]][i_fixel,]
  dat <- phenotypes
  dat[[scalar]] <- values
  onemodel <- stats::lm(formula, data = dat, ...)   
  onemodel.tidy <- onemodel %>% broom::tidy()
  onemodel.glance <- onemodel %>% broom::glance()
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
  onemodel.tidy <- dplyr::select(onemodel.tidy, -all_of(var.terms.remove))
  onemodel.glance <- dplyr::select(onemodel.glance, -all_of(var.model.remove))
  
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
  
  # now you can get the headers, # of columns, etc of the output results
  
  
  if (flag_initiate == TRUE) { # initiate the saving:
    # check if group "results" already exists!
    if (fn.output.h5$exists("results") == TRUE) { # group "results" exist
      results.grp <- fn.output.h5$open("results")
    } else {
      results.grp <- fn.output.h5$create_group("results")
    }
    
    # check if group "results\<analysis_name>" exists:
    if (results.grp$exists(analysis_name) == TRUE & overwrite == FALSE) {
      warning(paste0(analysis_name, " exists but not to overwrite!"))
      # TODO: add checker for exisiting analysis_name, esp the matrix size
      results.analysis.grp <- results.grp$open(analysis_name)
      results_matrix_ds <- results.analysis.grp[["results_matrix"]]
      
    } else {     # not exist; or exist & overwrite: to create
      if (results.grp$exists(analysis_name) == TRUE & overwrite == TRUE) {  # delete existing one first
        results.grp$link_delete(analysis_name)   # NOTE: the file size will not shrink after your deletion.. this is because of HDF5, regardless of package of hdf5r or rhdf5
      }
      
      # create:
      results.analysis.grp <- results.grp$create_group(analysis_name)  # create a subgroup called analysis_name under results.grp
      results.analysis.grp[["results_matrix"]] <- matrix(0, nrow=nrow(scalars(fa)[[scalar]]), ncol = ncol(onemodel.onerow))   # all 0
      results_matrix_ds <- results.analysis.grp[["results_matrix"]]   # name it
      # attach column names:
      h5attr(results.analysis.grp[["results_matrix"]], "colnames") <- colnames(onemodel.onerow)  
      
    }
    
    # return:
    output_list <- list(results.grp = results.grp,
                        results.analysis.grp = results.analysis.grp,
                        results_matrix_ds = results_matrix_ds)
    return(output_list)
    
  } else if (flag_initiate == FALSE) {  # to save this fixel:
    # assuming group "results", group "results\<analysis_name>", and its "results_matrix" dataset already exist and correctly created; they should be entered into this function

    
    # then flush into .h5 file: 
    results_matrix_ds[i_fixel,] <- as.numeric(onemodel.onerow)   # TODO: ask Matt if a column of fixel_id is needed?
    
    # # will return nothing (as results.* are just provided in the arguments....)
    # a = as.numeric(onemodel.onerow)
    # results_matrix_ds.afterflush <- results_matrix_ds
    # output_list <- list(a = a, 
    #                     results_matrix_ds.afterflush = results_matrix_ds.afterflush)
    # return(output_list)
  }
  
    
}

#' Analyse (fit linear model) and write the outputs for 1 fixel
#'
#' @param i_fixel The i_th fixel, starting from 1, integer. For initiating (flag_initiate = TRUE), use i_fixel=1
#' @param formula Formula (passed to `lm()`)
#' @param fa FixelArray class
#' @param phenotypes The cohort matrix with covariates to be added to the model  
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param var.terms The list of variables to save for terms (got from lm %>% tidy())
#' @param var.model The list of variables to save for the model (got from lm %>% glance())
#' @param flag_initiate Whether this is to initiate the new group (TRUE or FALSE) - if this is the first i_fixel, then TRUE and it will return column names.
#' 
#' @return if flag_initiate==TRUE, returns column names & list of terms of final results; if flag_initiate==FALSE, returns the final results for a fixel
#' @export
#' @import hdf5r
#' @import broom
#' @import dplyr

analyseOneFixel.lm <- function(i_fixel, 
                               formula, fa, phenotypes, scalar, 
                               var.terms, var.model, 
                               flag_initiate = FALSE, 
                               ...) {
  values <- scalars(fa)[[scalar]][i_fixel,]
  dat <- phenotypes
  dat[[scalar]] <- values
  
  # dots <- list(...)
  # dots_names <- names(dots)
  # if ("weights" %in% dots_names) {
  #   message(dots$weights)
  #   myWeights <- dots$weights
  #   dots$weights <- NULL  # remove weights from 
  #   
  #   arguments_lm <- dots
  #   
  # }
  arguments_lm <- list(...)
  arguments_lm$formula <- formula
  arguments_lm$data <- dat
  
  # onemodel <- stats::lm(formula, data = dat, ...)   
  # onemodel <- stats::lm(formula, data = dat, weights = myWeights,...)   
  onemodel <- do.call(lm, arguments_lm)   # explicitly passing arguments into lm, to avoid error of argument "weights"
  
  onemodel.tidy <- onemodel %>% broom::tidy()
  onemodel.glance <- onemodel %>% broom::glance()
  
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
  if (length(var.terms.remove) != 0) {    # if length=0, it's list(), nothing to remove
    onemodel.tidy <- dplyr::select(onemodel.tidy, -all_of(var.terms.remove))
  }
  if (length(var.model.remove) != 0) {
    onemodel.glance <- dplyr::select(onemodel.glance, -all_of(var.model.remove))
  }
  
  # adjust:
  onemodel.tidy$term[onemodel.tidy$term == "(Intercept)"] <- "Intercept"  # change the term name from "(Intercept)" to "Intercept"
  onemodel.glance <- onemodel.glance %>% mutate(term="model")   # add a column 
  
  # get the list of terms:
  list.terms <- onemodel.tidy$term
  
  # flatten .tidy results into one row:
  onemodel.tidy.onerow <- onemodel.tidy %>% tidyr::pivot_wider(names_from = term,
                                                               values_from = all_of(var.terms.orig),
                                                               names_glue = "{term}.{.value}")
  onemodel.glance.onerow <- onemodel.glance %>%  tidyr::pivot_wider(names_from = term, 
                                                                    values_from = all_of(var.model),
                                                                    names_glue = "{term}.{.value}")
  
  # combine the tables:
  onemodel.onerow <- bind_cols(onemodel.tidy.onerow, onemodel.glance.onerow)
  # add a column of fixel ids:
  colnames.temp <- colnames(onemodel.onerow)
  onemodel.onerow <- onemodel.onerow %>% tibble::add_column(fixel_id = i_fixel-1, .before = colnames.temp[1])   # add as the first column
  
  # now you can get the headers, # of columns, etc of the output results
  
  
  if (flag_initiate == TRUE) { # return the column names:
    
    # return:
    column_names = colnames(onemodel.onerow)
    toreturn <- list(column_names = column_names,
                     list.terms = list.terms)
    toreturn
    
  } else if (flag_initiate == FALSE) {  # return the one row results:

    # return: 
    onerow <- as.numeric(onemodel.onerow)   # change from tibble to numeric to save some space
    onerow
  }
  
  
}  
#' Analyse (fit mgcv::gam()) and write the outputs for 1 fixel
#'
#' @param i_fixel The i_th fixel, starting from 1, integer. For initiating (flag_initiate = TRUE), use i_fixel=1
#' @param formula Formula (passed to `mgcv::gam()`)
#' @param fa FixelArray class
#' @param phenotypes The cohort matrix with covariates to be added to the model  
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param var.smoothTerms The list of variables to save for smooth terms (got from gam %>% tidy(parametric = FALSE)). Example smooth term: age in formula "outcome ~ s(age)".
#' @param var.parametricTerms The list of variables to save for parametric terms (got from gam %>% tidy(parametric = TRUE)). Example parametric term: sex in formula "outcome ~ s(age) + sex"
#' @param var.model The list of variables to save for the model (got from gam %>% glance() and gam %>% summary())
#' @param flag_initiate Whether this is to initiate the new group (TRUE or FALSE) - if this is the first i_fixel, then TRUE and it will return column names.
#' 
#' @return if flag_initiate==TRUE, returns column names & list of terms of final results; if flag_initiate==FALSE, returns the final results for a fixel
#' @export
# #' @import hdf5r
#' @import mgcv
#' @import broom
#' @import dplyr

analyseOneFixel.gam <- function(i_fixel, formula, fa, phenotypes, scalar, 
                                var.smoothTerms, var.parametricTerms, var.model, 
                                flag_initiate = FALSE, 
                                ...) {
  values <- scalars(fa)[[scalar]][i_fixel,]
  dat <- phenotypes
  dat[[scalar]] <- values
  
  arguments <- list(...)
  arguments$formula <- formula
  arguments$data <- dat
  
  onemodel <- do.call(mgcv::gam, arguments)   # explicitly passing arguments into command, to avoid error of argument "weights"
  
  onemodel.tidy.smoothTerms <- onemodel %>% broom::tidy(parametric = FALSE)
  onemodel.tidy.parametricTerms <- onemodel %>% broom::tidy(parametric = TRUE)
  onemodel.glance <- onemodel %>% broom::glance()
  onemodel.summary <- onemodel %>% summary()
  # add additional model's stat to onemodel.glance():
  onemodel.glance[["adj.r.squared"]] <- onemodel.summary$r.sq
  onemodel.glance[["dev.expl"]] <- onemodel.summary$dev.expl

  sp.criterion.attr.name <- onemodel.summary$sp.criterion %>% attr(which = "name")
  onemodel.glance[["sp.criterion"]] <- onemodel.summary$sp.criterion[[ sp.criterion.attr.name ]]   # TODO: add this attr name as return, and write to somewhere in .h5 
  onemodel.glance[["scale"]] <- onemodel.summary$scale   # scale estimate

  num.smoothTerms <- onemodel.summary$m   # The number of smooth terms in the model.
  


  # delete columns you don't want:
  var.smoothTerms.full <- names(onemodel.tidy.smoothTerms)
  var.parametricTerms.full <- names(onemodel.tidy.parametricTerms)
  var.model.full <- names(onemodel.glance)

  # list to remove:
  var.smoothTerms.orig <- var.smoothTerms
  var.smoothTerms <- list("term", var.smoothTerms) %>% unlist()  # we will always keep "term" column
  var.smoothTerms.remove <- list()
  for (l in var.smoothTerms.full) {
    if (!(l %in% var.smoothTerms)) {
      var.smoothTerms.remove <- var.smoothTerms.remove %>% append(., l) %>% unlist()  # the order will still be kept
    }
  }

  var.parametricTerms.orig <- var.parametricTerms
  var.parametricTerms <- list("term", var.parametricTerms) %>% unlist()  # we will always keep "term" column
  var.parametricTerms.remove <- list()
  for (l in var.parametricTerms.full) {
    if (!(l %in% var.parametricTerms)) {
      var.parametricTerms.remove <- var.parametricTerms.remove %>% append(., l) %>% unlist()  # the order will still be kept
    }
  }

  var.model.remove <- list()
  for (l in var.model.full) {
    if (!(l %in% var.model)) {
      var.model.remove <- var.model.remove %>% append(., l) %>% unlist()  # the order will still be kept
    }
  }

  # remove those columns:
  if (length(var.smoothTerms.remove) != 0) {    # if length=0, it's list(), nothing to remove
    onemodel.tidy.smoothTerms <- dplyr::select(onemodel.tidy.smoothTerms, -all_of(var.smoothTerms.remove))
  }
  if (length(var.parametricTerms.remove) != 0) {    # if length=0, it's list(), nothing to remove
    onemodel.tidy.parametricTerms <- dplyr::select(onemodel.tidy.parametricTerms, -all_of(var.parametricTerms.remove))
  }
  if (length(var.model.remove) != 0) {
    onemodel.glance <- dplyr::select(onemodel.glance, -all_of(var.model.remove))
  }

  # adjust:
  if (num.smoothTerms > 0) {   # if there is any smooth term
    onemodel.tidy.smoothTerms$term[onemodel.tidy.smoothTerms$term == "(Intercept)"] <- "Intercept"  # change the term name from "(Intercept)" to "Intercept"  
  }
  if (nrow(onemodel.tidy.parametricTerms) > 0) {  # if there is any parametric term
    onemodel.tidy.parametricTerms$term[onemodel.tidy.parametricTerms$term == "(Intercept)"] <- "Intercept"  # change the term name from "(Intercept)" to "Intercept"
  }
  
    # change from s(age) to s_age: (could be s, te, etc)
  if (num.smoothTerms > 0) {   # if there is any smooth term
    for (i_row in 1:nrow(onemodel.tidy.smoothTerms)) {  # change from s(age) to s_age
      term_name <- onemodel.tidy.smoothTerms$term[i_row]
      str_list <- strsplit(term_name, split="[()]")[[1]]
      
      str <- str_list[2]   # extract string between ()
      smooth_name <- str_list[1]   # "s" or some other smooth method type such as "te"
      onemodel.tidy.smoothTerms$term[i_row] <- paste0(smooth_name, "_",str)
    }
  }
  
  
  onemodel.glance <- onemodel.glance %>% mutate(term="model")   # add a column 

  # get the list of terms:
  if (num.smoothTerms >0) {
    list.smoothTerms <- onemodel.tidy.smoothTerms$term   # if empty, gives warning
  } else {
    list.smoothTerms <- NULL
  }
  
  if (nrow(onemodel.tidy.parametricTerms)>0) {
    list.parametricTerms <- onemodel.tidy.parametricTerms$term  
  } else {
    list.parametricTerms <- NULL
  }
  

  # flatten .tidy results into one row:
  if (all(dim(onemodel.tidy.smoothTerms))) {   # not empty | if any dim is 0, all=FALSE
    onemodel.tidy.smoothTerms.onerow <- onemodel.tidy.smoothTerms %>% tidyr::pivot_wider(names_from = term,
                                                                                         values_from = all_of(var.smoothTerms.orig),
                                                                                         names_glue = "{term}.{.value}")
  } else {
    onemodel.tidy.smoothTerms.onerow <- onemodel.tidy.smoothTerms
  }
  
  if (all(dim(onemodel.tidy.parametricTerms))) {  # not empty
    onemodel.tidy.parametricTerms.onerow <- onemodel.tidy.parametricTerms %>% tidyr::pivot_wider(names_from = term,
                                                                                                 values_from = all_of(var.parametricTerms.orig),
                                                                                                 names_glue = "{term}.{.value}")
  } else {
    onemodel.tidy.parametricTerms.onerow <- onemodel.tidy.parametricTerms
  }
  
  if (all(dim(onemodel.glance))) {  # not empty
    onemodel.glance.onerow <- onemodel.glance %>%  tidyr::pivot_wider(names_from = term, 
                                                                      values_from = all_of(var.model),
                                                                      names_glue = "{term}.{.value}")
  } else {
    onemodel.glance.onerow <- onemodel.glance
  }
  

  # combine the tables:
  if ( ! all(dim(onemodel.tidy.smoothTerms.onerow)) ) {  # empty
    onemodel.onerow <- onemodel.tidy.parametricTerms.onerow
  } else {   # combine
    onemodel.onerow <- bind_cols(onemodel.tidy.smoothTerms.onerow, 
                                 onemodel.tidy.parametricTerms.onerow)
  }
  if ( ! all(dim(onemodel.onerow))   ){   # empty
    onemodel.onerow <- onemodel.glance.onerow
  } else {   # combine
    onemodel.onerow <- bind_cols(onemodel.onerow,
                                 onemodel.glance.onerow)
  }
  

  # add a column of fixel ids:
  colnames.temp <- colnames(onemodel.onerow)
  onemodel.onerow <- onemodel.onerow %>% tibble::add_column(fixel_id = i_fixel-1, .before = colnames.temp[1])   # add as the first column
  
  # now you can get the headers, # of columns, etc of the output results


  if (flag_initiate == TRUE) { # return the column names:
    
    # return:
    column_names = colnames(onemodel.onerow)
    toreturn <- list(column_names = column_names,
                     list.smoothTerms = list.smoothTerms,
                     list.parametricTerms = list.parametricTerms,
                     sp.criterion.attr.name = sp.criterion.attr.name)
    toreturn

  } else if (flag_initiate == FALSE) {  # return the one row results:

    # return: 
    onerow <- as.numeric(onemodel.onerow)   # change from tibble to numeric to save some space
    onerow
  }
}


  
#' Write outputs from fixel-based analysis out to the h5 file. Write one results (i.e. for one analysis) at a time. This is ".old": for writing results with multiple rows for one fixel
#' 
#' @param fa FixelArray object
#' @param data A data.frame object with model results at each fixel
#' @param analysis_name The subfolder name in results, holding the analysis results 
#' @param flag_overwrite If same analysis_name exists, whether overwrite or not
#'

writeResults.old <- function(fa, data, analysis_name = "myAnalysis", flag_overwrite=TRUE){ 

  # check if analysis_name subfolder already exists in the .h5 file:
  h5closeAll()
  flag_analysis_exist <- flagAnalysisExistInh5(fa@path, analysis_name)
  
  subfolder <- paste0("results/", analysis_name)
  
  if ((flag_analysis_exist == TRUE) && (flag_overwrite == TRUE)) {   # exist & want to overwrite
    rhdf5::h5delete(file=fa@path, name=subfolder)   # remove this subdirectory if it exists; otherwise, hd5f cannot be overwritten and it will throw out an error message  
  }
  
  # create results group (if needed) and create subfolder:
  flag_results_exist <- flagResultsGroupExistInh5(fa@path) # check if h5 group "results" exist
  if (flag_results_exist==FALSE) {   # results group does not exist in .h5 file
    suppressMessages(rhdf5::h5createGroup(fa@path, paste0("results"))) # first, create results group
  }
  suppressMessages(rhdf5::h5createGroup(fa@path, paste0("results/", analysis_name)))  # create subfolder, e.g. "/results/ttest"

  
  names <- names(data)  # column names
  
  ## check "data"
  if(!("data.frame" %in% class(data))) {
    stop("Results dataset is not correct; must be data of type `data.frame`")
  }
  
  
  ## check "data": make sure all columns are floats (i.e. numeric)
  for (i_col in seq(1, ncol(data), by=1)) {     # for each column of data
    col_class <- as.character(sapply(data, class)[i_col])    # class of this column
    if (col_class != "numeric") {    # the column class is not numeric
      message(paste0("the column #", as.character(i_col)," of data to save: data class is not numeric...fixing it"))
      # turn into numeric & write the notes in .h5 file...:
      factors <- data %>% pull(., var=i_col) %>% factor
      data[,i_col] <- data %>% pull(., var=i_col) %>% factor %>% as.numeric(.)    # change into numeric of 1,2,3....
      # write a LUT for this column:
      rhdf5::h5write(levels(factors), fa@path, paste0(subfolder, "/", "lut_forcol", as.character(i_col))) # save lut to .h5/results/<myAnalysis>/lut_col<?>
      
    }
  }
  

  rhdf5::h5write(t(as.matrix(data)), fa@path, paste0(subfolder, "/","results_matrix"))  # save results to .h5/results/<myAnalysis>/results_matrix

  rhdf5::h5write(names, fa@path, paste0(subfolder, "/has_names"))  # save variable names of results to .h5/results/<myAnalysis>/has_names

  message("Results file written!")
}

#' Write outputs from fixel-based analysis out to the h5 file. Write one results (i.e. for one analysis) at a time. This is ".enh": 1) change to hdf5r; 2) write results with only one row for one fixel
#' debug tip: For "Error in H5File.open(filename, mode, file_create_pl, file_access_pl)", check if there is message 'No such file or directory'. Try absolute .h5 filename.
#' 
#' @param fn.output The .h5 filename for the output, including folder directory
#' @param df.output A data.frame object with model results at each fixel, returned from FixelArray.lm() etc
#' @param analysis_name The subfolder name in results, holding the analysis results 
#' @param overwrite If same analysis_name exists, whether overwrite (TRUE) or not (FALSE)
#' @import hdf5r
#' @export

writeResults <- function(fn.output, df.output, analysis_name = "myAnalysis", overwrite=TRUE){ 
  
  # check "df.output"
  if(!("data.frame" %in% class(df.output))) {
    stop("Results dataset is not correct; must be data of type `data.frame`")
  }
  
  fn.output.h5 <- hdf5r::H5File$new(fn.output, mode="a")    # open; "a": creates a new file or opens an existing one for read/write

  # check if group "results" already exists!
  if (fn.output.h5$exists("results") == TRUE) { # group "results" exist
    results.grp <- fn.output.h5$open("results")
  } else {
    results.grp <- fn.output.h5$create_group("results")
  }
  
  # check if group "results\<analysis_name>" exists:
  if (results.grp$exists(analysis_name) == TRUE & overwrite == FALSE) {
    warning(paste0(analysis_name, " exists but not to overwrite!"))
    # TODO: add checker for exisiting analysis_name, esp the matrix size
    results.analysis.grp <- results.grp$open(analysis_name)
    results_matrix_ds <- results.analysis.grp[["results_matrix"]]
    
  } else {     # not exist; or exist & overwrite: to create
    if (results.grp$exists(analysis_name) == TRUE & overwrite == TRUE) {  # delete existing one first
      results.grp$link_delete(analysis_name)   # NOTE: the file size will not shrink after your deletion.. this is because of HDF5, regardless of package of hdf5r or rhdf5
      # TODO: add a garbage collector after saving the results
    }
    
    # create:
    results.analysis.grp <- results.grp$create_group(analysis_name)  # create a subgroup called analysis_name under results.grp
    
    # check "df.output": make sure all columns are floats (i.e. numeric)
    for (i_col in seq(1, ncol(df.output), by=1)) {     # for each column of df.output
      col_class <- as.character(sapply(df.output, class)[i_col])    # class of this column

      if ((col_class != "numeric") & (col_class != "integer")) {    # the column class is not numeric or integer
        message(paste0("the column #", as.character(i_col)," of df.output to save: data class is not numeric or integer...fixing it"))
        
        # turn into numeric & write the notes in .h5 file...:
        factors <- df.output %>% pull(., var=i_col) %>% factor
        df.output[,i_col] <- df.output %>% pull(., var=i_col) %>% factor %>% as.numeric(.)    # change into numeric of 1,2,3....
        
        # write a LUT for this column:
        results.analysis.grp[[paste0("lut_forcol", as.character(i_col))]] <- levels(factors)   # save lut to .h5/results/<myAnalysis>/lut_col<?>
        
      }
    }
    
    # save: 
    results.analysis.grp[["results_matrix"]] <- as.matrix(df.output)
    # results_matrix_ds <- results.analysis.grp[["results_matrix"]]   # name it
    
    # attach column names:
    hdf5r::h5attr(results.analysis.grp[["results_matrix"]], "colnames") <- colnames(df.output)   # NOTES: update ConFixel correspondingly
    
  }
  
  # # return:   # will not work if fn.output.h5$close_all()
  # output_list <- list(results.grp = results.grp,
  #                     results.analysis.grp = results.analysis.grp,
  #                     results_matrix_ds = results_matrix_ds)
  # return(output_list)
  
  fn.output.h5$close_all()
  

  # message("Results file written!")
  
  
}