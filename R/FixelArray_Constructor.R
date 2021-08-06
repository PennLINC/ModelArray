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
      filepath, name = name, type = type)
    
    seed
    
  } else {
    
    stop("Improperly formatted Fixel data")
    
  }
  
}


#' Load fixel data output from mrtrix as an h5 file into R as a FixelArray object
#' TODO: IN THE FUTURE, THE SCALAR_TYPES AND ANALYSIS_NAMES SHOULD BE AUTOMATICALLY DETECTED!
#' @param filepath file
#' @param scalar_types expected scalars
#' @param analysis_names the subfolder names for results in .h5 file
#' @return FixelArray object
#' @export
#' @import methods
#'

FixelArray <- function(filepath, scalar_types = c("FD"), analysis_names = c("myAnalysis")) {
  ## fixel_data: 
  
  # TODO: try and use hdf5r instead of rhdf5 and delayedarray here
  fixel_data <- FixelArraySeed(filepath, name = "fixels", type = NA) %>%
    DelayedArray::DelayedArray()

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
        names_results_matrix <- FixelArraySeed(filepath, name = sprintf("results/%s/has_names", analysis_name), type = NA) %>%
          DelayedArray::DelayedArray()
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
      h5attr(results.analysis.grp[["results_matrix"]], "colnames") <- colnames(onemodel.onerow)   # TODO: confirm with Matt that this is fine for ConFixel
      
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

#' Analyse (fit statistical model) and write the outputs for 1 fixel
#'
#' @param formula Formula (passed to `lm()`)
#' @param fa FixelArray class
#' @param phenotypes The cohort matrix with covariates to be added to the model  
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param i_fixel The i_th fixel, starting from 1, integer. For initiating (flag_initiate = TRUE), use i_fixel=1
#' @param var.terms The list of variables to save for terms (got from lm %>% tidy())
#' @param var.model The list of variables to save for the model (got from lm %>% glance())
#' @param flag_initiate Whether this is to initiate the new group (TRUE or FALSE) - if this is the first i_fixel, then TRUE and it will return column names.

#' 
#' @param verbose Print progress messages
#' @import hdf5r

analyseOneFixel.lm <- function(i_fixel, 
                               formula, fa, phenotypes, scalar, 
                               var.terms, var.model, 
                               flag_initiate = FALSE, 
                               verbose = TRUE, ...) {
  values <- scalars(fa)[[scalar]][i_fixel,]
  dat <- phenotypes
  dat[[scalar]] <- values
  onemodel <- stats::lm(formula, data = dat, ...)   
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
  if (length(var.terms.remove) != 0) {    # if length=0, it's list(), nothing to remove
    onemodel.tidy <- select(onemodel.tidy, -all_of(var.terms.remove))
  }
  if (length(var.model.remove) != 0) {
    onemodel.glance <- select(onemodel.glance, -all_of(var.model.remove))
  }
  
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
  
  
  if (flag_initiate == TRUE) { # return the column names:
    
    # return:
    column_names = colnames(onemodel.onerow)
    column_names
    
  } else if (flag_initiate == FALSE) {  # return the one row results:

    # return: 
    onerow <- as.numeric(onemodel.onerow)   # change from tibble to numeric to save a lot of space!
    onerow
  }
  
  
}  
  
  
  
#' Write outputs from fixel-based analysis out to the h5 file. Write one results (i.e. for one analysis) at a time.
#' 
#' @param fa FixelArray object
#' @param data A data.frame object with model results at each fixel
#' @param analysis_name The subfolder name in results, holding the analysis results 
#' @param flag_overwrite If same analysis_name exists, whether overwrite or not
#'

writeResults <- function(fa, data, analysis_name = "myAnalysis", flag_overwrite=TRUE){ 

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

