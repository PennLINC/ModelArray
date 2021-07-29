FixelArraySeed <- function(
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

### setClass of "FixelArray" #####
setClass(
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

#' check if an object in .h5 exists
#' @param fn_h5 filename of the .h5 file
#' @param group_name full directory of this object in .h5 name
#' @param object_name name of the object, should be a string without "/"
flagObjectExistInh5 <- function(fn_h5, group_name="/results",object_name="myAnalysis") {
  h5closeAll()
  h5 <- h5ls(fn_h5)
  h5 %>% filter(.$group==group_name & .$name==object_name) %>% nrow() -> h5.nrow
  if (h5.nrow==0) {
    flagObjectExistInh5 <- FALSE
  } else {
    flagObjectExistInh5 <- TRUE
  }
}



#' check if h5 group "results" exist in current .h5 file
#' @param fn_h5 filename of the .h5 file

flagResultsGroupExistInh5 <- function(fn_h5) {
  h5closeAll()
  h5 <- h5ls(fn_h5)
  h5 %>% filter(.$group=="/" & .$name=="results") %>% nrow() -> h5.nrow
  if (h5.nrow==0) {
    flagResultsGroupExistInh5 <- FALSE
  } else {
    flagResultsGroupExistInh5 <- TRUE
  }
}

#' check if a subfolder of results exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @param analysis_name The subfolder name in "results" in .h5 file 
#' 
flagAnalysisExistInh5 <- function(fn_h5, analysis_name) {
  h5closeAll()
  h5 <- h5ls(fn_h5)
  h5 %>% filter(.$group=="/results" & .$name==analysis_name) %>% nrow() -> h5.nrow
  if (h5.nrow==0) {
    flagAnalysisExistInh5 <- FALSE
  } else {
    flagAnalysisExistInh5 <- TRUE
  }
}


#' Load fixel data output from mrtrix as an h5 file into R as a FixelArray object
#' IN THE FUTURE, THE SCALAR_TYPES AND ANALYSIS_NAMES SHOULD BE AUTOMATICALLY DETECTED!
#' @param filepath file
#' @param scalar_types expected scalars
#' @param analysis_names the subfolder names for results in .h5 file
#' @return FixelArray object
#'

FixelArray <- function(filepath, scalar_types = c("FD"), analysis_names = c("myAnalysis")) {
  ## fixel_data: 
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

    fixel_data <- t(fixel_data)   # Chenying's note: @Tinashe, you wanted to transpose voxel_data?

  }

  colnames(voxel_data) <- c("Voxel_id", "x", "y", "z")

  ids <- vector("list", length(scalar_types))

  ## scalar_data:
  scalar_data <- vector("list", length(scalar_types))

  for(x in 1:length(scalar_types)){
    
    # IT'S BETTER TO CHECK IF THIS SCALAR_TYPE EXISTS OR NOT..... - Chenying
    
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
        
        colnames(results_data[[x]]$results_matrix) <- as.character(realize(names_results_matrix))    # designate the column names
        

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

# FixelMatrix <- function(fa){
#
#
# }




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

setMethod("show", "FixelArray", function(object) {  # , group_name_results="results"
  
  # # check if there is a group of results:
  # flag_results_exist <- flagResultsExist(object, group_name_results)
  # if (flag_results_exist==TRUE) {
  #   str_results <- paste0("There is ", group_name_results, " in this FixelArray")
  # } else {
  #   str_results <- paste0("There is no ", group_name_results, " in this FixelArray")
  # }

  cat(is(object)[[1]], " located at ", object@path, "\n\n",
      format("  Fixel data:", justify = "left", width = 20), dim(fixels(object))[1], " fixels\n",
      format("  Voxel data:", justify = "left", width = 20), dim(voxels(object))[1], " voxels\n",
      format("  Subjects:", justify = "left", width = 20), dim(subjects(object)[[1]])[1], "\n",
      format("  Scalars:", justify = "left", width = 20), paste0(names(scalars(object)), collapse = ", "), "\n",
      # format("  Results:", justify = "left", width = 20), str_results, "\n",
      format("  Analyses:", justify = "left", width = 20), paste0(names(results(object)), collapse = ", "), "\n",
      sep = ""
  
  )
})


setGeneric("fixels", function(x) standardGeneric("fixels"))
setMethod("fixels", "FixelArray", function(x) x@fixels)

setGeneric("voxels", function(x) standardGeneric("voxels"))
setMethod("voxels", "FixelArray", function(x) x@voxels)

setGeneric("subjects", function(x) standardGeneric("subjects"))
setMethod("subjects", "FixelArray", function(x) x@subjects)

setGeneric("scalars", function(x, ...) standardGeneric("scalars"))
setMethod(
  "scalars",
  "FixelArray",
  function(x, ...) {

    dots <- list(...)

    if(length(dots) == 1) {

      scalar <- dots[[1]]
      x@scalars[[scalar]]

    } else {

      x@scalars

    }
  }
)

setGeneric("results", function(x, ...) standardGeneric("results"))
setMethod(
  "results", "FixelArray", function(x, ...) {
    dots <- list(...)
    
    if(length(dots) == 1) {
      
      analysis_name <- dots[[1]]
      x@results[[analysis_name]]
      
    } else {
      
      x@results
      # message: if the type of $results_matrix is character, this may not be the case for all columns, but for columns that involve look-up table (lut)
    }
  }
)
