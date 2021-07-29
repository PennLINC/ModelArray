### Methods of "FixelArray" ###

#' @export
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
