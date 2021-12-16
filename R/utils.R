#' check if an object in .h5 exists
#' @param fn_h5 filename of the .h5 file
#' @param group_name full directory of this object in .h5 name
#' @param object_name name of the object, should be a string without "/"
#' @noRd
#' @importFrom rhdf5 h5ls h5closeAll
#' @importFrom dplyr filter
#' @importFrom rlang .data
flagObjectExistInh5 <- function(fn_h5, group_name="/results",object_name="myAnalysis") {
  
  rhdf5::h5closeAll()
  
  h5 <- rhdf5::h5ls(fn_h5)
  
  h5 %>% 
    dplyr::filter(.data$group==group_name & .data$name==object_name) %>% 
    nrow() -> h5.nrow
  
  if (h5.nrow==0) {
    object_exists <- FALSE
  } else {
    object_exists <- TRUE
  }
  
  object_exists
}



#' check if h5 group "results" exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @noRd
#' @importFrom rhdf5 h5ls h5closeAll
#' @importFrom dplyr filter
#' @importFrom rlang .data
flagResultsGroupExistInh5 <- function(fn_h5) {
  
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  
  h5 %>% 
    dplyr::filter(.data$group=="/" & .data$name=="results") %>% 
    nrow() -> h5.nrow
  
  if (h5.nrow==0) {
    object_exists <- FALSE
  } else {
    object_exists <- TRUE
  }
  
  object_exists
}

#' check if a subfolder of results exist in current .h5 file
#' @param fn_h5 filename of the .h5 file
#' @param analysis_name The subfolder name in "results" in .h5 file 
#' @noRd
#' @importFrom rhdf5 h5ls h5closeAll
#' @importFrom dplyr filter
#' @importFrom rlang .data
flagAnalysisExistInh5 <- function(fn_h5, analysis_name) {
  
  rhdf5::h5closeAll()
  h5 <- rhdf5::h5ls(fn_h5)
  
  h5 %>% 
    dplyr::filter(.data$group=="/results" & .data$name==analysis_name) %>% 
    nrow() -> h5.nrow
  
  if (h5.nrow==0) {
    object_exists <- FALSE
  } else {
    object_exists <- TRUE
  }
  
  object_exists
}