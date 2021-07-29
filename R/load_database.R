# TODO: remove this file

# this script contains code that initialises the Fixel HDF5 backend

# the goal is to build and execute this docker command
# docker run -ti --rm --name=fixeldb -v $PWD/:/inputs pennbbl/fixeldb:latest --index-file /inputs/FD/index.mif --directions-file /inputs/FD/directions.mif --cohort-file /inputs/fd_inputs_all.csv


CreateFixelArrayFile <- function(index_file, directions_file, cohort_file, output_h5 = 'fixels.h5',
                          fixel_directory='/', img_name = "pennbbl/fixeldb", remove_img = TRUE) {

  # check if the image already exists
  docker <- stevedore::docker_client()

  img_exists <- docker_available(img_name, tag="latest", docker)

  if(!img_exists){

    message("Docker image not found locally! Pulling from Dockerhub")
    docker$image$pull(name = img_name, tag = "latest")

  }
  
  img_exists <- docker_available(img_name, tag="latest", docker)
  if(!img_exists){
    stop("Unable to pull docker image to create FixelArray Backend!")
  }
  
  # create the docker command
  
  command = glue::glue(
    "docker run",
    " {ifelse(remove_img, '--rm', '')}",
    " --name=fixelarray",
    " -v {fixel_directory}:/inputs",
    " {img_name}:latest",
    " --index-file /inputs/{index_file}",
    " --directions-file /inputs/{directions_file}",
    " --cohort-file /inputs/{cohort_file}",
    " --output-hdf5 {output_h5}",
    sep = " ")

  # run the docker command
  message("Trying the Docker image script...")
  print(command)
  out <- system(command)

  # ensure it worked
  if(out != 0){
    message("Error creating FixelArray File!")
  } else {
    message(glue::glue("FixelArray file created in fixel directory as {fixel_directory}/{output_h5}. You can now read this into R with FixelArray()"))
  }

}

docker_available <- function(name, tag, docker_client){
  
  # this function makes sure the docker image is available 
  
  search_result <- docker_client$image$list() %>%
    dplyr::as_tibble() %>%
    tidyr::unnest(repo_tags) %>%
    filter(stringr::str_detect(repo_tags, glue::glue("{name}:{tag}")))
  
  ifelse(
    nrow(search_result) > 0, 
    return(TRUE), 
    return(FALSE)
  )
  
}

# 
# # to read attributes -> rhdf5::h5readAttributes("/storage/fixel_stats_testing/fixel_components.h5", "results/has_names")
