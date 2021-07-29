#' Run a linear model at each fixel location
#'
#' @param formula Formula (passed to `lm()`)
#' @param data FixelArray dataset
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param subset A vector of fixel IDs to subset
#' @param verbose Print progress messages
#' @param n_cores The number of cores to run on
#' @param pbar Print progress bar
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.lm <- function(formula, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, write = TRUE, ...){
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  # ensure we can write to fixelarray$results
  
  
  if(verbose){
    message(glue::glue("Fitting fixel-wise linear models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        lm(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, mc.cores = n_cores)
      
    } else {
      
      fits <- parallel::mclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        lm(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
        
      }, mc.cores = n_cores)
      
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        lm(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      })
      
    }
    else {
      
      fits <- lapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        lm(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      })
      
    }
  }
  
  df_out <- do.call(rbind, fits)
  
  # if(write){
  #   WriteResult(data, df_out, glue::glue("scalars/{scalar}/results"))
  # }
  
  df_out
  
}

#' Run a t.test at each fixel location
#'
#' @param formula Formula (passed to `lm()`)
#' @param data FixelArray dataset
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param subset A vector of fixel IDs to subset
#' @param verbose Print progress messages
#' @param n_cores The number of cores to run on
#' @param pbar Print progress bar
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.t.test <- function(formula, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, write = TRUE, ...){
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  # ensure we can write to fixelarray$results
  
  
  if(verbose){
    message(glue::glue("Running t.test models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        t.test(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, mc.cores = n_cores)
      
    } else {
      
      fits <- parallel::mclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        t.test(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
        
      }, mc.cores = n_cores)
      
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        t.test(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      })
      
    }
    else {
      
      fits <- lapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        t.test(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      })
      
    }
  }
  
  df_out <- do.call(rbind, fits)
  
  # if(write){
  #   WriteResult(data, df_out, glue::glue("scalars/{scalar}/results"))
  # }
  
  df_out
  
}

#' Run a GAMM4 model at each fixel location
#'
#' @param formula Formula (passed to `gamm4()`)
#' @param data FixelArray dataset
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param subset A vector of fixel IDs to subset
#' @param verbose Print progress messages
#' @param n_cores The number of cores to run on
#' @param pbar Print progress bar
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.gamm4 <- function(formula, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, write = TRUE, ...){
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  # ensure we can write to fixelarray$results
  
  
  if(verbose){
    message(glue::glue("Fitting fixel-wise linear models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        gamm4::gamm4(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, mc.cores = n_cores)
      
    } else {
      
      fits <- parallel::mclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        gamm4::gamm4(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
        
      }, mc.cores = n_cores)
      
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        gamm4::gamm4(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      })
      
    }
    else {
      
      fits <- lapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        gamm4::gamm4(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      })
      
    }
  }
  
  df_out <- do.call(rbind, fits)
  
  # if(write){
  #   WriteResult(data, df_out, glue::glue("scalars/{scalar}/results"))
  # }
  
  df_out
  
}

#' Run a GAMM4 model at each fixel location
#'
#' @param formula Formula (passed to `gamm4()`)
#' @param data FixelArray dataset
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param subset A vector of fixel IDs to subset
#' @param verbose Print progress messages
#' @param n_cores The number of cores to run on
#' @param pbar Print progress bar
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.gam <- function(formula, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, write = TRUE, ...){
# TODO: how to pass arguments to gam like in this example https://broom.tidyverse.org/reference/mgcv_tidy_gam.html  
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  # ensure we can write to fixelarray$results
  
  
  if(verbose){
    message(glue::glue("Fitting fixel-wise linear models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        mgcv::gam(formula, data = dat, ...) %>%
          broom::tidy(parametric = TRUE) %>%
          dplyr::mutate(fixel_id = i-1)
        
      }, mc.cores = n_cores)
      
    } else {
      
      fits <- parallel::mclapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        mgcv::gam(formula, data = dat, ...) %>%
          broom::tidy(parametric = TRUE) %>%
          dplyr::mutate(fixel_id = i-1)
        
        
      }, mc.cores = n_cores)
      
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        mgcv::gam(formula, data = dat, ...) %>%
          broom::tidy() %>%
          dplyr::mutate(fixel_id = i-1)
        
      })
      
    }
    else {
      
      fits <- lapply(ids, function(i, ...){
        
        values <- scalars(data)[[scalar]][i,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        mgcv::gam(formula, data = dat, ...) %>%
          broom::tidy(parametric = TRUE) %>%
          dplyr::mutate(fixel_id = i-1)
        
      })
      
    }
  }
  
  df_out <- do.call(rbind, fits)
  
  # if(write){
  #   WriteResult(data, df_out, glue::glue("scalars/{scalar}/results"))
  # }
  
  df_out
  
}

#' Run a user-defined function on each fixel location
#'
#' @param formula Formula (passed to `lm()`)
#' @param FUN User-defined modelling function; must return a 1-row named vector or data.frame
#' @param data FixelArray dataset
#' @param scalar The name of the scalar to be analysed fixel-wise
#' @param phenotypes The cohort file with covariates to be added to the model
#' @param subset A vector of fixel IDs to subset
#' @param verbose Print progress messages
#' @param n_cores The number of cores to run on
#' @param pbar Print progress bar
#' @return Tibble with the summarised model statistics at each fixel location
#' 
FixelArray.model <- function(formula, FUN, data, phenotypes, scalar, verbose = TRUE, idx = NULL, pbar = TRUE, n_cores = 1, ...){
  
  # data type assertions
  if(class(data) != "FixelArray") {
    stop("Not a fixel array for analysis")
  }
  
  # ensure we can write to fixelarray$results
  
  FUN <- match.fun(FUN)
  
  if(verbose){
    message(glue::glue("Fitting fixel-wise linear models for {scalar}", ))
  }
  
  n_models <- length(fixels(data)[,1])
  
  if(is.null(idx)){
    ids <- 1:n_models
  } else {
    ids <- idx
  }
  
  # is it a multicore process?
  if(n_cores > 1){
    
    if(pbar){
      
      fits <- pbmcapply::pbmclapply(ids, function(x, ...){
        
        values <- scalars(data)[[scalar]][x,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        FUN(formula, data = dat, ...) %>%
          broom::tidy()
        
      }, mc.cores = n_cores)
      
    } else {
      
      fits <- parallel::mclapply(ids, function(x, ...){
        
        values <- scalars(data)[[scalar]][x,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        FUN(formula, data = dat, ...) %>%
          broom::tidy()
        
        
      }, mc.cores = n_cores)
      
    }
  } else {
    
    if(pbar){
      
      fits <- pbapply::pblapply(ids, function(x, ...){
        
        values <- scalars(data)[[scalar]][x,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        FUN(formula, data = dat, ...) %>%
          broom::tidy()
        
      })
      
    }
    else {
      
      fits <- lapply(ids, function(x, ...){
        
        values <- scalars(data)[[scalar]][x,]
        dat <- phenotypes
        dat[[scalar]] <- values
        
        FUN(formula, data = dat, ...) %>%
          broom::tidy()
        
      })
      
    }
  }
  
  return(do.call(rbind, fits))
  
}