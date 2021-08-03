# testing lm's weights. 

i<-1
scalar <- "FD"
formula <- FD~age
values <- scalars(fa)[[scalar]][i,]
dat <- phenotypes
dat[[scalar]] <- values
stats::lm(formula, data = dat, weights=rep(1,50)) %>%   # works..
  broom::tidy() %>%
  dplyr::mutate(fixel_id = i-1)


# example from: https://stackoverflow.com/questions/38683076/ellipsis-trouble-passing-to-lm

LmWrapper <- function(df, fmla, ...) {
  # get names of stuff in ...
  argNames = sapply(substitute(list(...))[-1L], deparse)
  # look for identical names in df
  m = match(names(df), argNames, 0L)
  # store other arguments from ... in a list
  args = list(eval(parse(text = argNames[-m])))
  # name the list
  names(args) = names(argNames[-m])
  # store complete values in args, instead of just references to columns
  # the unlist code is rather ugly, the goal is to create a list where every
  # element is a column of interest
  args[names(argNames)[m]] = unlist(apply(df[, as.logical(m), drop = FALSE], 
                                          2, list), recursive = FALSE)
  # also put other stuff in there
  args$formula = fmla
  args$data = df
  # do lm
  est = do.call(lm, args)
  list(model = est)
}

data(airquality)

airquality$subset = airquality$Solar.R > 200
lm_airquality <- LmWrapper(airquality, Ozone ~ Wind, weights = Temp, subset = subset, 
          method = 'qr')