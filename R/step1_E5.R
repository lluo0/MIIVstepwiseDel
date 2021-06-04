step1_E5 <- function(data, threshold){
  r2 <- order_r2(data)
  models <- list()
  fits <- list()
  for (i in 1:length(r2)){
    models[[i]] <- paste("f1=~", paste(r2[[i]], collapse = "+"), sep = "")
    fits[[i]] <- miive(models[[i]], data, var.cov = T)
  }
  #extract the variable names that had sig sargans in each of the model fit above
  badvar_list <- lapply(fits, function(x) getbadvar(x, threshold))
  #choose the one that has less sig sargan
  #when there's multiple model with the same number of sig sargan
  #choose whoever that's first aka scaling indicator with a bigger r2
  best_num <- which.min(lapply(badvar_list, length))
  #if = 1/0 print the output
  #if problematic variable length (badvar output) > =2 create a new latent factor
  length_best <- length(badvar_list[[best_num]])
  num_badvar <- length_best
  if(num_badvar == 0){
    finalobj <- list(model = models[[best_num]],
                     fit  = fits[[best_num]],
                     num_factor = 1,
                     num_badvar = 0)
  }
  if(!num_badvar ==0){
    #save the problematic variables
    badvar <- badvar_list[[best_num]]
    #save the good variables
    goodvar <- list()
    goodvar[[1]] <- setdiff(r2[[best_num]], badvar)
    #save the part of the model that is good and shoule be untouched for the next step
    goodmodelpart <- paste("f1=~", paste(goodvar[[1]], collapse = "+"), sep = "")
    ##this will become a list once we have multiple factors, so does the goodvar list.
    finalobj <- list(model = models[[best_num]],
                     fit  = fits[[best_num]],
                     num_factor = 1,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart)
  }
  return(finalobj)
}
