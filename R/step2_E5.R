step2_E5 <- function(stepPrev, data, threshold, criterion){
  r2 <- order_r2(subset(data, select = stepPrev$badvar))
  models <- list()
  fits <- list()
  num_factor <- stepPrev$num_factor+1
  for (i in 1:length(r2)){
    models[[i]] <- paste(stepPrev$goodmodelpart,
                         paste(paste0("f",num_factor), "=~",paste(r2[[i]], collapse = "+"), sep = ""),
                         sep = "\n")
    fits[[i]] <- miive(models[[i]], data, var.cov = T)
  }
  #extract the variable names that had sig sargans in each of the model fit above
  badvar_list <- lapply(fits, function(x) getbadvar(x, threshold))
  #choose the one that has less sig sargan
  #when there's multiple model with the same number of sig sargan
  #choose whoever that's first aka scaling indicator with a bigger r2
  ##added: if criterion is 1, consider both number of bad vars and r2; if 2, only consider r2. 120720.
  if(criterion=='1'){
    best_num <- which.min(lapply(badvar_list, length))
  }
  if(criterion=='2'){
    best_num <- 1
  }

  #if = 1/0 print the output
  #if problematic variable length (badvar output) > =2 create a new latent factor
  length_best <- length(badvar_list[[best_num]])
  num_badvar <- length_best
  if(num_badvar==0){
    finalobj <- list(model = models[[best_num]],
                     fit  = fits[[best_num]],
                     num_factor = num_factor,
                     num_badvar = 0)
  }
  if(!num_badvar==0){
    #save the problematic variables
    badvar <- badvar_list[[best_num]]
    #save the good variables - update the ones from the previous output
    goodvar <- stepPrev$goodvar
    #goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)
    goodvar[[num_factor]] <- setdiff(r2[[best_num]], badvar)

    # goodmodelpart <- paste(stepPrev$goodmodelpart,
    #                        paste(paste0("f", num_factor), "=~",
    #                              paste(goodvar[[num_factor]], collapse = "+"), sep = ""),
    #                        sep = "\n")
    goodmodelpart <- list(stepPrev$goodmodelpart,
                          paste(paste0("f", num_factor), "=~",
                                                        paste(goodvar[[num_factor]], collapse = "+"), sep = ""))
    finalobj <- list(model = models[[best_num]],
                     fit  = fits[[best_num]],
                     num_factor = num_factor,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart)
  }
  return(finalobj)
}
