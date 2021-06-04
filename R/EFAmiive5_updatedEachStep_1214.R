step1_E5 <- function(data, threshold, priority = 'order'){
  scalingindicator <- select_scalingind(data, threshold, priority)
  order_scalingind <- which(colnames(data)==scalingindicator)
  model <- paste('f1=~', paste0(colnames(data)[order_scalingind]), '+',
                 paste0(colnames(data)[-order_scalingind], collapse = '+'))
  fit <- miive(model, data, var.cov = T)
  badvar <- getbadvar(fit, threshold)
  num_badvar <- length(badvar)

  # r2 <- order_r2(data)
  # models <- list()
  # fits <- list()
  # for (i in 1:length(r2)){
  #   models[[i]] <- paste("f1=~", paste(r2[[i]], collapse = "+"), sep = "")
  #   fits[[i]] <- miive(models[[i]], data, var.cov = T)
  # }
  # #extract the variable names that had sig sargans in each of the model fit above
  # badvar_list <- lapply(fits, function(x) getbadvar(x, threshold))
  # #choose the one that has less sig sargan
  # #when there's multiple model with the same number of sig sargan
  # #choose whoever that's first aka scaling indicator with a bigger r2
  # best_num <- which.min(lapply(badvar_list, length))
  # #if = 1/0 print the output
  # #if problematic variable length (badvar output) > =2 create a new latent factor
  # length_best <- length(badvar_list[[best_num]])
  # num_badvar <- length_best
  if(num_badvar <=1){
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = 1,
                     num_badvar = num_badvar)
  }
  if(!num_badvar <=1){
    goodvar <- list()
    goodvar[[1]] <- setdiff(colnames(data), badvar)
    #save the part of the model that is good and shoule be untouched for the next step
    goodmodelpart <- paste("f1=~", paste(goodvar[[1]], collapse = "+"), sep = "")
    ##this will become a list once we have multiple factors, so does the goodvar list.
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = 1,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart)
  }
  return(finalobj)
}

step2_E5 <- function(stepPrev, data, threshold, priority){
  # r2 <- order_r2(subset(data, select = stepPrev$badvar))
  # models <- list()
  # fits <- list()
  # num_factor <- stepPrev$num_factor+1
  # for (i in 1:length(r2)){
  #   models[[i]] <- paste(stepPrev$goodmodelpart,
  #                        paste(paste0("f",num_factor), "=~",paste(r2[[i]], collapse = "+"), sep = ""),
  #                        sep = "\n")
  #   fits[[i]] <- miive(models[[i]], data, var.cov = T)
  # }
  # #extract the variable names that had sig sargans in each of the model fit above
  # badvar_list <- lapply(fits, function(x) getbadvar(x, threshold))
  # #choose the one that has less sig sargan
  # #when there's multiple model with the same number of sig sargan
  # #choose whoever that's first aka scaling indicator with a bigger r2
  # ##added: if criterion is 1, consider both number of bad vars and r2; if 2, only consider r2. 120720.
  # if(criterion=='1'){
  #   best_num <- which.min(lapply(badvar_list, length))
  # }
  # if(criterion=='2'){
  #   best_num <- 1
  # }
  #
  # #if = 1/0 print the output
  # #if problematic variable length (badvar output) > =2 create a new latent factor
  # length_best <- length(badvar_list[[best_num]])
  # num_badvar <- length_best
  badvar <- stepPrev$badvar
  goodmodelpart <- stepPrev$goodmodelpart
  num_factor <- stepPrev$num_factor+1

  scalingindicator <- select_scalingind_stepN(data, threshold, priority, goodmodelpart, badvar, num_factor)

  order_scalingind <- which(badvar==scalingindicator)

  model <- paste(paste0(goodmodelpart, collapse = '\n'),
                 paste(paste0("f",num_factor), "=~",paste0(badvar[order_scalingind]), '+',
                       paste0(badvar[-order_scalingind], collapse = "+"), sep = ""),
                 sep = "\n")
  fit <- miive(model, data, var.cov = T)
  badvar <- getbadvar(fit, threshold)
  num_badvar <- length(badvar)

  if(num_badvar==0){
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = num_factor,
                     num_badvar = 0)
  }
  if(!num_badvar==0){
    #save the problematic variables
    badvar <- getbadvar(fit, threshold)
    #save the good variables - update the ones from the previous output
    goodvar <- stepPrev$goodvar
    #goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)
    goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)

    # goodmodelpart <- paste(stepPrev$goodmodelpart,
    #                        paste(paste0("f", num_factor), "=~",
    #                              paste(goodvar[[num_factor]], collapse = "+"), sep = ""),
    #                        sep = "\n")
    goodmodelpart <- list(stepPrev$goodmodelpart,
                          paste(paste0("f", num_factor), "=~",
                                paste(goodvar[[num_factor]], collapse = "+"), sep = ""))
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = num_factor,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart)
  }
  return(finalobj)
}


stepN_E5 <- function(stepPrev, data, threshold, priority){
  ##first crossload the bad variables
  crossloadmodel <- lapply(stepPrev$goodmodelpart, function(x)
    paste0(x, '+',paste0(stepPrev$badvar, collapse = '+')))
  crossloadfit <- miive(paste0(crossloadmodel, collapse = '\n'), data, var.cov = T)
  ##then see if any variables actually crossload
  ##update the bad variable list
  newbadvar <- getbadvar_crossload(crossloadfit, .05, stepPrev$num_factor, stepPrev$badvar)
  ##add the crossloaded variables to the model
  # newgoodvar_addon <- lapply(newbadvar, function(x) setdiff(x, stepPrev$badvar))
  newgoodvar_addon <- list()
  for(p in 1:length(newbadvar)){
    if(length(newbadvar[[p]])==0){
      newgoodvar_addon[[p]] <- stepPrev$badvar
    }
    if(identical(newbadvar[[p]], stepPrev$badvar)){
      newgoodvar_addon[[p]] <- NA
    }
    if(!length(newbadvar[[p]])==0 & !identical(newbadvar[[p]], stepPrev$badvar) ){
      newgoodvar_addon[[p]] <- setdiff(stepPrev$badvar,newbadvar[[p]])
    }
  }
  #decide if needs to create a new factor
  newgoodmodelpart <- stepPrev$goodmodelpart
  #if the length if zero, means no crossloaded variable exists.
  newfactor_decision <- vector()
  if(length(newgoodvar_addon)==0){
    newfactor_decision <- 'yes'
  }else{
    ##add the crossloaded variables to the goodmodelpart syntax
    for(p in 1:length(newgoodvar_addon)){
      if(length(na.omit(newgoodvar_addon[[p]]))!=0){
        newgoodmodelpart[[p]] <- paste0(newgoodmodelpart[[p]], '+',
                                        paste0(newgoodvar_addon[[p]], collapse='+'))
      }
    }
    newgoodvar_vector <- na.omit(unique(unlist(newgoodvar_addon)))
    if(length(stepPrev$badvar)-length(newgoodvar_vector)<=1){
      newfactor_decision <- 'no'
    }else{newfactor_decision <- 'yes'}
  }
  ##then combine the badvar for each factor together
  newbadvar <- unique(unlist(newbadvar))
  if(newfactor_decision == 'no'){
    #then fit again
    newmodel <- paste0(newgoodmodelpart, collapse='\n')
    newfit <- miive(newmodel, data, var.cov = T)
    newbadvar <- getbadvar_crossload(newfit, threshold, stepPrev$num_factor, stepPrev$badvar)
    newbadvar <- unique(unlist(newbadvar))

    if(length(newbadvar)==0){
      finalobj <- list(model = newmodel,
                       fit  = newfit,
                       num_factor = stepPrev$num_factor,
                       num_badvar = 0,
                       nextstep = 'no')
    }else{
      if(length(newbadvar)==1 & !newbadvar == stepPrev$badvar){
        newmodel <- list()
        for(p in 1:length(newgoodmodelpart)){
          newmodel[[p]] <- paste0(newgoodmodelpart[[p]], '+',
                                  newbadvar)
        }
        newmodel <- paste0(newmodel, collapse = '\n')
        newfit <- miive(newmodel, data, var.cov = T)
        finalobj <- list(model = newmodel,
                         fit  = newfit,
                         num_factor = stepPrev$num_factor,
                         num_badvar = 1,
                         nextstep = 'no')
      }

      if(newbadvar == stepPrev$badvar){
        finalobj <- c(stepPrev[1:4],
                      nextstep = 'no')
      }
    }

  }


  if(newfactor_decision == 'yes'){
    # ##then same as before
    # r2 <- order_r2(subset(data, select = newbadvar))
    # models <- list()
    # fits <- list()
    # num_factor <- stepPrev$num_factor+1
    # for (i in 1:length(r2)){
    #   models[[i]] <- paste(paste0(newgoodmodelpart, collapse = '\n'),
    #                        paste(paste0("f",num_factor), "=~",paste(r2[[i]], collapse = "+"), sep = ""),
    #                        sep = "\n")
    #   fits[[i]] <- miive(models[[i]], data, var.cov = T)
    # }
    # #extract the variable names that had sig sargans in each of the model fit above
    # badvar_list <- lapply(fits, function(x) getbadvar(x, threshold))
    # #choose the one that has less sig sargan
    # #when there's multiple model with the same number of sig sargan
    # #choose whoever that's first aka scaling indicator with a bigger r2
    # ##added: if criterion is 1, consider both number of bad vars and r2; if 2, only consider r2. 120720.
    # if(criterion=='1'){
    #   best_num <- which.min(lapply(badvar_list, length))
    # }
    # if(criterion=='2'){
    #   best_num <- 1
    # }
    #
    # #if = 1/0 print the output
    # #if problematic variable length (badvar output) > =2 create a new latent factor
    # length_best <- length(badvar_list[[best_num]])
    # num_badvar <- length_best
    num_factor <- stepPrev$num_factor+1
    scalingindicator <- select_scalingind_stepN(data, threshold, priority, newgoodmodelpart, newbadvar, num_factor)


    order_scalingind <- which(newbadvar==scalingindicator)

    model <- paste(paste0(newgoodmodelpart, collapse = '\n'),
                   paste(paste0("f",num_factor), "=~",paste0(newbadvar[order_scalingind]), '+',
                         paste0(newbadvar[-order_scalingind], collapse = "+"), sep = ""),
                   sep = "\n")
    fit <- miive(model, data, var.cov = T)
    badvar <- getbadvar(fit, threshold)
    num_badvar <- length(badvar)

    if(num_badvar==0){
      finalobj <- list(model = model,
                       fit  = fit,
                       num_factor = num_factor,
                       num_badvar = num_badvar,
                       nextstep = 'no')
    }
    if(num_badvar>0){

      #save the good variables - update the ones from the previous output
      goodvar <- stepPrev$goodvar
      #goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)
      goodvar[[num_factor]] <- setdiff(newbadvar, badvar)

      goodmodelpart <- c(newgoodmodelpart,
                         paste(paste0("f", num_factor), "=~",
                               paste(goodvar[[num_factor]], collapse = "+"), sep = ""))
      finalobj <- list(model = model,
                       fit  = fit,
                       num_factor = num_factor,
                       num_badvar = num_badvar,
                       goodvar = goodvar,
                       badvar = badvar,
                       goodmodelpart = goodmodelpart,
                       nextstep = 'yes')
    }
  }

  return(finalobj)
}


EFAmiive5 <- function(data, threshold = .05, priority = 'order'){
  step1 <- step1_E5(data, threshold, priority)
  finalobj <- step1
  if(!length(step1$badvar)==0){
    step2 <- step2_E5(step1, data, threshold, priority)
    finalobj <- step2
    if(!length(step2$badvar)==0){
      stepN <- stepN_E5(step2, data, threshold, priority)
      finalobj <- stepN
      while(stepN$nextstep == 'yes'){
        stepN <- stepN_E5(stepN, data, threshold, priority)
        finalobj <- stepN
      }
    }
  }
  return(finalobj[1:4])
}


