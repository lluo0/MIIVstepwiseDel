step1_E5 <- function(data, sigLevel, scalingCrit = 'order', correlatedErrors  = NULL){
  scalingindicator <- select_scalingind(data, sigLevel, scalingCrit)
  order_scalingind <- which(colnames(data)==scalingindicator)
  model <- paste0('f1=~', paste0(colnames(data)[order_scalingind]), '+',
                  paste0(colnames(data)[-order_scalingind], collapse = '+'))

  #add in provided list of correlated errors
  if(!is.null(correlatedErrors)){
    model <- paste0(model, '\n', correlatedErrors)
  }
  fit <- miive(model, data, var.cov = T)
  badvar <- getbadvar(fit, sigLevel)
  num_badvar <- length(badvar)



  # r2 <- order_r2(data)
  # models <- list()
  # fits <- list()
  # for (i in 1:length(r2)){
  #   models[[i]] <- paste("f1=~", paste(r2[[i]], collapse = "+"), sep = "")
  #   fits[[i]] <- miive(models[[i]], data, var.cov = T)
  # }
  # #extract the variable names that had sig sargans in each of the model fit above
  # badvar_list <- lapply(fits, function(x) getbadvar(x, sigLevel))
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
    #reorder the good variables so the scaling indicator appears first
    goodvar[[1]] <- goodvar[[1]][c(match(scalingindicator, goodvar[[1]]), setdiff(order(goodvar[[1]]),match(scalingindicator, goodvar[[1]])))]
    #save the part of the model that is good and should be untouched as the first factor for the next step
    goodmodelpart <- list()
    goodmodelpart[[1]] <- paste("f1=~", paste(goodvar[[1]], collapse = "+"), sep = "")
    ##this will become a list once we have multiple factors, so does the goodvar list.
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = 1,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart,
                     nextstep = ifelse(length(badvar!=0), 'yes', 'no'),
                     correlatedErrors = correlatedErrors)
  }
  return(finalobj)
}

step2_E5 <- function(stepPrev, data, sigLevel, scalingCrit){

  badvar <- stepPrev$badvar
  goodmodelpart <- stepPrev$goodmodelpart
  # num_factor <- stepPrev$num_factor+1
  num_factor <- stepPrev$num_factor
  correlatedErrors <- stepPrev$correlatedErrors

  # scalingindicator <- select_scalingind_stepN(data, sigLevel, scalingCrit, goodmodelpart, badvar, num_factor)
  scalingindicator <- select_scalingind_stepN(data, sigLevel, scalingCrit, stepPrev)

  order_scalingind <- which(badvar==scalingindicator)

  model <- paste(paste0(goodmodelpart, collapse = '\n'),
                 paste(paste0("f",num_factor+1), "=~",paste0(badvar[order_scalingind]), '+',
                       paste0(badvar[-order_scalingind], collapse = "+"), sep = ""),
                 sep = "\n")
  #add in provided list of correlated errors
  if(!is.null(correlatedErrors)){
    model <- paste0(model, '\n', correlatedErrors)
  }
  fit <- miive(model, data, var.cov = T)
  badvar <- getbadvar(fit, sigLevel)
  num_badvar <- length(badvar)

  if(num_badvar==0){
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = num_factor+1,
                     num_badvar = 0,
                     nextstep = 'no',
                     correlatedErrors = correlatedErrors)
  }
  if(num_badvar!=0){
    # #save the problematic variables
    # badvar <- getbadvar(fit, sigLevel)
    #save the good variables - update the ones from the previous output
    goodvar <- stepPrev$goodvar
    #goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)
    goodvar[[num_factor+1]] <- setdiff(stepPrev$badvar, badvar)

    #reorder the goodvar part again so the new scaling indicator appears first on the new list
    goodvar[[num_factor+1]] <- goodvar[[num_factor+1]][c(match(scalingindicator, goodvar[[num_factor+1]]),
                                                         setdiff(order(goodvar[[num_factor+1]]),match(scalingindicator, goodvar[[num_factor+1]])))]

    # goodmodelpart <- paste(stepPrev$goodmodelpart,
    #                        paste(paste0("f", num_factor), "=~",
    #                              paste(goodvar[[num_factor]], collapse = "+"), sep = ""),
    #                        sep = "\n")

    # goodmodelpart <- list(stepPrev$goodmodelpart,
    #                       paste(paste0("f", num_factor+1), "=~",
    #                             paste(goodvar[[num_factor+1]], collapse = "+"), sep = ""))
    goodmodelpart[[num_factor+1]] <- paste(paste0("f", num_factor+1), "=~",
                                           paste(goodvar[[num_factor+1]], collapse = "+"), sep = "")
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = stepPrev$num_factor+1,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart,
                     nextstep = ifelse(length(badvar!=0), 'yes', 'no'),
                     correlatedErrors = correlatedErrors)
  }
  return(finalobj)
}


stepN_E5 <- function(stepPrev, data, sigLevel, scalingCrit){

  correlatedErrors <- stepPrev$correlatedErrors

  ##first crossload the bad variables
  crossloadmodel <- lapply(stepPrev$goodmodelpart, function(x)
    paste0(x, '+',paste0(stepPrev$badvar, collapse = '+')))

  #then add in list of correlated errors
  crossloadmodel <- lapply(crossloadmodel, function(x)
    paste0(x, '\n', correlatedErrors))

  crossloadfit <- miive(paste0(crossloadmodel, collapse = '\n'), data, var.cov = T)
  ##then see if any variables actually crossload
  ##update the bad variable list
  newbadvar <- getbadvar_crossload(crossloadfit, sigLevel, stepPrev$num_factor, stepPrev$badvar)


  # #check if crossloads doesn't work at all
  # crossloadcheck <- mapply(function(x) all(x == stepPrev$badvar), newbadvar, SIMPLIFY = T)
  # crossloadcheck_TF <- all(crossloadcheck == T)
  # #if crossloadcheck_TF == T, means


  ##add the crossloaded variables to the model
  # newgoodvar_addon <- lapply(newbadvar, function(x) setdiff(x, stepPrev$badvar))
  newgoodvar_addon <- list()
  for(p in 1:length(newbadvar)){
    if(length(newbadvar[[p]])==0){
      newgoodvar_addon[[p]] <- stepPrev$badvar
    }
    # if(identical(newbadvar[[p]], stepPrev$badvar) | length(newbadvar[[p]]) == length(stepPrev$badvar)){
    #   newgoodvar_addon[[p]] <- NULL
    # }
    # if(!length(newbadvar[[p]])==0 & !identical(newbadvar[[p]], stepPrev$badvar) ){
    #   newgoodvar_addon[[p]] <- setdiff(stepPrev$badvar,newbadvar[[p]])
    # }
    if(length(newbadvar[[p]])!=0){
      newgoodvar_addon[[p]] <- setdiff(stepPrev$badvar,newbadvar[[p]])
    }
  }




  # #create new badvar
  # newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))

  newgoodmodelpart <- stepPrev$goodmodelpart
  for(p in 1:length(newgoodvar_addon)){
    if(length(newgoodvar_addon[[p]]!=0)){
      newgoodmodelpart[[p]] <- paste0(newgoodmodelpart[[p]], '+',
                                      paste0(newgoodvar_addon[[p]], collapse = '+'))
    }
  }
  newgoodvar <- mapply(c, stepPrev$goodvar, newgoodvar_addon, SIMPLIFY  = F)

  #check if need a new factor
  newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))


  if(length(unique(unlist(newbadvar)))== 0){
    if(length(setdiff(newgoodvar[[stepPrev$num_factor]],stepPrev$goodvar[[stepPrev$num_factor]])) == 0){
      model <- paste0(crossloadmodel, collapse = '\n')
      #add in provided list of correlated errors
      if(!is.null(correlatedErrors)){
        model <- paste0(model, '\n', correlatedErrors)
      }
    }
    if(length(setdiff(newgoodvar[[stepPrev$num_factor]],stepPrev$goodvar[[stepPrev$num_factor]])) != 0){
      model <- paste0(newgoodmodelpart, collapse = '\n')
      #add in provided list of correlated errors
      if(!is.null(correlatedErrors)){
        model <- paste0(model, '\n', correlatedErrors)
      }
    }
    fit <- miive(model, data, var.cov = T)
    badvar <- getbadvar(fit, sigLevel)
    num_badvar <- length(badvar)

    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
                     num_badvar = num_badvar,
                     #goodvar = newgoodvar,
                     badvar = badvar,
                     #goodmodelpart = newgoodmodelpart,
                     nextstep = 'no',
                     correlatedErrors = correlatedErrors)
  }
  if(length(unique(unlist(newbadvar)))== 1){
    if(all(unlist(newgoodmodelpart) == unlist(stepPrev$goodmodelpart))){ #then we keep everything from the last step
      finalobj <- stepPrev
      finalobj$nextstep <- 'no' #need to stop the function from running
    }else{
      #we keep the single problematic variable on the latest created factor.
      newgoodmodelpart[[length(newgoodmodelpart)]] <- paste0(newgoodmodelpart[[length(newgoodmodelpart)]], '+', newbadvar)
      model <- paste0(newgoodmodelpart, collapse = '\n')
      #add in provided list of correlated errors
      if(!is.null(correlatedErrors)){
        model <- paste0(model, '\n', correlatedErrors)
      }
      fit <- miive(model, data, var.cov = T)
      badvar <- getbadvar(fit, sigLevel)
      num_badvar <- length(badvar)

      finalobj <- list(model = model,
                       fit  = fit,
                       num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
                       num_badvar = num_badvar,
                       #goodvar = newgoodvar,
                       badvar = badvar,
                       #goodmodelpart = newgoodmodelpart,
                       nextstep = 'no',
                       correlatedErrors = correlatedErrors)
    }

  }

  if(length(unique(unlist(newbadvar))) >1){
    # model <- paste0(paste0(newgoodmodelpart, collapse = '\n'), '\n',
    #                 paste0('f',stepPrev$num_factor+1, '=~' ,paste0(newbadvar, collapse = '+')))
    #
    # #model <- paste0(newgoodmodelpart, collapse = '\n')
    # fit <- miive(model, data, var.cov = T)
    # badvar <- getbadvar(fit, sigLevel)
    # num_badvar <- length(badvar)

    # newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))
    # if(length(newbadvar) <= 1){
    #   model <- paste0(newgoodmodelpart, collapse = '\n')
    #   fit <- miive(model, data, var.cov = T)
    #   badvar <- getbadvar(fit, sigLevel)
    #   num_badvar <- length(badvar)
    #
    #   finalobj <- list(model = model,
    #                    fit  = fit,
    #                    num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
    #                    num_badvar = num_badvar,
    #                    #goodvar = newgoodvar,
    #                    badvar = badvar,
    #                    #goodmodelpart = newgoodmodelpart,
    #                    nextstep = 'no')
    # }
    #if(length(newbadvar) > 1){
    stepPrev$badvar <- newbadvar
    stepPrev$goodmodelpart <- newgoodmodelpart

    scalingindicator <- select_scalingind_stepN(data, sigLevel, scalingCrit, stepPrev)
    order_scalingind <- which(newbadvar==scalingindicator)

    model <- paste(paste0(newgoodmodelpart, collapse = '\n'),
                   paste(paste0("f",stepPrev$num_factor+1), "=~",paste0(newbadvar[order_scalingind]), '+',
                         paste0(newbadvar[-order_scalingind], collapse = "+"), sep = ""),
                   sep = "\n")
    #add in provided list of correlated errors
    if(!is.null(correlatedErrors)){
      model <- paste0(model, '\n', correlatedErrors)
    }
    fit <- miive(model, data, var.cov = T)
    badvar <- getbadvar(fit, sigLevel)
    num_badvar <- length(badvar)

    num_factor <- stepPrev$num_factor
    #update goodvar and goodmodelpart
    newgoodvar[[num_factor+1]] <- setdiff(newbadvar, badvar)

    newgoodvar[[num_factor+1]] <- newgoodvar[[num_factor+1]][c(match(scalingindicator, newgoodvar[[num_factor+1]]),
                                                               setdiff(order(newgoodvar[[num_factor+1]]),match(scalingindicator, newgoodvar[[num_factor+1]])))]


    newgoodmodelpart[[num_factor+1]] <- paste(paste0("f", num_factor+1), "=~",
                                              paste(newgoodvar[[num_factor+1]], collapse = "+"), sep = "")

    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = stepPrev$num_factor+1,
                     num_badvar = num_badvar,
                     goodvar = newgoodvar,
                     badvar = badvar,
                     goodmodelpart = newgoodmodelpart,
                     nextstep = ifelse(length(badvar!=0), 'yes', 'no'),
                     correlatedErrors = correlatedErrors)
    # }


  }




  # if(length(newbadvar)<=1){ #we do not need to create a new factor
  # if(num_badvar <= 1){
  #aka after crossloading we either no longer have problematic variables or have only one.

  # #add badvar back to the whole model, if any
  # if(length(newbadvar==1)){
  #
  #   #create new goodvar and goodmodelvar part based on new crossload addition. 060821
  #
  #   newgoodvar <- mapply(c, stepPrev$goodvar, newgoodvar_addon, SIMPLIFY  = F)
  #   newgoodmodelpart <- mapply(function(x,y) paste0(x, '+', y), stepPrev$goodmodelpart,
  #                              lapply(newgoodvar_addon, function(i) paste0(i, collapse = '+')), SIMPLIFY=FALSE)
  #   # model <- paste0(paste0(newgoodmodelpart[[-length(newgoodmodelpart)]], collapse = '\n'),
  #   #                 paste0(newgoodmodelpart[[length(newgoodmodelpart)]], "+", newbadvar))
  #   model <- paste0(newgoodmodelpart, collapse = '\n')
  # }
  # if(length(newbadvar==0)){
  #   model <- paste0(newgoodmodelpart, collapse = '\n')
  # }
  # model <- paste0(crossloadmodel, collapse = '\n')
  # fit <- miive(model, data, var.cov = T)
  # badvar <- getbadvar(fit, sigLevel)
  # num_badvar <- length(badvar)
  #
  #     finalobj <- list(model = model,
  #                      fit  = fit,
  #                      num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
  #                      num_badvar = num_badvar,
  #                      #goodvar = newgoodvar,
  #                      badvar = badvar,
  #                      #goodmodelpart = newgoodmodelpart,
  #                      nextstep = 'no')
  #   }

  # if(length(newbadvar)>1){ #aka we still need to create a new factor
  # if(num_badvar > 1){

  #create new goodvar and goodmodelvar part based on new crossload addition. 060821

  # newgoodmodelpart <- mapply(function(x,y) paste0(x, '+', y), stepPrev$goodmodelpart,
  #        lapply(newgoodvar_addon, function(i) paste0(i, collapse = '+')), SIMPLIFY=FALSE)
  # ^the above code will be problematic when we have a NULL for any element in the list of newgoodvar_addon

  # newgoodmodelpart <- stepPrev$goodmodelpart
  # for(p in 1:length(newgoodvar_addon)){
  #   if(length(newgoodvar_addon[[p]]!=0)){
  #     newgoodmodelpart[[p]] <- paste0(newgoodmodelpart[[p]], '+',
  #                                     paste0(newgoodvar_addon[[p]], collapse = '+'))
  #   }
  # }
  # newgoodvar <- mapply(c, stepPrev$goodvar, newgoodvar_addon, SIMPLIFY  = F)

  #update stepPrev for scaling indicator selection
  #   stepPrev$badvar <- newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))
  #   stepPrev$goodmodelpart <- newgoodmodelpart
  #
  #   scalingindicator <- select_scalingind_stepN(data, sigLevel, scalingCrit, stepPrev)
  #   order_scalingind <- which(newbadvar==scalingindicator)
  #
  #   model <- paste(paste0(newgoodmodelpart, collapse = '\n'),
  #                  paste(paste0("f",stepPrev$num_factor+1), "=~",paste0(newbadvar[order_scalingind]), '+',
  #                        paste0(newbadvar[-order_scalingind], collapse = "+"), sep = ""),
  #                  sep = "\n")
  #   fit <- miive(model, data, var.cov = T)
  #   badvar <- getbadvar(fit, sigLevel)
  #   num_badvar <- length(badvar)
  #
  #   num_factor <- stepPrev$num_factor
  #   #update goodvar and goodmodelpart
  #   newgoodvar[[num_factor+1]] <- setdiff(newbadvar, badvar)
  #
  #   newgoodvar[[num_factor+1]] <- newgoodvar[[num_factor+1]][c(match(scalingindicator, newgoodvar[[num_factor+1]]),
  #                                                        setdiff(order(newgoodvar[[num_factor+1]]),match(scalingindicator, newgoodvar[[num_factor+1]])))]
  #
  #
  #   newgoodmodelpart[[num_factor+1]] <- paste(paste0("f", num_factor+1), "=~",
  #                                             paste(newgoodvar[[num_factor+1]], collapse = "+"), sep = "")
  #
  #   finalobj <- list(model = model,
  #                    fit  = fit,
  #                    num_factor = stepPrev$num_factor+1,
  #                    num_badvar = num_badvar,
  #                    goodvar = newgoodvar,
  #                    badvar = badvar,
  #                    goodmodelpart = newgoodmodelpart,
  #                    nextstep = ifelse(length(badvar!=0), 'yes', 'no'))
  # }


  #
  #
  #   # #decide if needs to create a new factor
  #   # newgoodmodelpart <- stepPrev$goodmodelpart
  #   #
  #   #if the length if zero, means no crossloaded variable exists.
  #   newfactor_decision <- vector()
  #   if(length(newgoodvar_addon)==0){
  #     newfactor_decision <- 'yes'
  #   }else{
  #     ##add the crossloaded variables to the goodmodelpart syntax
  #     for(p in 1:length(newgoodvar_addon)){
  #       if(length(na.omit(newgoodvar_addon[[p]]))!=0){
  #         newgoodmodelpart[[p]] <- paste0(newgoodmodelpart[[p]], '+',
  #                                         paste0(newgoodvar_addon[[p]], collapse='+'))
  #       }
  #     }
  #     newgoodvar_vector <- na.omit(unique(unlist(newgoodvar_addon)))
  #     if(length(stepPrev$badvar)-length(newgoodvar_vector)<=1){
  #       newfactor_decision <- 'no'
  #     }else{newfactor_decision <- 'yes'}
  #   }
  #   ##then combine the badvar for each factor together
  #   newbadvar <- unique(unlist(newbadvar))
  #   if(newfactor_decision == 'no'){
  #     #then fit again
  #     newmodel <- paste0(newgoodmodelpart, collapse='\n')
  #     newfit <- miive(newmodel, data, var.cov = T)
  #     newbadvar <- getbadvar_crossload(newfit, sigLevel, stepPrev$num_factor, stepPrev$badvar)
  #     newbadvar <- unique(unlist(newbadvar))
  #
  #     if(length(newbadvar)==0){
  #       finalobj <- list(model = newmodel,
  #                        fit  = newfit,
  #                        num_factor = stepPrev$num_factor,
  #                        num_badvar = 0,
  #                        nextstep = 'no')
  #     }else{
  #       if(length(newbadvar)==1 & !newbadvar == stepPrev$badvar){
  #         newmodel <- list()
  #         for(p in 1:length(newgoodmodelpart)){
  #           newmodel[[p]] <- paste0(newgoodmodelpart[[p]], '+',
  #                                   newbadvar)
  #         }
  #         newmodel <- paste0(newmodel, collapse = '\n')
  #         newfit <- miive(newmodel, data, var.cov = T)
  #         finalobj <- list(model = newmodel,
  #                          fit  = newfit,
  #                          num_factor = stepPrev$num_factor,
  #                          num_badvar = 1,
  #                          nextstep = 'no')
  #       }
  #
  #       if(newbadvar == stepPrev$badvar){
  #         finalobj <- c(stepPrev[1:4],
  #                       nextstep = 'no')
  #       }
  #     }
  #
  #   }
  #
  #
  #   if(newfactor_decision == 'yes'){
  #     # ##then same as before
  #     # r2 <- order_r2(subset(data, select = newbadvar))
  #     # models <- list()
  #     # fits <- list()
  #     # num_factor <- stepPrev$num_factor+1
  #     # for (i in 1:length(r2)){
  #     #   models[[i]] <- paste(paste0(newgoodmodelpart, collapse = '\n'),
  #     #                        paste(paste0("f",num_factor), "=~",paste(r2[[i]], collapse = "+"), sep = ""),
  #     #                        sep = "\n")
  #     #   fits[[i]] <- miive(models[[i]], data, var.cov = T)
  #     # }
  #     # #extract the variable names that had sig sargans in each of the model fit above
  #     # badvar_list <- lapply(fits, function(x) getbadvar(x, sigLevel))
  #     # #choose the one that has less sig sargan
  #     # #when there's multiple model with the same number of sig sargan
  #     # #choose whoever that's first aka scaling indicator with a bigger r2
  #     # ##added: if criterion is 1, consider both number of bad vars and r2; if 2, only consider r2. 120720.
  #     # if(criterion=='1'){
  #     #   best_num <- which.min(lapply(badvar_list, length))
  #     # }
  #     # if(criterion=='2'){
  #     #   best_num <- 1
  #     # }
  #     #
  #     # #if = 1/0 print the output
  #     # #if problematic variable length (badvar output) > =2 create a new latent factor
  #     # length_best <- length(badvar_list[[best_num]])
  #     # num_badvar <- length_best
  #     num_factor <- stepPrev$num_factor+1
  #     scalingindicator <- select_scalingind_stepN(data, sigLevel, scalingCrit, newgoodmodelpart, newbadvar, num_factor)
  #
  #
  #     order_scalingind <- which(newbadvar==scalingindicator)
  #
  #     model <- paste(paste0(newgoodmodelpart, collapse = '\n'),
  #                    paste(paste0("f",num_factor), "=~",paste0(newbadvar[order_scalingind]), '+',
  #                          paste0(newbadvar[-order_scalingind], collapse = "+"), sep = ""),
  #                    sep = "\n")
  #     fit <- miive(model, data, var.cov = T)
  #     badvar <- getbadvar(fit, sigLevel)
  #     num_badvar <- length(badvar)
  #
  #     if(num_badvar==0){
  #       finalobj <- list(model = model,
  #                        fit  = fit,
  #                        num_factor = num_factor,
  #                        num_badvar = num_badvar,
  #                        nextstep = 'no')
  #     }
  #     if(num_badvar>0){
  #
  #       #save the good variables - update the ones from the previous output
  #       goodvar <- stepPrev$goodvar
  #       #goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)
  #       goodvar[[num_factor]] <- setdiff(newbadvar, badvar)
  #
  #       goodmodelpart <- c(newgoodmodelpart,
  #                          paste(paste0("f", num_factor), "=~",
  #                                paste(goodvar[[num_factor]], collapse = "+"), sep = ""))
  #       finalobj <- list(model = model,
  #                        fit  = fit,
  #                        num_factor = num_factor,
  #                        num_badvar = num_badvar,
  #                        goodvar = goodvar,
  #                        badvar = badvar,
  #                        goodmodelpart = goodmodelpart,
  #                        nextstep = 'yes')
  #     }
  #   }

  return(finalobj)
}


EFAmiive5 <- function(data, sigLevel = .05, scalingCrit = 'order', correlatedErrors = NULL){
  step1 <- step1_E5(data, sigLevel, scalingCrit, correlatedErrors)
  # finalobj <- step1
  # if(!length(step1$badvar)==0){
  #   step2 <- step2_E5(step1, data, sigLevel, scalingCrit)
  #   finalobj <- step2
  #   if(!length(step2$badvar)==0){
  #     stepN <- stepN_E5(step2, data, sigLevel, scalingCrit)
  #     finalobj <- stepN
  #     while(stepN$nextstep == 'yes'){
  #       stepN <- stepN_E5(stepN, data, sigLevel, scalingCrit)
  #       finalobj <- stepN
  #     }
  #   }
  # }
  if(step1$nextstep == 'no'){
    finalobj <- step1
  }
  if(step1$nextstep == 'yes'){
    step2 <- step2_E5(step1, data, sigLevel, scalingCrit )
    if(step2$nextstep == 'no'){
      finalobj <- step2
    }
    if(step2$nextstep =='yes'){
      stepN <- stepN_E5(step2, data, sigLevel, scalingCrit)
      finalobj <- stepN
      while(stepN$nextstep =='yes'){
        stepN <- stepN_E5(stepN, data, sigLevel, scalingCrit)
        finalobj <- stepN
      }

    }
  }
  return(finalobj[1:4])
}


