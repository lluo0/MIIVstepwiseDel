stepN_E5 <- function(stepPrev, data, sigLevel, scalingCrit){
  
  correlatedErrors <- stepPrev$correlatedErrors
  
  ##first crossload the bad variables
  crossloadmodel <- lapply(stepPrev$goodmodelpart, function(x)
    paste0(x, '+',paste0(stepPrev$badvar, collapse = '+')))
  
  # #then add in list of correlated errors
  # crossloadmodel <- lapply(crossloadmodel, function(x)
  #   paste0(x, '\n', correlatedErrors))
  # 
  if(!is.null(correlatedErrors)){
    crossloadmodel <- c(crossloadmodel, correlatedErrors)
  }
  
  crossloadfit <- miive(paste0(crossloadmodel, collapse = '\n'), data, var.cov = T)
  ##then see if any variables actually crossload
  ##update the bad variable list
  newbadvar <- getbadvar_crossload(crossloadfit, sigLevel, stepPrev$num_factor, stepPrev$badvar)
  
  #see if the crossload model would be the final model
  #aka if we have zero newbadvar
  if(length(unique(unlist(newbadvar))) == 0){ #we keep this as the final model
    # finalmodel <- list(model = crossloadmodel,
    #                    fit = crossloadfit,
    #                    num_factor = stepPrev$num_factor,
    #                    nextstep = 'no',
    #                    correlatedErrors = correlatedErrors)
    newmodel <- crossloadmodel
    newfit <- crossloadfit
    
  }
  if(length(unique(unlist(newbadvar))) == 1){ #need to further clean up the model
    #if this variable is bad on all factors, need to remove it from previous factors and only keep on the last one
    allsamebadvar <- all(mapply(identical, head(newbadvar,1), tail(newbadvar, -1)))
    
    if(allsamebadvar == F){
      newmodel <- list()
      newaddon <- lapply(newbadvar, function(x) setdiff(stepPrev$badvar, x))
      for(n in 1:length(stepPrev$goodmodelpart)){
        if(length(newaddon[[n]])==0){
          newmodel[[n]] <- stepPrev$goodmodelpart[[n]]
        }
        if(length(newaddon[[n]])!=0){
          newmodel[[n]] <- paste0(stepPrev$goodmodelpart[[n]], '+', paste0(newaddon[[n]], collapse = '+'))
        }
        
      }
      newmodel <- paste0(unlist(newmodel), collapse = '\n')
      if(!is.null(correlatedErrors)){
        newmodel <- c(newmodel, correlatedErrors)
      }
      newfit <- miive(paste0(unlist(newmodel), collapse ='\n'), data, var.cov = T)
      # finalobj <- list(model = newmodel,
      #                  fit = newfit,
      #                  nextstep = 'no')
    }
    if(allsamebadvar == T){
      otheraddon <- setdiff(stepPrev$badvar, unique(unlist(newbadvar)))
      # newmodel_firstpart <- lapply(stepPrev$goodmodelpart[-length(stepPrev$goodmodelpart)], function(x)
      #   paste0(x, '+',paste0(otheraddon, collapse = '+')))
      newmodel_firstpart <-stepPrev$goodmodelpart[-length(stepPrev$goodmodelpart)]
      newmodel_secondpart <- paste0(stepPrev$goodmodelpart[length(stepPrev$goodmodelpart)], '+',paste0(stepPrev$badvar, collapse = '+'))
      newmodel <- c(newmodel_firstpart, newmodel_secondpart)
      if(!is.null(correlatedErrors)){
        newmodel <- c(newmodel, list(correlatedErrors))
      }
      newfit <- miive(paste0(unlist(newmodel), collapse = '\n'), data, var.cov = T)
      # finalobj <- list(model = newmodel,
      #                  fit = newfit,
      #                  nextstep = 'no')
    }
  }
  if(length(unique(unlist(newbadvar))) > 1){
    newmodel <- list()
    newaddon <- lapply(newbadvar, function(x) setdiff(stepPrev$badvar, x))
    stillbadvar <- setdiff(stepPrev$badvar, unique(unlist(newaddon))) #variables that are still bad and will only be loaded on the last factor for now
    
    if(length(stillbadvar)!=0){ #then add these variables to the last factor
      newaddon[[length(newaddon)]] <- c( newaddon[[length(newaddon)]], stillbadvar)
    }
    for(n in 1:length(stepPrev$goodmodelpart)){
      if(length(newaddon[[n]])==0){
        newmodel[[n]] <- stepPrev$goodmodelpart[[n]]
      }
      if(length(newaddon[[n]])!=0){
        newmodel[[n]] <- paste0(stepPrev$goodmodelpart[[n]], '+', paste0(newaddon[[n]], collapse = '+'))
      }
      
    }
    
    # newmodel <- paste0(unlist(newmodel), collapse = '\n')
    if(!is.null(correlatedErrors)){
      newmodel <- c(newmodel, correlatedErrors)
    }
    newmodel <- paste0(unlist(newmodel), collapse = '\n')
    newfit <- miive(newmodel, data, var.cov = T)
    
    # badvar <- getbadvar(newfit, sigLevel)
  }
  
  badvar <- getbadvar(newfit, sigLevel)
  if(length(badvar)<=1){
    finalobj <- list(model = newmodel,
                     fit = newfit,
                     num_factor = stepPrev$num_factor,
                     num_badvar = length(newbadvar),
                     badvar = badvar,
                     nextstep = 'no', 
                     correlatedErrors = correlatedErrors)
  }
  if(length(badvar)>1){
    goodaddon <- lapply(newaddon, function(x) setdiff(x, badvar))
    newgoodmodelpart <- list()
    newgoodvar <- list()
    for(n in 1:length(stepPrev$goodmodelpart)){
      if(length(goodaddon[[n]])==0){
        newgoodmodelpart[[n]] <- stepPrev$goodmodelpart[[n]]
        newgoodvar[[n]] <- stepPrev$goodvar[[n]]
      }
      if(length(goodaddon[[n]])!=0){
        newgoodmodelpart[[n]] <- paste0(stepPrev$goodmodelpart[[n]], '+', paste0(goodaddon[[n]], collapse = '+'))
        newgoodvar[[n]] <- c(stepPrev$goodvar[[n]],goodaddon[[n]])
      }
      
    }
    
    stepPrev$badvar <- badvar
    stepPrev$goodmodelpart <- newgoodmodelpart
    
    scalingindicator <- select_scalingind_stepN(data, sigLevel, scalingCrit, stepPrev)
    
    order_scalingind <- which(badvar==scalingindicator)
    
    num_factor <- stepPrev$num_factor
    
    model <- paste(paste0(newgoodmodelpart, collapse = '\n'),
                   paste(paste0("f",num_factor+1), "=~",paste0(badvar[order_scalingind]), '+',
                         paste0(badvar[-order_scalingind], collapse = "+"), sep = ""),
                   sep = "\n")
    
    fit <-  miive(model, data, var.cov = T)
    newbadvar <- getbadvar(fit, sigLevel)
    
    #update finalobj
    newgoodvar <- c(newgoodvar, list(setdiff(c(badvar[order_scalingind], badvar[-order_scalingind]), newbadvar)))
    order_scalingind <- which(newgoodvar[[length(newgoodvar)]]==scalingindicator)
    if(length(newgoodvar[[length(newgoodvar)]])==1){
      newgoodmodelpart <- c(newgoodmodelpart, list(
        paste0(paste0("f",num_factor+1), "=~",paste0(newgoodvar[[length(newgoodvar)]][order_scalingind]))
      ))
    }
    if(length(newgoodvar[[length(newgoodvar)]])!=1){
      newgoodmodelpart <- c(newgoodmodelpart, list(
        paste(paste0("f",num_factor+1), "=~",paste0(newgoodvar[[length(newgoodvar)]][order_scalingind]), '+',
              paste0(newgoodvar[[length(newgoodvar)]][-order_scalingind], collapse = "+"), sep = "")
      ))
    }
  
    
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = stepPrev$num_factor+1,
                     num_badvar = length(newbadvar),
                     goodvar = newgoodvar,
                     badvar = newbadvar,
                     goodmodelpart = newgoodmodelpart,
                     nextstep = ifelse(length(newbadvar)>0, 'yes', 'no'),
                     correlatedErrors = correlatedErrors)
  }
  
  # #see if this crossload model would be the final model
  # diffnewbadvar <- setdiff(stepPrev$badvar,unique(unlist(lapply(newbadvar, function(x)
  #   setdiff(stepPrev$badvar, x)))))
  # if(length(diffnewbadvar) ==0){
  #   model <- crossloadmodel
  #   #add in provided list of correlated errors
  #   if(!is.null(correlatedErrors)){
  #     model <- paste0(model, '\n', correlatedErrors)
  #   }
  # }
  # finalobj <- list(model = model,
  #                  fit  = crossloadfit,
  #                  num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
  #                  nextstep = 'no',
  #                  correlatedErrors = correlatedErrors)
  # }
  # 
  
  # #check if crossloads doesn't work at all
  # crossloadcheck <- mapply(function(x) all(x == stepPrev$badvar), newbadvar, SIMPLIFY = T)
  # crossloadcheck_TF <- all(crossloadcheck == T)
  # #if crossloadcheck_TF == T, means
  
  # 
  # ##add the crossloaded variables to the model
  # # newgoodvar_addon <- lapply(newbadvar, function(x) setdiff(x, stepPrev$badvar))
  # newgoodvar_addon <- list()
  # for(p in 1:length(newbadvar)){
  #   if(length(newbadvar[[p]])==0){
  #     newgoodvar_addon[[p]] <- stepPrev$badvar
  #   }
  #   
  #   if(length(newbadvar[[p]])!=0){
  #     newgoodvar_addon[[p]] <- setdiff(stepPrev$badvar,newbadvar[[p]])
  #   }
  # }
  # 
  # 
  # 
  # 
  # # #create new badvar
  # # newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))
  # 
  # newgoodmodelpart <- stepPrev$goodmodelpart
  # for(p in 1:length(newgoodvar_addon)){
  #   if(length(newgoodvar_addon[[p]]!=0)){
  #     newgoodmodelpart[[p]] <- paste0(newgoodmodelpart[[p]], '+',
  #                                     paste0(newgoodvar_addon[[p]], collapse = '+'))
  #   }
  # }
  # newgoodvar <- mapply(c, stepPrev$goodvar, newgoodvar_addon, SIMPLIFY  = F)
  # 
  # #check if need a new factor
  # newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))
  # 
  # if(length(unique(unlist(newbadvar)))== 0 && length(diffnewbadvar) !=0){
  #   model <- paste0(newgoodmodelpart, collapse = '\n')
  #   #add in provided list of correlated errors
  #   if(!is.null(correlatedErrors)){
  #     model <- paste0(model, '\n', correlatedErrors)
  #   }
  # 
  # fit <- miive(model, data, var.cov = T)
  # badvar <- getbadvar(fit, sigLevel)
  # num_badvar <- length(badvar)
  # 
  # finalobj <- list(model = model,
  #                  fit  = fit,
  #                  num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
  #                  num_badvar = num_badvar,
  #                  #goodvar = newgoodvar,
  #                  badvar = badvar,
  #                  #goodmodelpart = newgoodmodelpart,
  #                  nextstep = 'no',
  #                  correlatedErrors = correlatedErrors)
  # }
  # 
  # # if(length(unique(unlist(newbadvar)))== 0){
  # #   if(length(setdiff(newgoodvar[[stepPrev$num_factor]],stepPrev$goodvar[[stepPrev$num_factor]])) == 0){
  # #     model <- paste0(crossloadmodel, collapse = '\n')
  # #     #add in provided list of correlated errors
  # #     if(!is.null(correlatedErrors)){
  # #       model <- paste0(model, '\n', correlatedErrors)
  # #     }
  # #   }
  # #   if(length(setdiff(newgoodvar[[stepPrev$num_factor]],stepPrev$goodvar[[stepPrev$num_factor]])) != 0 && 
  # #      length(unique(unlist(newbadvar)))!= 1){
  # #     model <- paste0(newgoodmodelpart, collapse = '\n')
  # #     #add in provided list of correlated errors
  # #     if(!is.null(correlatedErrors)){
  # #       model <- paste0(model, '\n', correlatedErrors)
  # #     }
  # #   }
  # #   fit <- miive(model, data, var.cov = T)
  # #   badvar <- getbadvar(fit, sigLevel)
  # #   num_badvar <- length(badvar)
  # #   
  # #   finalobj <- list(model = model,
  # #                    fit  = fit,
  # #                    num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
  # #                    num_badvar = num_badvar,
  # #                    #goodvar = newgoodvar,
  # #                    badvar = badvar,
  # #                    #goodmodelpart = newgoodmodelpart,
  # #                    nextstep = 'no',
  # #                    correlatedErrors = correlatedErrors)
  # # }
  # if(length(unique(unlist(newbadvar)))== 1){
  #   if(all(unlist(newgoodmodelpart) == unlist(stepPrev$goodmodelpart))){ #then we keep everything from the last step
  #     finalobj <- stepPrev
  #     finalobj$nextstep <- 'no' #need to stop the function from running
  #   }else{
  #     #we keep the single problematic variable on the latest created factor.
  #     newgoodmodelpart[[length(newgoodmodelpart)]] <- paste0(newgoodmodelpart[[length(newgoodmodelpart)]], '+', newbadvar)
  #     model <- paste0(newgoodmodelpart, collapse = '\n')
  #     #add in provided list of correlated errors
  #     if(!is.null(correlatedErrors)){
  #       model <- paste0(model, '\n', correlatedErrors)
  #     }
  #     fit <- miive(model, data, var.cov = T)
  #     badvar <- getbadvar(fit, sigLevel)
  #     num_badvar <- length(badvar)
  #     
  #     finalobj <- list(model = model,
  #                      fit  = fit,
  #                      num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
  #                      num_badvar = num_badvar,
  #                      #goodvar = newgoodvar,
  #                      badvar = badvar,
  #                      #goodmodelpart = newgoodmodelpart,
  #                      nextstep = 'no',
  #                      correlatedErrors = correlatedErrors)
  #   }
  #   
  # }
  # 
  # if(length(unique(unlist(newbadvar))) >1){
  #   
  #   stepPrev$badvar <- newbadvar
  #   stepPrev$goodmodelpart <- newgoodmodelpart
  #   
  #   scalingindicator <- select_scalingind_stepN(data, sigLevel, scalingCrit, stepPrev)
  #   order_scalingind <- which(newbadvar==scalingindicator)
  #   
  #   model <- paste(paste0(newgoodmodelpart, collapse = '\n'),
  #                  paste(paste0("f",stepPrev$num_factor+1), "=~",paste0(newbadvar[order_scalingind]), '+',
  #                        paste0(newbadvar[-order_scalingind], collapse = "+"), sep = ""),
  #                  sep = "\n")
  #   #add in provided list of correlated errors
  #   if(!is.null(correlatedErrors)){
  #     model <- paste0(model, '\n', correlatedErrors)
  #   }
  #   fit <- miive(model, data, var.cov = T)
  #   badvar <- getbadvar(fit, sigLevel)
  #   num_badvar <- length(badvar)
  #   
  #   num_factor <- stepPrev$num_factor
  #   #update goodvar and goodmodelpart
  #   newgoodvar[[num_factor+1]] <- setdiff(newbadvar, badvar)
  #   
  #   newgoodvar[[num_factor+1]] <- newgoodvar[[num_factor+1]][c(match(scalingindicator, newgoodvar[[num_factor+1]]),
  #                                                              setdiff(order(newgoodvar[[num_factor+1]]),match(scalingindicator, newgoodvar[[num_factor+1]])))]
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
  #                    nextstep = ifelse(length(badvar!=0), 'yes', 'no'),
  #                    correlatedErrors = correlatedErrors)
  #   
  #   
  #   
  # }
  # 
  
  
  
  return(finalobj)
}
