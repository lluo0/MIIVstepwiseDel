stepN_E5 <- function(stepPrev, data, threshold, criterion){
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
    ##then same as before
    r2 <- order_r2(subset(data, select = newbadvar))
    models <- list()
    fits <- list()
    num_factor <- stepPrev$num_factor+1
    for (i in 1:length(r2)){
      models[[i]] <- paste(paste0(newgoodmodelpart, collapse = '\n'),
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
                       num_badvar = num_badvar,
                       nextstep = 'no')
    }
    if(num_badvar>0){
      #save the problematic variables
      badvar <- badvar_list[[best_num]]
      #save the good variables - update the ones from the previous output
      goodvar <- stepPrev$goodvar
      #goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)
      goodvar[[num_factor]] <- setdiff(r2[[best_num]], badvar)

      goodmodelpart <- c(newgoodmodelpart,
                         paste(paste0("f", num_factor), "=~",
                               paste(goodvar[[num_factor]], collapse = "+"), sep = ""))
      finalobj <- list(model = models[[best_num]],
                       fit  = fits[[best_num]],
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



