
library(lavaan)
library(MIIVsem)
# getbadvar <- function(fit, sigLevel=.05){
#   v_list <- vector()
#   for (p in 1:length(fit$eqn))
#     if (fit$eqn[[p]]$sargan.p < sigLevel){
#       v_list <- append(v_list,fit$eqn[[p]]$DVobs)
#     }
#   if(length(v_list)==0){
#     v_list_final <- NULL
#   }else{
#     v_list2 <- vector()
#     table <- na.omit(estimatesTable(fit))
#     for (p in 1:length(fit$eqn))
#       if (table[p,7] > sigLevel){
#         v_list2 <- append(v_list2, table[p,3])
#       }
#     v_list_final <- unique(c(v_list, v_list2))
#   }
#
#   return(v_list_final)
# }

# all.identical <- function(list) {
#   all(mapply(identical, head(list, 1), tail(list, -1)))
# }
#
# all.identical <- function(list) {
#   mapply(identical, head(list, 1), tail(list, -1))
# }
#

getbadvar <- function(fit, sigLevel=.05){
  v_list <- vector()
  for (p in 1:length(fit$eqn))
    if (fit$eqn[[p]]$sargan.p < sigLevel){
      v_list <- append(v_list,fit$eqn[[p]]$DVobs)
    }

  v_list2 <- vector()
  table <- na.omit(estimatesTable(fit))
  for (p in 1:length(fit$eqn))
    if (table[p,7] > sigLevel){
      v_list2 <- append(v_list2, table[p,3])
    }
  v_list_final <- unique(c(v_list, v_list2))


  return(v_list_final)
}

getbadvar_crossload <- function(fit, sigLevel=.05, num_fac, badvar){
  newbadvar_coef <- list()
  coeftable <- estimatesTable(fit)[estimatesTable(fit)[,2] == '=~',]
  #the bad var with non-siginificant coefficients for each factor
  for(p in 1:num_fac){
    newbadvar_coef[[p]] <- vector()
    for(i in 1:nrow(coeftable)){
      if(paste0('f',p)%in%coeftable[i,1]
         & coeftable[i,3]%in% badvar
         & coeftable[i,7] > sigLevel){
        newbadvar_coef[[p]] <- append(newbadvar_coef[[p]],coeftable[i,3])
      }
    }
  }
  #the bad var with significant sargans.
  newbadvar_sargan <- vector()
  for (p in 1:length(fit$eqn))
    if (fit$eqn[[p]]$sargan.p < sigLevel & fit$eqn[[p]]$DVobs %in% badvar){
      newbadvar_sargan <- append(newbadvar_sargan,fit$eqn[[p]]$DVobs)
    }
  newbadvar <- lapply(newbadvar_coef, function(x)
    unique(c(x, newbadvar_sargan)))
  return(newbadvar)
}

order_r2 <- function(object){
  r2 <- matrix(NA, nrow = 1, ncol = dim(object)[2])
  colnames(r2) <- colnames(object)
  for (i in 1:dim(object)[2]){
    r2[,i] <- summary(lm(paste(colnames(object)[i], paste(colnames(object)[-i], collapse = "+"), sep = "~"), data = object))$r.squared
  }
  r2 <- as.matrix(t(r2[,order(r2[nrow(r2),],decreasing=TRUE)]))
  #return(r2)
  r2_list <- list()
  for (i in 1:length(colnames(r2))-1){
    r2_list[[1]] <- colnames(r2)
    r2_list[[i+1]] <- c(colnames(r2)[-c(1:i)], colnames(r2)[c(1:i)])
  }
  return(r2_list)
}

r2_order <- function(object){
  ##check the r2 for each variable and use the highest r2 as the initial scaling indicator
  r2 <- matrix(NA, nrow = 1, ncol = dim(object)[2])
  colnames(r2) <- colnames(object)
  for (i in 1:dim(object)[2]){
    r2[,i] <- summary(lm(paste(colnames(object)[i], paste(colnames(object)[-i], collapse = "+"), sep = "~"), data = object))$r.squared
  }
  r2 <- as.matrix(t(r2[,order(r2[nrow(r2),],decreasing=TRUE)]))
  return(r2)}

select_scalingind <- function(data, sigLevel = .05,
                              scalingCrit = "order"){

  scalingindicator <- character()

  num_sigsargan <- list()

  num_nonsigfactorloading <- list()

  #if scaling indicator selection is order, just use the first variable as the scaling indicator.
  if(scalingCrit == 'order'){
    scalingindicator <- colnames(data)[1]
  }
  #otherwise, need to run R2 values.
  if(scalingCrit != 'order'){
    ##order or R2
    R2_order <- colnames(r2_order(data))

    ##fit for each indicator as the scaling indicator
    model <- list()
    fit <- list()
    for(p in 1:ncol(data)){
      model[[p]] <- paste0('f1=~', paste0(colnames(data)[p]), '+', paste0(colnames(data)[-p], collapse = '+'))
      fit[[p]] <- miive(model[[p]], data, var.cov = T)
      names(model)[p] <- names(fit)[p] <- colnames(data)[p] #name the lists using the variable used as the scaling indicator
    }

    ##calculate the number of significant sargan for each variable as the scaling indicator
    for(p in 1:ncol(data)){
      num_sigsargan[[p]] <- 0
      names(num_sigsargan)[p] <- names(fit)[p]
      for(i in 1:length(fit[[p]]$eqn)){
        if(fit[[p]]$eqn[[i]]$sargan.p < sigLevel){
          num_sigsargan[[p]] <- num_sigsargan[[p]]+1
        }
      }
    }

    ##calculate the number of insignificant factor loading for each variable as the scaling indicator
    for(p in 1:dim(data)[2]){
      num_nonsigfactorloading[[p]] <- 0
      #estimatefittable[[p]] <- vector()
      names(num_nonsigfactorloading)[p] <-  names(fit)[p]
      #names(estimatefittable)[p] <-
      #estimatefittable[[p]] <- estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",]
      num_nonsigfactorloading[[p]] <- length(which(estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",7] > sigLevel))
    }

    var_sigsargan <- list()
    var_nonsigfactorloading <- list()
    ##sargan+scaling indicator
    ## see negfisher[[12]] for example - a variable with both significant sargan and non-significant factor loading will only be counted once
    for(p in 1:dim(data)[2]){
      var_sigsargan[[p]] <- var_nonsigfactorloading[[p]] <- vector()
      names(var_sigsargan)[p] <- names(var_nonsigfactorloading)[p] <-  names(fit)[p]
      for(i in 1:length(fit[[p]]$eqn)){
        if(fit[[p]]$eqn[[i]]$sargan.p < sigLevel){
          var_sigsargan[[p]] <- c(var_sigsargan[[p]], fit[[p]]$eqn[[i]]$DVobs)
        }
      }

      loadingtable <- estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",]
      for(i in 1:nrow(loadingtable)){
        if(loadingtable[i, 7] > sigLevel && !is.na(loadingtable[i,7])){
          var_nonsigfactorloading[[p]] <- c(var_nonsigfactorloading[[p]], loadingtable[i,3])
        }
      }
    }

    var_sarganANDfactorloading <- list()
    for(p in 1:length(var_sigsargan)){
      var_sarganANDfactorloading[[p]] <- vector()
      names(var_sarganANDfactorloading)[p] <- names(fit)[p]
      diffvar <- setdiff(var_sigsargan[[p]], var_nonsigfactorloading[[p]]) #this gets variables only with sig sargans and not with nonsig loadings
      var_sarganANDfactorloading[[p]] <- c(var_nonsigfactorloading[[p]], diffvar)
    }
    ##for sargan+factorloading_R2
    min_sarganANDloading <- min(sapply(var_sarganANDfactorloading, length))
    sarganANDloadingmin <- which(sapply(var_sarganANDfactorloading, length) == min_sarganANDloading)

    ##calculate the least number of significant sargans, and least number of non-significant factor loadings.
    ##then see which variables meet the criteria.
    min_sarganbad <- min(unlist(num_sigsargan))
    min_factorbad <- min(unlist(num_nonsigfactorloading))
    sarganallmin <- colnames(t(which(num_sigsargan==min_sarganbad)))
    factorallmin <- colnames(t(which(num_nonsigfactorloading==min_factorbad)))


  }

  ##return the scaling indicator
  if(scalingCrit == 'sargan'){
    scalingindicator <- sarganallmin[1]
  }
  if(scalingCrit == 'R2'){
    scalingindicator <- R2_order[1]
  }
  if(scalingCrit == 'factorloading'){
    scalingindicator <- factorallmin[1]
  }
  if(scalingCrit == 'sargan_R2'){
    scalingindicator <- sarganallmin[order(match(sarganallmin, R2_order))][1]
  }
  if(scalingCrit == 'sargan_factorloading'){
    # mixedlist <- num_nonsigfactorloading[names(num_nonsigfactorloading)==sarganallmin]
    mixedlist <- num_nonsigfactorloading[sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]

  }
  if(scalingCrit == 'sargan_factorloading_R2'){
    #mixedlist <- factor_loading_badvar[names(factor_loading_badvar)==sarganallmin]
    mixedlist <- num_nonsigfactorloading[sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))

  }
  if(scalingCrit == 'factorloading_R2'){
    scalingindicator <- factorallmin[order(match(factorallmin, R2_order))][1]
  }
  if(scalingCrit == 'factorloading_sargan'){
    #mixedlist <- num_sigsargan[names(num_sigsargan)==factorallmin]
    mixedlist <- num_sigsargan[factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
  }
  if(scalingCrit == 'factorloading_sargan_R2'){
    #mixedlist <- num_sigsargan[names(num_sigsargan)==factorallmin]
    mixedlist <- num_sigsargan[factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[order(match(mixedallmin, R2_order))][1]
  }
  if(scalingCrit =='sargan+factorloading'){
    # num_sum <- mapply("+", num_sigsargan, num_nonsigfactorloading, SIMPLIFY = FALSE)
    # scalingindicator <- colnames(t(which(num_sum==min(unlist(num_sum)))))[1]
    scalingindicator <- colnames(data)[sarganANDloadingmin][1]
  }
  if(scalingCrit =='sargan+factorloading_R2'){
    # num_sum <- mapply("+", num_sigsargan, num_nonsigfactorloading, SIMPLIFY = FALSE)
    # min_sum <- min(unlist(num_sum))
    # sumallmin <- colnames(t(which(num_sum==min_sum)))
    # scalingindicator <- sumallmin[order(match(sumallmin, R2_order))][1]
    scalingindicator <- colnames(data)[sarganANDloadingmin][order(match(colnames(data)[sarganANDloadingmin], R2_order))][1]
  }
  # else{
  #   stop('ERROR: please specify a valid order of criteria.')
  # }
  # if(is.null(scalingindicator)){
  #   scalingindicator <- colnames(data)[1]
  # }
  return(scalingindicator)
}


select_scalingind_stepN <- function(data, sigLevel = .05,
                                    scalingCrit = "order", stepPrev){

  scalingindicator <- character()

  num_sigsargan <- list()

  num_nonsigfactorloading <- list()
  #extract info from  the previous step finalobj, aka stepPrev
  goodmodelpart <- stepPrev$goodmodelpart
  badvar <- stepPrev$badvar
  num_factor <- stepPrev$num_factor
  #if scaling indicator selection is order, just use the first variable as the scaling indicator.
  if(scalingCrit == 'order'){
    scalingindicator <- badvar[1]
  }
  #otherwise, need to run R2 values.
  if(scalingCrit != 'order'){
    ##order or R2
    R2_order <- colnames(r2_order(data[,badvar]))

    ##fit for each indicator as the scaling indicator
    model <- list()
    fit <- list()
    # for(p in 1:ncol(data)){
    #   model[[p]] <- paste0('f1=~', paste0(colnames(data)[p]), '+', paste0(colnames(data)[-p], collapse = '+'))
    #   fit[[p]] <- miive(model[[p]], data, var.cov = T)
    #   names(model)[p] <- names(fit)[p] <- colnames(data)[p] #name the lists using the variable used as the scaling indicator
    # }
    for(p in 1:length(badvar)){
      model[[p]] <- paste(paste0(goodmodelpart, collapse = '\n'),
                          paste(paste0("f",num_factor+1), "=~",paste0(badvar[p]), '+',
                                paste0(badvar[-p], collapse = "+"), sep = ""),
                          sep = "\n")
      fit[[p]] <- miive(model[[p]], data, var.cov = T)
      names(model)[p] <- names(fit)[p] <- badvar[p]
    }

    ##calculate the number of significant sargan for each variable as the scaling indicator
    for(p in 1:length(badvar)){
      num_sigsargan[[p]] <- 0
      names(num_sigsargan)[p] <- names(fit)[p]
      for(i in 1:length(fit[[p]]$eqn)){
        if(fit[[p]]$eqn[[i]]$sargan.p < sigLevel){
          num_sigsargan[[p]] <- num_sigsargan[[p]]+1
        }
      }
    }

    ##calculate the number of insignificant factor loading for each variable as the scaling indicator
    for(p in 1:length(badvar)){
      num_nonsigfactorloading[[p]] <- 0
      #estimatefittable[[p]] <- vector()
      names(num_nonsigfactorloading)[p] <-  names(fit)[p]
      #names(estimatefittable)[p] <-
      #estimatefittable[[p]] <- estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",]
      num_nonsigfactorloading[[p]] <- length(which(estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",7] > .05))
    }

    var_sigsargan <- list()
    var_nonsigfactorloading <- list()
    ##sargan+scaling indicator
    ## see negfisher[[12]] for example - a variable with both significant sargan and non-significant factor loading will only be counted once
    for(p in 1:length(badvar)){
      var_sigsargan[[p]] <- var_nonsigfactorloading[[p]] <- vector()
      names(var_sigsargan)[p] <- names(var_nonsigfactorloading)[p] <-  names(fit)[p]
      for(i in 1:length(fit[[p]]$eqn)){
        if(fit[[p]]$eqn[[i]]$sargan.p < sigLevel){
          var_sigsargan[[p]] <- c(var_sigsargan[[p]], fit[[p]]$eqn[[i]]$DVobs)
        }
      }

      loadingtable <- estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",]
      for(i in 1:nrow(loadingtable)){
        if(loadingtable[i, 7] > sigLevel && !is.na(loadingtable[i,7])){
          var_nonsigfactorloading[[p]] <- c(var_nonsigfactorloading[[p]], loadingtable[i,3])
        }
      }
    }

    var_sarganANDfactorloading <- list()
    for(p in 1:length(var_sigsargan)){
      var_sarganANDfactorloading[[p]] <- vector()
      names(var_sarganANDfactorloading)[p] <- names(fit)[p]
      diffvar <- setdiff(var_sigsargan[[p]], var_nonsigfactorloading[[p]]) #this gets variables only with sig sargans and not with nonsig loadings
      var_sarganANDfactorloading[[p]] <- c(var_nonsigfactorloading[[p]], diffvar)
    }
    ##for sargan+factorloading_R2
    min_sarganANDloading <- min(sapply(var_sarganANDfactorloading, length))
    sarganANDloadingmin <- which(sapply(var_sarganANDfactorloading, length) == min_sarganANDloading)

    ##calculate the least number of significant sargans, and least number of non-significant factor loadings.
    ##then see which variables meet the criteria.
    min_sarganbad <- min(unlist(num_sigsargan))
    min_factorbad <- min(unlist(num_nonsigfactorloading))
    sarganallmin <- colnames(t(which(num_sigsargan==min_sarganbad)))
    factorallmin <- colnames(t(which(num_nonsigfactorloading==min_factorbad)))
  }

  ##return the scaling indicator
  if(scalingCrit == 'sargan'){
    scalingindicator <- sarganallmin[1]
  }
  if(scalingCrit == 'R2'){
    scalingindicator <- R2_order[1]
  }
  if(scalingCrit == 'factorloading'){
    scalingindicator <- factorallmin[1]
  }
  if(scalingCrit == 'sargan_R2'){
    scalingindicator <- sarganallmin[order(match(sarganallmin, R2_order))][1]
  }
  if(scalingCrit == 'sargan_factorloading'){
    # mixedlist <- num_nonsigfactorloading[names(num_nonsigfactorloading)==sarganallmin]
    # mixedlist <- num_nonsigfactorloading[sarganallmin]
    # mixedmin <- min(unlist(mixedlist))
    # mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    # scalingindicator <- mixedallmin[1]
    scalingindicator <- colnames(data)[]
  }
  if(scalingCrit == 'sargan_factorloading_R2'){
    #mixedlist <- factor_loading_badvar[names(factor_loading_badvar)==sarganallmin]
    mixedlist <- num_nonsigfactorloading[sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[order(match(mixedallmin, R2_order))][1]
  }
  if(scalingCrit == 'factorloading_R2'){
    scalingindicator <- factorallmin[order(match(factorallmin, R2_order))][1]
  }
  if(scalingCrit == 'factorloading_sargan'){
    #mixedlist <- num_sigsargan[names(num_sigsargan)==factorallmin]
    mixedlist <- num_sigsargan[factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
  }
  if(scalingCrit == 'factorloading_sargan_R2'){
    #mixedlist <- num_sigsargan[names(num_sigsargan)==factorallmin]
    mixedlist <- num_sigsargan[factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[order(match(mixedallmin, R2_order))][1]
  }
  if(scalingCrit =='sargan+factorloading'){
    # num_sum <- mapply("+", num_sigsargan, num_nonsigfactorloading, SIMPLIFY = FALSE)
    # scalingindicator <- colnames(t(which(num_sum==min(unlist(num_sum)))))[1]
    # scalingindicator <- colnames(data)[sarganANDloadingmin][1]

    scalingindicator <- rownames(as.data.frame(sarganANDloadingmin))[1] #can't just do colnames because now we are only considering from badvars
  }
  if(scalingCrit =='sargan+factorloading_R2'){
    # num_sum <- mapply("+", num_sigsargan, num_nonsigfactorloading, SIMPLIFY = FALSE)
    # min_sum <- min(unlist(num_sum))
    # sumallmin <- colnames(t(which(num_sum==min_sum)))
    # scalingindicator <- sumallmin[order(match(sumallmin, R2_order))][1]
    scalingindicator <- rownames(as.data.frame(sarganANDloadingmin))[order(match(colnames(data)[sarganANDloadingmin], R2_order))][1]
  }
  # else{
  #   stop('ERROR: please specify a valid order of criteria.')
  # }
  # if(is.null(scalingindicator)){
  #   scalingindicator <- colnames(data)[1]
  # }
  return(scalingindicator)
}

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

  if(num_badvar <=1){
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = 1,
                     num_badvar = num_badvar,
                     nextstep = 'no')
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


# stepN_E5 <- function(stepPrev, data, sigLevel, scalingCrit){
#
#   correlatedErrors <- stepPrev$correlatedErrors
#
#   ##first crossload the bad variables
#   crossloadmodel <- lapply(stepPrev$goodmodelpart, function(x)
#     paste0(x, '+',paste0(stepPrev$badvar, collapse = '+')))
#
#   #then add in list of correlated errors
#   crossloadmodel <- lapply(crossloadmodel, function(x)
#     paste0(x, '\n', correlatedErrors))
#
#   crossloadfit <- miive(paste0(crossloadmodel, collapse = '\n'), data, var.cov = T)
#   ##then see if any variables actually crossload
#   ##update the bad variable list
#   newbadvar <- getbadvar_crossload(crossloadfit, sigLevel, stepPrev$num_factor, stepPrev$badvar)
#
#
#   # #check if crossloads doesn't work at all
#   # crossloadcheck <- mapply(function(x) all(x == stepPrev$badvar), newbadvar, SIMPLIFY = T)
#   # crossloadcheck_TF <- all(crossloadcheck == T)
#   # #if crossloadcheck_TF == T, means
#
#
#   ##add the crossloaded variables to the model
#   # newgoodvar_addon <- lapply(newbadvar, function(x) setdiff(x, stepPrev$badvar))
#   newgoodvar_addon <- list()
#   for(p in 1:length(newbadvar)){
#     if(length(newbadvar[[p]])==0){
#       newgoodvar_addon[[p]] <- stepPrev$badvar
#     }
#
#     if(length(newbadvar[[p]])!=0){
#       newgoodvar_addon[[p]] <- setdiff(stepPrev$badvar,newbadvar[[p]])
#     }
#   }
#
#
#
#
#   # #create new badvar
#   # newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))
#
#   newgoodmodelpart <- stepPrev$goodmodelpart
#   for(p in 1:length(newgoodvar_addon)){
#     if(length(newgoodvar_addon[[p]]!=0)){
#       newgoodmodelpart[[p]] <- paste0(newgoodmodelpart[[p]], '+',
#                                       paste0(newgoodvar_addon[[p]], collapse = '+'))
#     }
#   }
#   newgoodvar <- mapply(c, stepPrev$goodvar, newgoodvar_addon, SIMPLIFY  = F)
#
#   #check if need a new factor
#   newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))
#
#
#   if(length(unique(unlist(newbadvar)))== 0){
#     if(length(setdiff(newgoodvar[[stepPrev$num_factor]],stepPrev$goodvar[[stepPrev$num_factor]])) == 0){
#       model <- paste0(crossloadmodel, collapse = '\n')
#       #add in provided list of correlated errors
#       if(!is.null(correlatedErrors)){
#         model <- paste0(model, '\n', correlatedErrors)
#       }
#     }
#     if(length(setdiff(newgoodvar[[stepPrev$num_factor]],stepPrev$goodvar[[stepPrev$num_factor]])) != 0){
#       model <- paste0(newgoodmodelpart, collapse = '\n')
#       #add in provided list of correlated errors
#       if(!is.null(correlatedErrors)){
#         model <- paste0(model, '\n', correlatedErrors)
#       }
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
#   if(length(unique(unlist(newbadvar)))== 1){
#     if(all(unlist(newgoodmodelpart) == unlist(stepPrev$goodmodelpart))){ #then we keep everything from the last step
#       finalobj <- stepPrev
#       finalobj$nextstep <- 'no' #need to stop the function from running
#     }else{
#       #we keep the single problematic variable on the latest created factor.
#       newgoodmodelpart[[length(newgoodmodelpart)]] <- paste0(newgoodmodelpart[[length(newgoodmodelpart)]], '+', newbadvar)
#       model <- paste0(newgoodmodelpart, collapse = '\n')
#       #add in provided list of correlated errors
#       if(!is.null(correlatedErrors)){
#         model <- paste0(model, '\n', correlatedErrors)
#       }
#       fit <- miive(model, data, var.cov = T)
#       badvar <- getbadvar(fit, sigLevel)
#       num_badvar <- length(badvar)
#
#       finalobj <- list(model = model,
#                        fit  = fit,
#                        num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
#                        num_badvar = num_badvar,
#                        #goodvar = newgoodvar,
#                        badvar = badvar,
#                        #goodmodelpart = newgoodmodelpart,
#                        nextstep = 'no',
#                        correlatedErrors = correlatedErrors)
#     }
#
#   }
#
#   if(length(unique(unlist(newbadvar))) >1){
#
#     stepPrev$badvar <- newbadvar
#     stepPrev$goodmodelpart <- newgoodmodelpart
#
#     scalingindicator <- select_scalingind_stepN(data, sigLevel, scalingCrit, stepPrev)
#     order_scalingind <- which(newbadvar==scalingindicator)
#
#     model <- paste(paste0(newgoodmodelpart, collapse = '\n'),
#                    paste(paste0("f",stepPrev$num_factor+1), "=~",paste0(newbadvar[order_scalingind]), '+',
#                          paste0(newbadvar[-order_scalingind], collapse = "+"), sep = ""),
#                    sep = "\n")
#     #add in provided list of correlated errors
#     if(!is.null(correlatedErrors)){
#       model <- paste0(model, '\n', correlatedErrors)
#     }
#     fit <- miive(model, data, var.cov = T)
#     badvar <- getbadvar(fit, sigLevel)
#     num_badvar <- length(badvar)
#
#     num_factor <- stepPrev$num_factor
#     #update goodvar and goodmodelpart
#     newgoodvar[[num_factor+1]] <- setdiff(newbadvar, badvar)
#
#     newgoodvar[[num_factor+1]] <- newgoodvar[[num_factor+1]][c(match(scalingindicator, newgoodvar[[num_factor+1]]),
#                                                                setdiff(order(newgoodvar[[num_factor+1]]),match(scalingindicator, newgoodvar[[num_factor+1]])))]
#
#
#     newgoodmodelpart[[num_factor+1]] <- paste(paste0("f", num_factor+1), "=~",
#                                               paste(newgoodvar[[num_factor+1]], collapse = "+"), sep = "")
#
#     finalobj <- list(model = model,
#                      fit  = fit,
#                      num_factor = stepPrev$num_factor+1,
#                      num_badvar = num_badvar,
#                      goodvar = newgoodvar,
#                      badvar = badvar,
#                      goodmodelpart = newgoodmodelpart,
#                      nextstep = ifelse(length(badvar!=0), 'yes', 'no'),
#                      correlatedErrors = correlatedErrors)
#
#
#
#   }
#
#
#
#
#   return(finalobj)
# }


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
      ##!! even for this one variable, its significance could still change and becomes problematic variable after removing it from other factors.
      ##an example would be negfisher[[13]]
      NEWbadvar <- getbadvar(newfit, sigLevel)
      if(NEWbadvar == unique(unlist(newbadvar)) && length(NEWbadvar)!=0){
        newmodel <- stepPrev$model
        newfit <- stepPrev$fit
      }
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
                     num_badvar = length(badvar),
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
    if(!is.null(correlatedErrors)){
      model <- c(model, correlatedErrors)
    }
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


EFAmiive5 <- function(data, sigLevel = .05, scalingCrit = 'order', correlatedErrors = NULL){
  step1 <- step1_E5(data, sigLevel, scalingCrit, correlatedErrors)

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

