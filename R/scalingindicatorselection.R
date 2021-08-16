#scaling indicator selection criterion

#1, order: uses the first appearing variable as the scaling indicator

#2, sargan: uses the variable with the least number of significant sargans as the scaling indicator.
#if multiple variables have the same least number of significant sargans, chooses the first appearing variable.

#3, R2: uses the variable with the highest R2 as the scaling indicator.
#if multiple variables have the same R2 value, chooses the first appearing variable.

#4, factor loading: uses the variable with the most number of significant factor loadings as the scaling indicator.
#if multiple variables have the same most number of significant factor loadings, chooses the first appearing variable.

#5, sargan_R2: uses the variable with the least number of significant sargans as the scaling indicator.
#if multiple variables have the same least number of significant sargans, chooses the one with higher R2.
#if still multiple options, chooses the first appearing variable.

#6, sargan_factorloading: uses the variable with the least number of significant sargans as the scaling indicator.
#if multiple variables have the same least number of significant sargans, chooses the one with more significant factor loadings.
#if still multiple options, chooses the first appearing variable.

#7, sargan_factorloading_R2: uses the variable with the least number of significant sargans as the scaling indicator.
#if multiple variables have the same least number of significant sargans, chooses the one with more significant factor loadings.
#if still multiple options, chooses the one with higher R2.

#8, factorloading_R2: uses the variable with the most number of significant factor loadings as the scaling indicator.
#if multiple variables have the same nist number of significant factor loadings, chooses the one with higher R2.
#if still multiple options, chooses the first appearing variable.

#9, factorloading_sargan: uses the variable with the most number of significant factor loadings as the scaling indicator.
#if multiple variables have the same nist number of significant factor loadings, chooses the one with less significant sargans.
#if still multiple options, chooses the first appearing variable.

#10, factorloading_sargan_R2: uses the variable with the most number of significant factor loadings as the scaling indicator.
#if multiple variables have the same nist number of significant factor loadings, chooses the one with less significant sargans.
#if still multiple options, chooses the one with higher R2.

#11, sargan+factorloading: uses the variable with the least sum of significant sargans and non-signficant factor loadings.
#if multiple variables have the same least sum of significnat sargans and non-significant factor loadings, chooses the first appearing variable.

#12, sargan+factorloading_R2: uses the variable with the least sum of significant sargans and non-signficant factor loadings.
#if multiple variables have the same least sum of significnat sargans and non-significant factor loadings, chooses the one with higher R2.


select_scalingind <- function(data, sigLevel = .05,
                              priority = "order"){

  scalingindicator <- character()

  num_sigsargan <- list()

  num_nonsigfactorloading <- list()

  #if scaling indicator selection is order, just use the first variable as the scaling indicator.
  if(priority == 'order'){
    scalingindicator <- colnames(data)[1]
  }
  #otherwise, need to run R2 values.
  if(priority != 'order'){
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
      num_nonsigfactorloading[[p]] <- length(which(estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",7] > .05))
    }

    ##calculate the least number of significant sargans, and least number of non-significant factor loadings.
    ##then see which variables meet the criteria.
    min_sarganbad <- min(unlist(num_sigsargan))
    min_factorbad <- min(unlist(num_nonsigfactorloading))
    sarganallmin <- colnames(t(which(num_sigsargan==min_sarganbad)))
    factorallmin <- colnames(t(which(num_nonsigfactorloading==min_factorbad)))
  }

  ##return the scaling indicator
  if(priority == 'sargan'){
    scalingindicator <- sarganallmin[1]
  }
  if(priority == 'R2'){
    scalingindicator <- R2_order[1]
  }
  if(priority == 'factorloading'){
    scalingindicator <- factorallmin[1]
  }
  if(priority == 'sargan_R2'){
    scalingindicator <- sarganallmin[order(match(sarganallmin, R2_order))][1]
  }
  if(priority == 'sargan_factorloading'){
    # mixedlist <- num_nonsigfactorloading[names(num_nonsigfactorloading)==sarganallmin]
    mixedlist <- num_nonsigfactorloading[sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
  }
  if(priority == 'sargan_factorloading_R2'){
    #mixedlist <- factor_loading_badvar[names(factor_loading_badvar)==sarganallmin]
    mixedlist <- num_nonsigfactorloading[sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[order(match(mixedallmin, R2_order))][1]
  }
  if(priority == 'factorloading_R2'){
    scalingindicator <- factorallmin[order(match(factorallmin, R2_order))][1]
  }
  if(priority == 'factorloading_sargan'){
    #mixedlist <- num_sigsargan[names(num_sigsargan)==factorallmin]
    mixedlist <- num_sigsargan[factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
  }
  if(priority == 'factorloading_sargan_R2'){
    #mixedlist <- num_sigsargan[names(num_sigsargan)==factorallmin]
    mixedlist <- num_sigsargan[factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[order(match(mixedallmin, R2_order))][1]
  }
  if(priority =='sargan+factorloading'){
    num_sum <- mapply("+", num_sigsargan, num_nonsigfactorloading, SIMPLIFY = FALSE)
    scalingindicator <- colnames(t(which(num_sum==min(unlist(num_sum)))))[1]
  }
  if(priority =='sargan+factorloading_R2'){
    num_sum <- mapply("+", num_sigsargan, num_nonsigfactorloading, SIMPLIFY = FALSE)
    min_sum <- min(unlist(num_sum))
    sumallmin <- colnames(t(which(num_sum==min_sum)))
    scalingindicator <- sumallmin[order(match(sumallmin, R2_order))][1]
  }
  # else{
  #   stop('ERROR: please specify a valid order of criteria.')
  # }
  # if(is.null(scalingindicator)){
  #   scalingindicator <- colnames(data)[1]
  # }
  return(scalingindicator)
}
#
#
# ##for step 2
# select_scalingind_stepN <- function(data, sigLevel = .05,
#                               priority = "order",prevStep){
#
#   scalingindicator <- character()
#
#   sargan_sig_badvar <- list()
#
#   factor_loading_badvar <- list()
#
#   ##order or R2
#   R2_order<- colnames(r2_order(data[,prevStep$badvar]))
#
#   ##new setup for newmodel
#   num_factor <- prevStep$num_factor+1
#   goodmodelpart <- prevStep$goodmodelpart
#   badvar <- prevStep$badvar
#
#   ##fit for each indicator as the scaling indicator
#   model <- list()
#   fit <- list()
#   for(p in 1:length(badvar)){
#     model[[p]] <- paste(paste0(goodmodelpart, collapse = '\n'),
#                         paste(paste0("f",num_factor), "=~",paste0(badvar[p]), '+',
#                               paste0(badvar[-p], collapse = "+"), sep = ""),
#                         sep = "\n")
#     fit[[p]] <- miive(model[[p]], data, var.cov = T)
#     names(model)[p] <- names(fit)[p] <- badvar[p]
#   }
#   ##number of significant sargan for each variable as the scaling indicator
#   for(p in 1:length(model)){
#     sargan_sig_badvar[[p]] <- 0
#     names(sargan_sig_badvar)[p] <- names(fit)[p]
#     for(i in 1:length(fit[[p]]$eqn)){
#       if(fit[[p]]$eqn[[i]]$sargan.p < sigLevel){
#         sargan_sig_badvar[[p]] <- sargan_sig_badvar[[p]]+1
#       }
#     }
#   }
#   #sargan_sig_badvar <- sargan_sig_badvar[badvar]
#   ##number of insignificant factor loading for each variable as the scaling indicator
#   #estimatefittable <- list()
#   for(p in 1:length(model)){
#     factor_loading_badvar[[p]] <- 0
#     #estimatefittable[[p]] <- vector()
#     names(factor_loading_badvar)[p] <-  names(fit)[p]
#     #names(estimatefittable)[p] <-
#     #estimatefittable[[p]] <- estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",]
#     factor_loading_badvar[[p]] <- length(which(estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",7] > .05))
#
#   }
#   #factor_loading_badvar <- factor_loading_badvar[badvar]
#   # return(sargan_sig_badvar)
#
#   ##then provide the scaling indicator based on  the criterion provided
#   ##some setup
#   min_sarganbad <- min(unlist(sargan_sig_badvar))
#   min_factorbad <- min(unlist(factor_loading_badvar))
#   sarganallmin <- colnames(t(which(sargan_sig_badvar==min_sarganbad)))
#   factorallmin <- colnames(t(which(factor_loading_badvar==min_factorbad)))
#   #candidate <- vector()
#   ##return the scaling indicator
#   if(priority == 'order'){
#     scalingindicator <- colnames(data)[1]
#   }
#   if(priority == 'sargan'){
#     scalingindicator <- sarganallmin[1]
#   }
#   if(priority == 'R2'){
#     scalingindicator <- R2_order[1]
#   }
#   if(priority == 'factorloading'){
#     scalingindicator <- factorallmin[1]
#   }
#   if(priority == 'sargan_R2'){
#     scalingindicator <- sarganallmin[order(match(sarganallmin, R2_order))][1]
#   }
#   if(priority == 'sargan_factorloading'){
#     mixedlist <- factor_loading_badvar[names(factor_loading_badvar)==sarganallmin]
#     mixedmin <- min(unlist(mixedlist))
#     mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
#     scalingindicator <- mixedallmin[1]
#   }
#   if(priority == 'sargan_factorloading_R2'){
#     mixedlist <- factor_loading_badvar[names(factor_loading_badvar)==sarganallmin]
#     mixedmin <- min(unlist(mixedlist))
#     mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
#     scalingindicator <- sarganallmin[order(match(mixedallmin, R2_order))][1]
#   }
#   if(priority == 'factorloading_R2'){
#     scalingindicator <- factorallmin[order(match(factorallmin, R2_order))][1]
#   }
#   if(priority == 'factorloading_sargan'){
#     mixedlist <- sargan_sig_badvar[names(sargan_sig_badvar)==factorallmin]
#     mixedmin <- min(unlist(mixedlist))
#     mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
#     scalingindicator <- mixedallmin[1]
#   }
#   if(priority == 'factorloading_sargan_R2'){
#     mixedlist <- sargan_sig_badvar[names(sargan_sig_badvar)==factorallmin]
#     mixedmin <- min(unlist(mixedlist))
#     mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
#     scalingindicator <- mixedallmin[1]
#     scalingindicator <- sarganallmin[order(match(mixedallmin, R2_order))][1]
#   }
#   # else{
#   #   stop('ERROR: please specify a valid order of criteria.')
#   # }
#
#   return(scalingindicator)
# }

select_scalingind_stepN <- function(data, sigLevel = .05,
                                    priority = "order", stepPrev){

  scalingindicator <- character()

  num_sigsargan <- list()

  num_nonsigfactorloading <- list()
  #extract info from  the previous step finalobj, aka stepPrev
  goodmodelpart <- stepPrev$goodmodelpart
  badvar <- stepPrev$badvar
  num_factor <- stepPrev$num_factor
  #if scaling indicator selection is order, just use the first variable as the scaling indicator.
  if(priority == 'order'){
    scalingindicator <- badvar[1]
  }
  #otherwise, need to run R2 values.
  if(priority != 'order'){
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

    ##calculate the least number of significant sargans, and least number of non-significant factor loadings.
    ##then see which variables meet the criteria.
    min_sarganbad <- min(unlist(num_sigsargan))
    min_factorbad <- min(unlist(num_nonsigfactorloading))
    sarganallmin <- colnames(t(which(num_sigsargan==min_sarganbad)))
    factorallmin <- colnames(t(which(num_nonsigfactorloading==min_factorbad)))
  }

  ##return the scaling indicator
  if(priority == 'sargan'){
    scalingindicator <- sarganallmin[1]
  }
  if(priority == 'R2'){
    scalingindicator <- R2_order[1]
  }
  if(priority == 'factorloading'){
    scalingindicator <- factorallmin[1]
  }
  if(priority == 'sargan_R2'){
    scalingindicator <- sarganallmin[order(match(sarganallmin, R2_order))][1]
  }
  if(priority == 'sargan_factorloading'){
    # mixedlist <- num_nonsigfactorloading[names(num_nonsigfactorloading)==sarganallmin]
    mixedlist <- num_nonsigfactorloading[sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
  }
  if(priority == 'sargan_factorloading_R2'){
    #mixedlist <- factor_loading_badvar[names(factor_loading_badvar)==sarganallmin]
    mixedlist <- num_nonsigfactorloading[sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[order(match(mixedallmin, R2_order))][1]
  }
  if(priority == 'factorloading_R2'){
    scalingindicator <- factorallmin[order(match(factorallmin, R2_order))][1]
  }
  if(priority == 'factorloading_sargan'){
    #mixedlist <- num_sigsargan[names(num_sigsargan)==factorallmin]
    mixedlist <- num_sigsargan[factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
  }
  if(priority == 'factorloading_sargan_R2'){
    #mixedlist <- num_sigsargan[names(num_sigsargan)==factorallmin]
    mixedlist <- num_sigsargan[factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[order(match(mixedallmin, R2_order))][1]
  }
  if(priority =='sargan+factorloading'){
    num_sum <- mapply("+", num_sigsargan, num_nonsigfactorloading, SIMPLIFY = FALSE)
    scalingindicator <- colnames(t(which(num_sum==min(unlist(num_sum)))))[1]
  }
  if(priority =='sargan+factorloading_R2'){
    num_sum <- mapply("+", num_sigsargan, num_nonsigfactorloading, SIMPLIFY = FALSE)
    min_sum <- min(unlist(num_sum))
    sumallmin <- colnames(t(which(num_sum==min_sum)))
    scalingindicator <- sumallmin[order(match(sumallmin, R2_order))][1]
  }
  # else{
  #   stop('ERROR: please specify a valid order of criteria.')
  # }
  # if(is.null(scalingindicator)){
  #   scalingindicator <- colnames(data)[1]
  # }
  return(scalingindicator)
}

##for step N
select_scalingind_stepN <- function(data, sigLevel = .05,
                                    priority = "order", goodmodelpart, badvar, num_factor){

  #reorder badvar based on their order in the original data
  badvar <- badvar[order(match(badvar, colnames(data)))]
  #setup
  scalingindicator <- character()

  sargan_sig_badvar <- list()

  factor_loading_badvar <- list()

  ##order or R2
  R2_order<- colnames(r2_order(data[badvar]))

  ##new setup for newmodel


  ##fit for each indicator as the scaling indicator
  model <- list()
  fit <- list()
  for(p in 1:length(badvar)){
    model[[p]] <- paste(paste0(goodmodelpart, collapse = '\n'),
                        paste(paste0("f",num_factor), "=~",paste0(badvar[p]), '+',
                              paste0(badvar[-p], collapse = "+"), sep = ""),
                        sep = "\n")
    fit[[p]] <- miive(model[[p]], data, var.cov = T)
    names(model)[p] <- names(fit)[p] <- badvar[p]
  }
  ##number of significant sargan for each variable as the scaling indicator
  for(p in 1:length(model)){
    sargan_sig_badvar[[p]] <- 0
    names(sargan_sig_badvar)[p] <- names(fit)[p]
    for(i in 1:length(fit[[p]]$eqn)){
      if(fit[[p]]$eqn[[i]]$sargan.p < sigLevel){
        sargan_sig_badvar[[p]] <- sargan_sig_badvar[[p]]+1
      }
    }
  }
  #sargan_sig_badvar <- sargan_sig_badvar[badvar]
  ##number of insignificant factor loading for each variable as the scaling indicator
  #estimatefittable <- list()
  for(p in 1:length(model)){
    factor_loading_badvar[[p]] <- 0
    #estimatefittable[[p]] <- vector()
    names(factor_loading_badvar)[p] <-  names(fit)[p]
    #names(estimatefittable)[p] <-
    #estimatefittable[[p]] <- estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",]
    factor_loading_badvar[[p]] <- length(which(estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",7] > .05))

  }
  #factor_loading_badvar <- factor_loading_badvar[badvar]
  # return(sargan_sig_badvar)

  ##then provide the scaling indicator based on  the criterion provided
  ##some setup
  min_sarganbad <- min(unlist(sargan_sig_badvar))
  min_factorbad <- min(unlist(factor_loading_badvar))
  sarganallmin <- colnames(t(which(sargan_sig_badvar==min_sarganbad)))
  factorallmin <- colnames(t(which(factor_loading_badvar==min_factorbad)))
  #candidate <- vector()
  ##return the scaling indicator
  if(priority == 'order'){
    scalingindicator <- colnames(data[,badvar])[1]
  }
  if(priority == 'sargan'){
    scalingindicator <- sarganallmin[1]
  }
  if(priority == 'R2'){
    scalingindicator <- R2_order[1]
  }
  if(priority == 'factorloading'){
    scalingindicator <- factorallmin[1]
  }
  if(priority == 'sargan_R2'){
    scalingindicator <- sarganallmin[order(match(sarganallmin, R2_order))][1]
  }
  if(priority == 'sargan_factorloading'){
    mixedlist <- factor_loading_badvar[names(factor_loading_badvar)==sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
  }
  if(priority == 'sargan_factorloading_R2'){
    mixedlist <- factor_loading_badvar[names(factor_loading_badvar)==sarganallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- sarganallmin[order(match(mixedallmin, R2_order))][1]
  }
  if(priority == 'factorloading_R2'){
    scalingindicator <- factorallmin[order(match(factorallmin, R2_order))][1]
  }
  if(priority == 'factorloading_sargan'){
    mixedlist <- sargan_sig_badvar[names(sargan_sig_badvar)==factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
  }
  if(priority == 'factorloading_sargan_R2'){
    mixedlist <- sargan_sig_badvar[names(sargan_sig_badvar)==factorallmin]
    mixedmin <- min(unlist(mixedlist))
    mixedallmin <- colnames(t(which(mixedlist==mixedmin)))
    scalingindicator <- mixedallmin[1]
    scalingindicator <- mixedallmin[order(match(mixedallmin, R2_order))][1]
  }
  # else{
  #   stop('ERROR: please specify a valid order of criteria.')
  # }
  if(is.null(scalingindicator)){
    scalingindicator <- colnames(data[,badvar])[1]
  }
  return(scalingindicator)
}



