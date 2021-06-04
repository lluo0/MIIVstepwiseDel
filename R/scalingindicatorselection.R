select_scalingind <- function(data, threshold = .05,
                            priority = "order"){

  scalingindicator <- character()

  sargan_sig_badvar <- list()

  factor_loading_badvar <- list()

  ##order or R2
  R2_order<- colnames(r2_order(data))
  ##fit for each indicator as the scaling indicator
  model <- list()
  fit <- list()
  for(p in 1:dim(data)[2]){
    model[[p]] <- paste0('f1=~', paste0(colnames(data)[p]), '+', paste0(colnames(data)[-p], collapse = '+'))
    fit[[p]] <- miive(model[[p]], data, var.cov = T)
    names(model)[p] <- names(fit)[p] <- colnames(data)[p]
  }
  ##number of significant sargan for each variable as the scaling indicator
  for(p in 1:dim(data)[2]){
    sargan_sig_badvar[[p]] <- 0
    names(sargan_sig_badvar)[p] <- names(fit)[p]
    for(i in 1:length(fit[[p]]$eqn)){
      if(fit[[p]]$eqn[[i]]$sargan.p < threshold){
        sargan_sig_badvar[[p]] <- sargan_sig_badvar[[p]]+1
      }
    }
  }
  ##number of insignificant factor loading for each variable as the scaling indicator
  #estimatefittable <- list()
  for(p in 1:dim(data)[2]){
    factor_loading_badvar[[p]] <- 0
    #estimatefittable[[p]] <- vector()
    names(factor_loading_badvar)[p] <-  names(fit)[p]
    #names(estimatefittable)[p] <-
    #estimatefittable[[p]] <- estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",]
    factor_loading_badvar[[p]] <- length(which(estimatesTable(fit[[p]])[estimatesTable(fit[[p]])[,2] == "=~",7] > .05))

  }
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
    scalingindicator <- colnames(data)[1]
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
    scalingindicator <- sarganallmin[order(match(mixedallmin, R2_order))][1]
  }
  # else{
  #   stop('ERROR: please specify a valid order of criteria.')
  # }
  if(is.null(scalingindicator)){
    scalingindicator <- colnames(data)[1]
  }
  return(scalingindicator)
}
#
#
# ##for step 2
# select_scalingind_stepN <- function(data, threshold = .05,
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
#       if(fit[[p]]$eqn[[i]]$sargan.p < threshold){
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


##for step N
select_scalingind_stepN <- function(data, threshold = .05,
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
      if(fit[[p]]$eqn[[i]]$sargan.p < threshold){
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
    scalingindicator <- sarganallmin[order(match(mixedallmin, R2_order))][1]
  }
  # else{
  #   stop('ERROR: please specify a valid order of criteria.')
  # }
if(is.null(scalingindicator)){
  scalingindicator <- colnames(data[,badvar])[1]
}
  return(scalingindicator)
}



