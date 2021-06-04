getbadvar <- function(fit, threshold=.05){
  v_list <- vector()
  for (p in 1:length(fit$eqn))
    if (fit$eqn[[p]]$sargan.p < threshold){
      v_list <- append(v_list,fit$eqn[[p]]$DVobs)
    }
  if(length(v_list)==0){
    v_list_final <- NULL
  }else{
    v_list2 <- vector()
    table <- na.omit(estimatesTable(fit))
    for (p in 1:length(fit$eqn))
      if (table[p,7] > threshold){
        v_list2 <- append(v_list2, table[p,3])
      }
    v_list_final <- unique(c(v_list, v_list2))
  }

  return(v_list_final)
}

getbadvar_crossload <- function(fit, threshold=.05, num_fac, badvar){
  newbadvar_coef <- list()
  coeftable <- estimatesTable(fit)[estimatesTable(fit)[,2] == '=~',]
  #the bad var with non-siginificant coefficients for each factor
  for(p in 1:num_fac){
    newbadvar_coef[[p]] <- vector()
    for(i in 1:nrow(coeftable)){
      if(paste0('f',p)%in%coeftable[i,1]
         & coeftable[i,3]%in% badvar
         & coeftable[i,7] > threshold){
        newbadvar_coef[[p]] <- append(newbadvar_coef[[p]],coeftable[i,3])
      }
    }
  }
  #the bad var with significant sargans.
  newbadvar_sargan <- vector()
  for (p in 1:length(fit$eqn))
    if (fit$eqn[[p]]$sargan.p < threshold & fit$eqn[[p]]$DVobs %in% badvar){
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
