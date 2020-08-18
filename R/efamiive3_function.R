###########function efamiive3##########
##function 1
##check the r2 for each variable and use the highest r2 as the initial scaling indicator
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
badVs <- function(fit){
  v_list <- vector()
  for (p in 1:length(fit$eqn))
    if (fit$eqn[[p]]$sargan.p < .05){
      v_list <- append(v_list,fit$eqn[[p]]$DVobs)
    }
  v_list2 <- vector()
  table <- na.omit(estimatesTable(fit))
  for (p in 1:length(fit$eqn))
    if (table[p,7] > .05){
      v_list2 <- append(v_list2, table[p,3])
    }
  v_list_final <- unique(c(v_list, v_list2))
  return(v_list_final)
}
EFAmiive3 <- function(object){
  #step 1

  #calculate the r2 for each of the variables and order them from high to low
  #this in ietegrated in the function "order_r2"
  r2 <- order_r2(object)
  #then run one factor model on the r2 list - aka each variable used once as the scaling indicator
  model_1 <- list()
  fit_1 <- list()
  for (i in 1:length(r2)){
    model_1[[i]] <- paste("f1=~", paste(r2[[i]], collapse = "+"), sep = "")
    fit_1[[i]] <- miive(model_1[[i]], object, var.cov = T)
  }
  #extract the variable names that had sig sargans in each of the model fit above
  v_list_1 <- lapply(fit_1, badVs)
  #choose the one that has less sig sargan
  #when there's multiple model with the same number of sig sargan
  #choose whoever that's first aka scaling indicator with a bigger r2
  best_1 <- which.min(lapply(v_list_1, length))
  #if best_1 = 2, correlate their errors
  #if best_1 <= 1, keep the model
  #if best_1 > 2, move on to a two factor model
  length_best1 <- length(v_list_1[[best_1]])
  #return(v_list_1[[best_1]])
  if(length_best1 <= 1){
    cat("EFAmiive suggested a one factor model. \n")
    #return(fit_1[[best_1]])
    final_model <- fit_1[[best_1]]
  }
  if(length_best1 >= 2){
    #repeat the above process
    r2_2 <- order_r2(subset(object, select = v_list_1[[best_1]]))
    model_2 <- list()
    fit_2 <- list()
    for (i in 1:length(r2_2)){
      model_2[[i]] <- paste( paste("f1=~", paste(setdiff(r2[[best_1]], r2_2[[1]]), collapse = "+"), sep = ""),
                             paste("f2=~", paste(r2_2[[i]], collapse = "+"), sep = ""), sep = " \n ")
      fit_2[[i]] <- miive(model_2[[i]], object, var.cov = T)
    }
    v_list_2 <- lapply(fit_2, badVs)
    best_2 <- which.min(lapply(v_list_2, length))

    if(length(v_list_2[[best_2]]) <= 1){
      cat("EFAmiive suggested a two factor model. \n")
      final_model <- fit_2[[best_2]]
    }
    if(length(v_list_2[[best_2]]) >= 2){
      #repeat the process again...
      #cat("needs to be furthur broken down to 3 factor model ?")
      r2_3 <- order_r2(subset(object, select = v_list_2[[best_2]]))
      model_3 <- list()
      fit_3 <- list()
      for (i in 1:length(r2_3)){
        model_3[[i]] <- paste( paste("f1=~", paste(setdiff(r2[[best_1]], r2_2[[1]]), collapse = "+"), sep = ""),
                               paste("f2=~", paste(setdiff(r2_2[[best_2]], r2_3[[1]]), collapse = "+"), sep = ""),
                               paste("f3=~", paste(r2_3[[i]], collapse = "+"), sep = ""), sep = " \n ")
        fit_3[[i]] <- miive(model_3[[i]], object, var.cov = T)
      }
      v_list_3 <- lapply(fit_3, badVs)
      best_3 <- which.min(lapply(v_list_3, length))

      if(length(v_list_3[[best_3]]) <= 1){
        cat("EFAmiive suggested a three factor model. \n")
        final_model <- fit_3[[best_3]]
      }
      if(length(v_list_3[[best_3]] >=2)){
        print(model_3[[best_3]])
        print(v_list_3[[best_3]])
        cat("neet to furthur break into a four factor model. \n")
        final_model <- fit_3[[best_3]]
        # print(fit_3[[best_3]])
      }
    }
  }
  return(final_model)
}
