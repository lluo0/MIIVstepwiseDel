

##3 algorithms

##########alg 1##########
EFAmiive <- function(object){
  ##check the r2 for each variable and use the highest r2 as the initial scaling indicator
  r2 <- matrix(NA, nrow = 1, ncol = dim(object)[2])
  colnames(r2) <- colnames(object)
  for (i in 1:dim(object)[2]){
    r2[,i] <- summary(lm(paste(colnames(object)[i], paste(colnames(object)[-i], collapse = "+"), sep = "~"), data = object))$r.squared
  }
  r2 <- as.matrix(t(r2[,order(r2[nrow(r2),],decreasing=TRUE)]))
  model_1 <- paste("f1=~", paste(colnames(r2), collapse = "+"), sep = "")
  fit_1 <- miive(model_1, object, var.cov = T)
  ##store sig sargans into one vector
  v_list_1 <- vector()
  for (p in 1:length(fit_1$eqn))
    if (fit_1$eqn[[p]]$sargan.p < .05){
      v_list_1 <- append(v_list_1, fit_1$eqn[[p]]$DVobs)
    }

  ##determine the next step
  if(length(v_list_1) <= 1){
    cat("EFAmiive suggested a one factor model. \n \n")
    summary(fit_1)} #print the fit if none or only one sig sargan
  if(length(setdiff(colnames(r2), v_list_1)) == 1 && !length(v_list_1) == 1){
    cat("temp_error 1. \n All sargan on the factor are significant. \n")
    summary(fit_1)
    stop
  }
  ##if more than 1 sig sargan, load the sig ones on a second factor
  ##repeat the r2 calculation stage
  if(length(v_list_1) > 1 && !length(setdiff(colnames(r2), v_list_1)) == 1){
    ##repeat the r2 calculation process
    r2_1 <- matrix(NA, nrow = 1, ncol = length(v_list_1))
    colnames(r2_1) <- paste(v_list_1, sep = ",")
    for (i in 1:length(v_list_1)){
      r2_1[,i] <- summary(lm(paste(v_list_1[i], paste(v_list_1[-i], collapse = "+"), sep = "~"), data = object))$r.squared
    }
    r2_1 <- as.matrix(t(r2_1[,order(r2_1[nrow(r2_1),],decreasing=TRUE)]))

    model_2 <- paste(paste("f1=~", paste(setdiff(colnames(r2), colnames(r2_1)), collapse = "+"), sep = ""),
                     paste("f2=~", paste(colnames(r2_1), collapse = "+"), sep = ""), sep = " \n ")
    fit_2 <- miive(model_2, object, var.cov = T)
    v_list_2 <- vector()
    for (p in 1:length(fit_2$eqn))
      if (fit_2$eqn[[p]]$sargan.p < .05){
        v_list_2 <- append(v_list_2, fit_2$eqn[[p]]$DVobs)
      }
    ##repeat the process
    if(length(v_list_2) <= 1){
      cat("EFAmiive suggested a two factor model. \n \n")
      summary(fit_2)}
    if(length(setdiff(colnames(r2_1), v_list_2)) == 1 && !length(v_list_2) == 1){
      cat("temp_error 2. \n All sargan on the second factor are significant. \n")
      summary(fit_2)
      # model_2a <- list()
      # while (p <= 1:length(colnames(r2_1))-1){
      # model_2a[[p]] <- paste(paste("f1=~", paste(setdiff(colnames(r2), colnames(r2_1)), collapse = "+"), sep = ""),
      #          paste("f2=~", paste(c(colnames(r2_1)[-p], colnames(r2_1)[p]), collapse = "+"), sep = ""), sep = " \n ")
      # p+1
      # }
      # print(model_2a)
      # print(paste(paste("f1=~", paste(setdiff(colnames(r2), colnames(r2_1)), collapse = "+"), sep = ""),
      #            paste("f2=~", paste(c(colnames(r2_1)[-1], colnames(r2_1)[1]), collapse = "+"), sep = ""), sep = " \n "))
    }
    if(length(v_list_2) > 1 && !length(setdiff(colnames(r2_1), v_list_2)) == 1){
      r2_2 <- matrix(NA, nrow = 1, ncol = length(v_list_2))
      colnames(r2_2) <- paste(v_list_2, sep = ",")
      for (i in 1:length(v_list_2)){
        r2_2[,i] <- summary(lm(paste(v_list_2[i], paste(v_list_2[-i], collapse = "+"), sep = "~"), data = object))$r.squared
      }
      r2_2 <- as.matrix(t(r2_2[,order(r2_2[nrow(r2_2),],decreasing=TRUE)]))
      model_3 <- paste(paste("f1=~", paste(setdiff(colnames(r2), colnames(r2_1)), collapse = "+"), sep = ""),
                       paste("f2=~", paste(setdiff(colnames(r2_1), colnames(r2_2)), collapse = "+"), sep = ""),
                       paste("f3=~", paste(colnames(r2_2), collapse = "+"), sep = ""), sep = " \n ")
      fit_3 <- miive(model_3, object, var.cov = T)
      cat("EFAmiive suggested a three factor model. \n \n")
      summary(fit_3)
    }
  }
}

########alg 2##########
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
  #return(v_list)
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


###########alg 3############

############demonstration examples##############
#m3: 2 factor model w/ x7 crossload on f1
sm3 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2'
sim3 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim3[[p]] <- simulateData(sm3, sample.nobs = 1000)
}
