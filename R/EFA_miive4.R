##SETUP FUNCTIONS FOR EFAMIIVE4
#######order_r2 function########
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
#########function to get problematic variables########
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
####step1 function#####
step1_EFAmiive <- function(data, threshold){
  r2 <- order_r2(data)
  models <- list()
  fits <- list()
  for (i in 1:length(r2)){
    models[[i]] <- paste("f1=~", paste(r2[[i]], collapse = "+"), sep = "")
    fits[[i]] <- miive(models[[i]], data, var.cov = T)
  }
  #extract the variable names that had sig sargans in each of the model fit above
  badvar_list <- lapply(fits, function(x) getbadvar(x, threshold))
  #choose the one that has less sig sargan
  #when there's multiple model with the same number of sig sargan
  #choose whoever that's first aka scaling indicator with a bigger r2
  best_num <- which.min(lapply(badvar_list, length))
  #if = 1/0 print the output
  #if problematic variable length (badvar output) > =2 create a new latent factor
  length_best <- length(badvar_list[[best_num]])
  num_badvar <- length_best
  if(num_badvar == 0){
    finalobj <- list(model = models[[best_num]],
                     fit  = fits[[best_num]],
                     num_factor = 1,
                     num_badvar = 0)
  }
  if(!num_badvar ==0){
    #save the problematic variables
    badvar <- badvar_list[[best_num]]
    #save the good variables
    goodvar <- list()
    goodvar[[1]] <- setdiff(r2[[best_num]], badvar)
    #save the part of the model that is good and shoule be untouched for the next step
    goodmodelpart <- paste("f1=~", paste(goodvar[[1]], collapse = "+"), sep = "")
    ##this will become a list once we have multiple factors, so does the goodvar list.
    finalobj <- list(model = models[[best_num]],
                     fit  = fits[[best_num]],
                     num_factor = 1,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart)
  }
  return(finalobj)
}
######stepN function########
stepN_EFAmiive <- function(stepPrev, data, threshold){
  r2 <- order_r2(subset(data, select = stepPrev$badvar))
  models <- list()
  fits <- list()
  num_factor <- stepPrev$num_factor+1
  for (i in 1:length(r2)){
    models[[i]] <- paste(stepPrev$goodmodelpart,
                         paste(paste0("f",num_factor), "=~",paste(r2[[i]], collapse = "+"), sep = ""),
                         sep = "\n")
   fits[[i]] <- miive(models[[i]], data, var.cov = T)
  }
  #extract the variable names that had sig sargans in each of the model fit above
  badvar_list <- lapply(fits, function(x) getbadvar(x, threshold))
  #choose the one that has less sig sargan
  #when there's multiple model with the same number of sig sargan
  #choose whoever that's first aka scaling indicator with a bigger r2
  best_num <- which.min(lapply(badvar_list, length))
  #if = 1/0 print the output
  #if problematic variable length (badvar output) > =2 create a new latent factor
  length_best <- length(badvar_list[[best_num]])
  num_badvar <- length_best
  if(num_badvar==0){
    finalobj <- list(model = models[[best_num]],
                     fit  = fits[[best_num]],
                     num_factor = num_factor,
                     num_badvar = 0)
  }
  if(!num_badvar==0){
    #save the problematic variables
    badvar <- badvar_list[[best_num]]
    #save the good variables - update the ones from the previous output
    goodvar <- stepPrev$goodvar
    #goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)
    goodvar[[num_factor]] <- setdiff(r2[[best_num]], badvar)

    goodmodelpart <- paste(stepPrev$goodmodelpart,
                           paste(paste0("f", num_factor), "=~",
                                 paste(goodvar[[num_factor]], collapse = "+"), sep = ""),
                           sep = "\n")
    finalobj <- list(model = models[[best_num]],
                     fit  = fits[[best_num]],
                     num_factor = num_factor,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart)
  }
  return(finalobj)
}

####when have one single problematic variable####
#function to transfer variables into model syntax
getsyntax <- function(inputlist){
  syntaxlist <- vector()
  for(i in 1:length(inputlist)){
    syntaxlist[i] <- paste(paste0("f",i), "=~",
                           paste(inputlist[[i]], collapse = "+"), sep = "")
  }
  syntax <- paste(syntaxlist, collapse="\n")
  return(syntax)
}

# getsyntax(model_list[[1]])
#the actual function
singlebadvarpruning <- function(stepNobj, data, threshold){
  trys <- stepNobj$num_factor - 1
  model_list <- list()
  for(i in 1:trys){
    model_list[[i]] <- stepNobj$goodvar[1:trys]
  }
  #then attach the badvar to each of the previous factors
  for(i in 1:trys){
    model_list[[i]][[i]] <-  c(model_list[[i]][[i]], stepNobj$badvar)
  }
  models <- lapply(model_list, getsyntax)
  #then attach the last factor to the other factors for the final model syntax
  for(i in 1:trys){
    models[[i]] <- paste(models[[i]],
                         paste(paste0("f", trys+1),"=~",
                               paste(stepNobj$goodvar[[trys+1]], collapse = "+"), sep = ""),
                         sep = "\n")
  }
  #fit the model
  fits <- list()
  for(i in 1:trys){
    fits[[i]] <- miive(models[[i]], data, var.cov = T)
  }
  #see if any fits better
  badvar_list <- lapply(fits, function(x) getbadvar(x, threshold))
  best_num <- which.min(lapply(badvar_list, length))
  num_badvar <- length(badvar_list[[best_num]])
  if(num_badvar == 0){
    finalobj <- list( num_factor = stepNobj$num_factor,
                      model = models[[best_num]],
                      fit = fits[[best_num]])
  }
  if(!num_badvar == 0){
    newmodel <- getsyntax(stepNobj$goodvar)
    newline <- paste(paste(paste0('f',stepNobj$num_factor+1), "=~", stepNobj$badvar, sep = ""),
                     paste(stepNobj$badvar, "~~0*", stepNobj$badvar), sep = "\n")
    newmodel <- paste(newmodel, newline, sep = "\n")
    newfit <- miive(newmodel, data, var.cov = T)
    finalobj <- list(num_factor = stepNobj$num_factor+1,
                     model = newmodel,
                     fit = newfit
                     )
  }
  return(finalobj)
}


#model_list <- singlebadvarpruning(fiveffinal, fivefsim[[1]], .05)




# singleFactor <- function(model, factornum, singlebadvar){
#   newfactor <- paste0("f", factornum+1)
#   newmodel <- paste(model,
#                     paste(newfactor, "=~", singlebadvar),
#                     paste(singlebadvar, "~~ 0*", singlebadvar, sep = ""),
#                     sep = "\n"
#                     )
#   return(newmodel)
# }

#######tests#########

fivefmodel <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3
        f2=~ 1*x4 + .7*x5 + .6*x6
        f3=~ 1*x7 + .8*x8
        f4=~ 1*x9 + .8*x10 + .7*x11
        f5=~ 1*x12 + .8*x13
        f1 ~~ .5*f2
        f1 ~~ .4*f3
        f1~~ .4*f4
        f1~~.4*f5
        f2~~.5*f3
        f2~~.4*f4
        f2~~.4*f5
        f3~~.5*f4
        f3~~.4*f5
        f4~~.4*f5
'
fivefsim <- list()
for (p in 1:30){
  set.seed(123.4+p)
  fivefsim[[p]] <- simulateData(fivefmodel, sample.nobs = 1000)
}
step1 <- step1_EFAmiive(fivefsim[[1]], .05)
step2 <- stepN_EFAmiive(step1, fivefsim[[1]], .05)
step3 <- stepN_EFAmiive(step2, fivefsim[[1]], .05)
step4 <- stepN_EFAmiive(step3, fivefsim[[1]], .05)
step5 <- stepN_EFAmiive(step4, fivefsim[[1]], .05)


onefsim <- list()
for (p in 1:30){
  set.seed(123.4+p)
  onefsim[[p]] <- simulateData('f1=~1*x1+.8*x2+.8*x3+.7*x4+.7*x5', sample.nobs = 1000)
}
step1_EFAmiive(onefsim[[1]], .05)
onef_fit <- miive('f1=~x1+x2+x3+x4+x5', onefsim[[1]], var.cov = T)
getbadvar(onef_fit)

EFAmiive4(sim4[[1]],.05)

step1A <- step1_EFAmiive(sim4[[1]], .05)
step2A <- stepN_EFAmiive(step1A,sim4[[1]], .05)
singlebadvarpruning(step2A, sim4[[1]], .05)

sm4b <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 + .65*x5 + .6*x6 + .6*x7 + .55*x8
          x4 ~~ .5*x5
          x4 ~~ .4*x2
          x2 ~~ .5*x5'
sim4b <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim4b[[p]] <- simulateData(sm4b, sample.nobs = 1000)
}
step1A <- step1_EFAmiive(sim4b[[1]], .05)
step2A <- stepN_EFAmiive(step1A,sim4b[[1]], .05)
singlebadvarpruning(step2A, sim4b[[1]], .05)
EFAmiive4(sim4b[[1]], .05)

sm4c <- 'f1 =~ 1*x1 +.7*x3 + .6*x6 + .6*x7 + .55*x8
          f2=~ 1*x2 + .8*x4 + .7*x5
f1~~.4*f2'
sim4c <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim4c[[p]] <- simulateData(sm4c, sample.nobs = 1000)
}
EFAmiive4(sim4c[[1]], .05)
