

###to get the model based on the scaling indicator and number of factor
getcfamodel <- function(num_factor, scalingindicator, nonscaling){
  modelparts <- list()
  for(i in 1:num_factor){
    modelparts[[i]] <- paste0('f', i, '=~', scalingindicator[i], '+', paste0(nonscaling, collapse = '+'))
  }
  model <- paste0(modelparts, collapse ='\n')
  return(model)
}




######using EFAmiive in a more confirmatory fashion######

CFAmiive_step1 <- function(data, num_factor, criterion, threshold = .05){ #criterion is the criterion in selecting the best scaling indicators
  #prep
  allvar <- colnames(data)
  allscaling <- combn(allvar, num_factor, FUN = NULL, simplify = F) #a list of all possible combinations of scaling indicators
  allnonscaling <- lapply(allscaling, function(x) setdiff(allvar,x))

  allmodels <- allfit <- list()
  #create all possible models by load nonscaling variables on all factors and estimate these models
  for(p in 1:length(allscaling)){
    allmodels[[p]] <- getcfamodel(num_factor, allscaling[[p]], allnonscaling[[p]])
    allfit[[p]] <- miive(allmodels[[p]], data, var.cov = T)
  }

  #choose the 'best' model
  allbadvar <- list()
  if(criterion == 'sargan+factorloading'){
    for(p in 1:length(allscaling)){
      allbadvar[[p]] <- getbadvar(allfit[[p]], threshold)
    }
  }

  return(allfit)
}


#######testing#########

#m2: 2 factor model
sm2 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
      f1 ~~ .5*f2'
sim2 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim2[[p]] <- simulateData(sm2, sample.nobs = 1000)
}
data <- sim2[[1]]
num_factor <- 2


combn(allvar, 2, FUN = NULL, simplify = F)
