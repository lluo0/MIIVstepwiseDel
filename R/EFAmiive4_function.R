library(MIIVsem)
library(lavaan)

########EFAmiive4#######
##ignore single problematic variable
EFAmiive4 <- function(data, threshold){
  step1 <- step1_EFAmiive(data, threshold)
  if(step1$num_badvar <= 1){
    int_obj <- step1
  }
  if(step1$num_badvar >= 2){
    stepN <- stepN_EFAmiive(step1,data, threshold)
    int_obj <- stepN
    while(stepN$num_badvar >= 2){
      stepN <- stepN_EFAmiive(stepN, data, threshold)
      int_obj <- stepN
    }
  }
  return(int_obj)
}
##deal with single problematic variable##
EFAmiive4 <- function(data, threshold){
  step1 <- step1_EFAmiive(data, threshold)
  if(step1$num_badvar <= 1){
    int_obj <- step1
  }
  if(step1$num_badvar >= 2){
    stepN <- stepN_EFAmiive(step1,data, threshold)
    int_obj <- stepN
    while(stepN$num_badvar >= 2){
      stepN <- stepN_EFAmiive(stepN, data, threshold)
      int_obj <- stepN
    }
  }
  # ##this moves any single variable factor to the badvar before moves on to the pruning stage
  # if(any(lengths(int_obj$goodvar)==1)){
  # }
  ##pruning stage for single problematic variable
  if(int_obj$num_badvar == 1 && length(int_obj$goodvar[[length(int_obj$goodvar)]])>1){
    finalobj <- singlebadvarpruning(int_obj, data, threshold)
  }else{
    finalobj <- list(num_factor = int_obj$num_factor,
                     model = int_obj$model,
                     fit = int_obj$fit)
  }
  return(finalobj)
}

########tests########
EFAmiive4(onefsim[[1]], .05)
fiveffinal <- EFAmiive4(fivefsim[[1]], .05)
