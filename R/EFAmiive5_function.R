EFAmiive5 <- function(data, threshold = .05, criterion = 1){
  step1 <- step1_E5(data, threshold)
  finalobj <- step1
  if(!length(step1$badvar)==0){
    step2 <- step2_E5(step1, data, threshold, criterion)
    finalobj <- step2
    if(!length(step2$badvar)==0){
      stepN <- stepN_E5(step2, data, threshold, criterion)
      finalobj <- stepN
      while(stepN$nextstep == 'yes'){
        stepN <- stepN_E5(stepN, data, threshold, criterion)
        finalobj <- stepN
      }
    }
  }
    return(finalobj[1:4])
}

