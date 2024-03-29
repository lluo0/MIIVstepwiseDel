# MIIVdel_new <- function(model, data, threshold=.05){
#   fit <- miive(model, data, var.cov = T)
#   sarganlist <- vector()
#   for(p in 1:length(fit$eqn)){
#     sarganlist <- c(sarganlist, fit$eqn[[p]]$sargan.p)
#   }
#     if(!any(sarganlist <= threshold))
#     {#cat('MIIVdel STOP: no significant sargan found. \n')
#       problematicMIIV <- NA
#   }
#   else
# {  #setup the initial step
#   setup <- ProbMiivs(model, data, threshold)
#   #create the table for deleted MIIVs for each variables
#   #variables that do not have sig sargan will have NA for their MIIVs
#   problematicMIIV <- matrix(NA, nrow = 1, ncol = length(setup$probVs))
#   colnames(problematicMIIV) <- setup$probVs
#   rownames(problematicMIIV) <- "Problematic MIIVs"
#   #a list that stores all sargan p values for the problematic variables one ach step of diff deleted MIIVs
#   sarganplist <- list()
#   #step 1
#   step1 <- step1del(model, data, setup)
#   #update the deleted MIIVs table accordingly
#   for(p in 1:ncol(problematicMIIV)){
#     if(step1$drop_MIIV[[p]][3,] > threshold){
#       problematicMIIV[,p] <- colnames(step1$drop_MIIV[[p]])
#     }
#     sarganplist[[p]] <- step1$drop_MIIV[[p]][3,]
#   }
#
#
#   ##then decide if we need to go to the next step
#   if(any(sapply(step1$drop_MIIV, function(i) i[3,]<threshold))){
#     #run the next step
#     stepN <- stepNdel(model, data, setup, step1)
#     #update the deleted MIIVs table
#     for(p in 1:ncol(problematicMIIV)){
#       #check the previous sargan p values, if non is > .05 then update
#       #this way we are not including extra MIIVs
#       #eg., for x2, the only problematic MIIV is x3. but in step 2, when we delete x3+y for x2, likely all sargan p will be > .05
#       #and we don't want to update the deleted MIIIV for x2 because it was already nonsignificant at a previous step
#      # if(!any(sarganplist[[p]]>.05) &&  stepN$drop_MIIV[[p]][3,] >= .05){
#       # changed to check only the max value instead of an any loop so its faster 7.16.20
#         if(max(sarganplist[[p]]<threshold) &&  !stepN$drop_MIIV[[p]][3,] < threshold){
#         problematicMIIV[,p] <- colnames(stepN$drop_MIIV[[p]])
#         }
#       sarganplist[[p]] <- c(sarganplist[[p]], stepN$drop_MIIV[[p]][3,])
#     }
#     #keep going to the next step until either no degrees of freedom left or newsargan is no longer significant
#     while(any(sapply(stepN$drop_MIIV, function(i) i[3,]< threshold))){
#       stepN <- stepNdel(model, data, setup, stepN)
#       for(p in 1:ncol(problematicMIIV)){
#         if(max(sarganplist[[p]]< threshold) &&  !stepN$drop_MIIV[[p]][3,] < threshold){
#           problematicMIIV[,p] <- colnames(stepN$drop_MIIV[[p]])
#         }
#         sarganplist[[p]] <- c(sarganplist[[p]], stepN$drop_MIIV[[p]][3,])
#       }
#     }
#   }
#   }
#   return(problematicMIIV)
# }

MIIVdel_new <- function(model, data, threshold=.05){
  fit <- miive(model, data, var.cov = T)
  sarganlist <- vector()
  for(p in 1:length(fit$eqn)){
    sarganlist <- c(sarganlist, fit$eqn[[p]]$sargan.p)
  }
  if(!any(sarganlist <= threshold))
  {#cat('MIIVdel STOP: no significant sargan found. \n')
    problematicMIIV <- NA
  }
  else
  {  #setup the initial step
    setup <- ProbMiivs(model, data)
    #create the table for deleted MIIVs for each variables
    #variables that do not have sig sargan will have NA for their MIIVs
    problematicMIIV <- matrix(NA, nrow = 1, ncol = length(setup$probVs))
    colnames(problematicMIIV) <- setup$probVs
    rownames(problematicMIIV) <- "Problematic MIIVs"
    #a list that stores all sargan p values for the problematic variables one each step of diff deleted MIIVs
    sarganplist <- list()
    #step 1
    step1 <- step1del(model, data, setup)
    #update the deleted MIIVs table accordingly
    for(p in 1:ncol(problematicMIIV)){
      if(step1$drop_MIIV[[p]][3,] > threshold){
        problematicMIIV[,p] <- colnames(step1$drop_MIIV[[p]])
      }
      sarganplist[[p]] <- step1$drop_MIIV[[p]][3,]
    }


    ##then decide if we need to go to the next step
    if(any(sapply(step1$drop_MIIV, function(i) i[3,]<threshold))){
      #run the next step
      stepN <- stepNdel(model, data, setup, step1)
      #update the deleted MIIVs table
      for(p in 1:ncol(problematicMIIV)){
        #check the previous sargan p values, if non is > .05 then update
        #this way we are not including extra MIIVs
        #eg., for x2, the only problematic MIIV is x3. but in step 2, when we delete x3+y for x2, likely all sargan p will be > .05
        #and we don't want to update the deleted MIIIV for x2 because it was already nonsignificant at a previous step
        # if(!any(sarganplist[[p]]>.05) &&  stepN$drop_MIIV[[p]][3,] >= .05){
        # changed to check only the max value instead of an any loop so its faster 7.16.20
        if(max(sarganplist[[p]]<threshold) &&  !stepN$drop_MIIV[[p]][3,] < threshold){
          problematicMIIV[,p] <- colnames(stepN$drop_MIIV[[p]])
        }
        sarganplist[[p]] <- c(sarganplist[[p]], stepN$drop_MIIV[[p]][3,])
      }

      #this makes sure we don't run out of df for sargan test
      nextlength <- vector()
      for(p in 1:length(stepN$badmiivs_nextstep)){
        nextlength[p] <- length(stepN$badmiivs_nextstep[[p]][[1]])
      }
      #keep going to the next step until either no degrees of freedom left or newsargan is no longer significant
      while(any(sapply(stepN$drop_MIIV, function(i) i[3,]< threshold)) && any(nextlength>2)){
        stepN <- stepNdel(model, data, setup, stepN)
        for(p in 1:ncol(problematicMIIV)){
          if(max(sarganplist[[p]]< threshold) &&  !stepN$drop_MIIV[[p]][3,] < threshold){
            problematicMIIV[,p] <- colnames(stepN$drop_MIIV[[p]])
          }
          sarganplist[[p]] <- c(sarganplist[[p]], stepN$drop_MIIV[[p]][3,])
          for(p in 1:length(stepN$badmiivs_nextstep)){
            nextlength[p] <- length(stepN$badmiivs_nextstep[[p]][[1]])
          }
        }
      }
      #add back the delMIIVs for nas in problematicMIIV
      for(p in 1:ncol(problematicMIIV)){
        if(is.na(problematicMIIV[1,p])){
          problematicMIIV[1,p] <- colnames(stepN$drop_MIIV[[p]])
        }
      }
    }
  }
  return(problematicMIIV)
}

