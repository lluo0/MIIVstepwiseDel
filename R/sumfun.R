MIIVdel <- function(model, data){
  #setup the initial step
  setup <- ProbMiivs_c(model, data)
  #create the table for deleted MIIVs for each variables
  #variables that do not have sig sargan will have NA for their MIIVs
  problematicMIIV <- matrix(NA, nrow = 1, ncol = length(setup$probVs))
  colnames(problematicMIIV) <- setup$probVs
  rownames(problematicMIIV) <- "Problematic MIIVs"
  #step 1
  step1 <- step1del_combo(model, data, setup)
  #update the deleted MIIVs table
  for(p in 1:ncol(problematicMIIV)){
    problematicMIIV[,p] <- colnames(step1$drop_MIIV[[p]])
  }
  #run the next step
  stepN <- stepNdel_combo(model, data, setup, step1)
  #update the deleted MIIVs table
  for(p in 1:ncol(problematicMIIV)){
    if(stepN$drop_MIIV[[p]][2,] < 4)
    problematicMIIV[,p] <- colnames(stepN$drop_MIIV[[p]])
  }
  #keep going to the next step until either no degrees of freedom left or newsargan is no longer significant
  while(any(sapply(stepN$drop_MIIV, function(i) i[2,]<4))){
    stepN <- stepNdel_combo(model, data, setup, stepN)
    for(p in 1:ncol(problematicMIIV)){
      if(stepN$drop_MIIV[[p]][2,] < 4)
        problematicMIIV[,p] <- colnames(stepN$drop_MIIV[[p]])
    }
  }


  return(problematicMIIV)
}
