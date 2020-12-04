sugmod <- function(model, data, MIIVdelobj){
  fit <- miive(model, data, var.cov = T)
  allmiivs <- getMIIVs_noplussign(fit)

  if(!is.na(MIIVdelobj)){
    for(p in 1:length(allmiivs)){
      for(q in 1:dim(MIIVdelobj)[2]){
        if(names(allmiivs)[p]==colnames(MIIVdelobj)[q]){
          allmiivs[[p]] <- setdiff(allmiivs[[p]], MIIVdelobj[,q])
        }
      }
    }

    for (p in 1:length(allmiivs)){
      allmiivs[[p]] <- paste0(names(allmiivs)[p], sep="~",
                              paste0(allmiivs[[p]], collapse = "+"))
    }

    allmiivs_mod <- paste0(allmiivs,collapse = "\n")

    newfit <- miive(model, data, var.cov = T, miiv.check = F, instruments = allmiivs_mod)

    sugmodobj <- list(newMIIVs = allmiivs_mod,
                      fit = newfit)
  }else{
    cat('the original model is recommended. \n')
    sugmodobj <- NA
  }
  return(sugmodobj)
}

sugmod(model1, sim1[[1]], sim1all[[1]])

MIIVdel_new(model1, sim1[[1]], .05)

