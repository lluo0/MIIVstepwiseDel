ProbMiivs_c <- function(model, data){
  ##first run the specified model and get variables that have significant sargan
  fit_int <- miive(model = model, data = data, var.cov = T)
  sargantable_int <- getSarganTable(fit_int)
  #probVs <- sigSargan(fit_int)
  sig_sargantable_int <- sargantable_int[which(sargantable_int[2,]<.05)]
  probVs <- colnames(sig_sargantable_int)

  ##get MIIVs for each variables
  miivs_all <- getMIIVs(fit_int)

  ##then step delete
  ##first get the good miivs that do not need to change for each problematic Vs
  goodmiivs <- getgoodmiivs(probVs, miivs_all)
  ##then get the bad miivs that need to be changed
  badmiivs <- lapply(probVs, function(i) getbadmiivs_combo(i,fit_int))
  names(badmiivs) <- probVs

  badmiivs_int <- getbadmiivs_all(probVs, fit_int)

  tentobj <- list(probVs = probVs,
                  badmiivs = badmiivs,
                  badmiivs_int = badmiivs_int,
                  goodmiivs = goodmiivs,
                  sig_sargantable_int = sig_sargantable_int)
  return(tentobj)
}
