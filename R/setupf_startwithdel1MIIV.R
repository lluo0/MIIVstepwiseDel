ProbMiivs <- function(model, data){
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
  badmiivs <- lapply(probVs, function(i) getbadmiivs(i,fit_int))
  names(badmiivs) <- probVs


  # fit_del1 <- list()
  # newmiivs <- list()
  # delMiivsSargan <- list()
  # for (p in 1:length(probVs)){
  #   ##p is the name of the problematic variable
  #  fit_del1[[p]] <- list()
  #  newmiivs[[p]] <- list()
  #  delMiivsSargan[[p]] <- as.data.frame(matrix(NA, nrow = 2, ncol = length(badmiivs[[p]])) )
  #  names(delMiivsSargan)[p] <- probVs[p]
  #  rownames(delMiivsSargan[[p]]) <- c("sargandef", "newsargan")
  #  for (i in 1:length(badmiivs[[p]])){
  #    ##i for each deleted MIIV
  #    ##then get the list of MIIVs
  #    newmiivs[[p]][[i]] <- paste0(goodmiivs[[p]], sep = "\n",
  #                             paste0(names(badmiivs)[p], sep = "~",
  #                                    paste0(badmiivs[[p]][[i]], collapse = "+")))
  #    ##get the fit for each MIIV deletion
  #    fit_del1[[p]][[i]] <- miive(model = model, data = data, var.cov = T, miiv.check = F,
  #                                instruments = newmiivs[[p]][[i]])
  #    colnames(delMiivsSargan[[p]])[i] <- names(badmiivs[[p]])[i]
  #
  #    ##then get the sargan for each miiv deletion
  #    delMiivsSargan[[p]][1,i] <- sig_sargantable_int[1,p] - fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
  #    delMiivsSargan[[p]][2,i] <- fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
  #  }
  # }
  # df1 <- qchisq(.95, df = 1)
  # #store the max chi square drop
  # maxdiff <- sapply(delMiivsSargan, function(i) max(i[1,]))
  # #max needs to be greater than df1
  # maxdiff <- maxdiff[maxdiff>df1]
  # #find list name
  # drop_location <- names(delMiivsSargan)[sapply(seq_along(delMiivsSargan),
  #                                               function(x) {any(delMiivsSargan[[x]]==maxdiff[x])})]
  # drop_location_num <- sapply(drop_location, function(x) match(x, names(delMiivsSargan)))
  #
  # drop_MIIV <- sapply(drop_location_num, function(i)
  #   colnames(delMiivsSargan[[i]])[delMiivsSargan[[i]][1,] == maxdiff[i]])
  # step1obj <- step1del(model = "f1=~x1+x2+x3+x4+x5+x6+x7+x8",
  #                      data = simD3[[1]],
  #                      probVs, badmiivs, goodmiivs, sig_sargantable_int)

  tentobj <- list(probVs = probVs,
                  badmiivs = badmiivs,
                  goodmiivs = goodmiivs,
                  sig_sargantable_int = sig_sargantable_int)
  return(tentobj)
}


getbadmiivs <- function(single_probVs, fit){
  miivs_all_nosign <- getMIIVs_noplussign(fit)
  badmiivs <- list()
  single_probV_miivs <- miivs_all_nosign[which(names(miivs_all_nosign)==single_probVs)][[1]]
  for (p in 1:length(single_probV_miivs)){
    badmiivs[[p]] <- single_probV_miivs[-p]
    names(badmiivs)[p] <- single_probV_miivs[p]
  }
  return(badmiivs)
}
