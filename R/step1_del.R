step1del <- function(model, data, setup){
  probVs <- setup$probVs
  badmiivs <- setup$badmiivs
  goodmiivs <- setup$goodmiivs
  #badmiivs_int <- setup$badmiivs_int
  sig_sargantable_int <- setup$sig_sargantable_int

  fit_del1 <- list()
  newmiivs <- list()
  delMiivsSargan <- list()
  for (p in 1:length(probVs)){
    ##p is the name of the problematic variable
    fit_del1[[p]] <- list()
    newmiivs[[p]] <- list()
    delMiivsSargan[[p]] <- as.data.frame(matrix(NA, nrow = 3, ncol = length(badmiivs[[p]])) )
    names(delMiivsSargan)[p] <- probVs[p]
    rownames(delMiivsSargan[[p]]) <- c("sargandef", "newsargan", "newsargan.p")
    for (i in 1:length(badmiivs[[p]])){
      ##i for each deleted MIIV
      ##then get the list of MIIVs
      newmiivs[[p]][[i]] <- paste0(goodmiivs[[p]], sep = "\n",
                                   paste0(names(badmiivs)[p], sep = "~",
                                          paste0(badmiivs[[p]][[i]], collapse = "+")))
      ##get the fit for each MIIV deletion
      fit_del1[[p]][[i]] <- miive(model = model, data = data, var.cov = T, miiv.check = F,
                                  instruments = newmiivs[[p]][[i]])
      colnames(delMiivsSargan[[p]])[i] <- names(badmiivs[[p]])[i]

      ##then get the sargan for each miiv deletion
      delMiivsSargan[[p]][1,i] <- sig_sargantable_int[1,p] - fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
      delMiivsSargan[[p]][2,i] <- fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
      delMiivsSargan[[p]][3,i] <- fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan.p
    }
  }

  drop_MIIV <- lapply(delMiivsSargan, function(i) i[which(i[2,]==min(i[2,]))])

  badmiivs_nextstep <- list()
  for(p in 1:length(drop_MIIV)){
    badmiivs_nextstep[[p]] <- list()
    names(badmiivs_nextstep)[p] <- names(drop_MIIV)[p]
    # for(i in 1:length(drop_MIIV[[p]])){
    #   badmiivs_nextstep[[p]][[i]] <- strsplit(colnames(drop_MIIV[[p]])[i], split = "+", fixed = T)[[1]]
    # }
    badmiivs_nextstep[[p]] <- badmiivs[[p]][which(names(badmiivs[[p]])==colnames(drop_MIIV[[p]]))]
  }

  finalobj <- list(#maxdiff = maxdiff,
    #drop_location = drop_location,
    #drop_location_num = drop_location_num,
    drop_MIIV = drop_MIIV,
    delMiivsSargan = delMiivsSargan,
    badmiivs = badmiivs,
    badmiivs_nextstep = badmiivs_nextstep
    #sargan_eqnorder = sargan_eqnorder
    #badmiivs = badmiivs_2,
    #delMiivsSargan_nextstep = delMiivsSargan_nextstep
  )
  return(finalobj)
}
