step1del_c <- function(model, data, setup){
  probVs <- setup$probVs
  badmiivs <- setup$badmiivs
  goodmiivs <- setup$goodmiivs
  badmiivs_int <- setup$badmiivs_int
  sig_sargantable_int <- setup$sig_sargantable_int

  fit_del1 <- list()
  newmiivs <- list()
  delMiivsSargan <- list()
  for (p in 1:length(probVs)){
    ##p is the name of the problematic variable
    fit_del1[[p]] <- list()
    newmiivs[[p]] <- list()
    delMiivsSargan[[p]] <- as.data.frame(matrix(NA, nrow = 2, ncol = ncol(badmiivs[[p]])) )
    names(delMiivsSargan)[p] <- probVs[p]
    rownames(delMiivsSargan[[p]]) <- c("sargandef", "newsargan")
    for (i in 1:ncol(badmiivs[[p]])){
      ##i for each deleted MIIV
      ##then get the list of MIIVs
      newmiivs[[p]][[i]] <- paste0(goodmiivs[[p]], sep = "\n",
                                   paste0(names(badmiivs)[p], sep = "~",
                                          paste0(paste0(badmiivs[[p]][,i]), collapse = "+")))
      ##get the fit for each MIIV deletion
      fit_del1[[p]][[i]] <- miive(model = model, data = data, var.cov = T, miiv.check = F,
                                  instruments = newmiivs[[p]][[i]])
      colnames(delMiivsSargan[[p]])[i] <- paste0(setdiff(badmiivs_int[[p]], paste0(badmiivs[[p]][,i])), collapse = "+")
      #paste0(paste0(badmiivs[[p]][,i]), collapse = "+")
      #paste0(setdiff(badmiivs_int[[p]], paste0(badmiivs[[p]][,i])), collapse = "+")

      ##then get the sargan for each miiv deletion
      delMiivsSargan[[p]][1,i] <- sig_sargantable_int[1,p] - fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
      delMiivsSargan[[p]][2,i] <- fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
    }
  }
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
  # ##get the equation order for these variables (the fit$eqn[[p]])
  # sargan_eqnorder <- as.numeric(sig_sargantable_int[3,])
  #
  # badmiivs_candidate <- list()
  # for (p in 1:length(drop_MIIV)){
  #   badmiivs_candidate[[p]] <- strsplit(drop_MIIV[p],split ='+', fixed = T)[[1]]
  #   names(badmiivs_candidate)[p] <- probVs[p]
  # }
  #
  drop_MIIV <- lapply(delMiivsSargan, function(i) i[which(i[2,]<4)])

  badmiivs_nextstep <- list()
  for(p in 1:length(drop_MIIV)){
    badmiivs_nextstep[[p]] <- list()
    names(badmiivs_nextstep)[p] <- names(drop_MIIV)[p]
    for(i in 1:length(drop_MIIV[[p]])){
      badmiivs_nextstep[[p]][[i]] <- strsplit(colnames(drop_MIIV[[p]])[i], split = "+", fixed = T)[[1]]
    }
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
