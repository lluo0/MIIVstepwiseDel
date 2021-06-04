stepNdel_c <- function(model, data, setup, stepPrevious_final){
  probVs <- setup$probVs
  badmiivs <- setup$badmiivs
  goodmiivs <- setup$goodmiivs
  badmiivs_int <- setup$badmiivs_int
  sig_sargantable_int <- setup$sig_sargantable_int

  # badmiivs_new <- stepPrevious_final$badmiivs_nextstep
  # ##first set up the new combos of miivs to be deleted/retained
  # for (p in 1:length(badmiivs_int)){
  #   new_colnum <- length(badmiivs_new[[p]])
  #   new_rownum <- length(badmiivs_int[[p]])-length(badmiivs_new[[p]])+1
  #   badmiivs[[p]] <- matrix(NA, nrow = new_rownum, ncol = new_colnum)
  #   for (i in 1:new_rownum){
  #     badmiivs[[p]][i,] <- setdiff(badmiivs_int[[p]],badmiivs_new[[p]])[i]
  #   }
  #   for (q in 1:new_colnum){
  #     badmiivs[[p]][new_rownum,q] <- badmiivs_new[[p]][q]
  #   }
  # }
  #

  ##add more candidates
  badmiivs_new <- stepPrevious_final$badmiivs_nextstep
  newbadmiivs <- list()
  for(p in 1:length(badmiivs_new)){
    newbadmiivs[[p]] <- list()
    new_colnum <- length(badmiivs_new[[1]][[1]])
    new_rownum <- dim(stepPrevious_final$badmiivs[[1]])[1]+1
    for(i in 1:length(badmiivs_new[[p]])){
      newbadmiivs[[p]][[i]] <- matrix(NA, nrow = new_rownum, ncol = new_colnum)
      for(r in 1:new_rownum){
        newbadmiivs[[p]][[i]][r,] <- setdiff(badmiivs_int[[p]], badmiivs_new[[p]][[i]])[r]
      }
      for(q in 1:new_colnum){
        newbadmiivs[[p]][[i]][new_rownum, q] <- badmiivs_new[[p]][[i]][q]
      }
    }

    # badmiivs_new[[p]] <- matrix(NA, nrow = new_rownum, ncol = new_colnum)
  }
  names(newbadmiivs) <- names(badmiivs)

  for(p in 1:length(badmiivs)){
    badmiivs[[p]] <- t(unique(t(apply(t(do.call(cbind, newbadmiivs[[p]])), 1, sort))))
  }

  ##then repeat the steps in step1del_combo
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


  drop_MIIV <- lapply(delMiivsSargan, function(i) i[which(i[2,]<4)])
  #drop_MIIV <- drop_MIIV[sapply(drop_MIIV, function(i) !ncol(i)==0)]

  # maxdiff <- sapply(delMiivsSargan, function(i) max(i[1,]))
  # drop_location <- names(delMiivsSargan)[sapply(seq_along(delMiivsSargan),
  #                                                function(x) {any(delMiivsSargan[[x]]==maxdiff[x])})]
  # drop_location_num <- sapply(drop_location, function(x) match(x, names(delMiivsSargan)))
  # # drop_MIIV_2 <- sapply(drop_location_num, function(i)
  # #    colnames(delMiivsSargan[[i]])[delMiivsSargan[[i]][1,] == maxdiff[i]])
  # drop_MIIV_2 <- sapply(drop_location_num, function(i)
  #   delMiivsSargan[[i]][delMiivsSargan[[i]][1,] == maxdiff[i]])

  minsargan <- sapply(delMiivsSargan, function(i) min(i[2,]))

  # drop_num <- vector()
  # for(p in 1:length(minsargan)){
  #   drop_num[p] <- which(delMiivsSargan[[p]][2,]==minsargan[p])
  # }
  drop_MIIV_2 <- list()
  for(p in 1:length(minsargan)){
    drop_MIIV_2[[p]] <- delMiivsSargan[[p]][which(delMiivsSargan[[p]][2,]==minsargan[p])]
  }

  droptest <- list()
  for(p in 1:length(drop_MIIV)){
    if(!ncol(drop_MIIV[[p]])==0){
      droptest[[p]] <-drop_MIIV[[p]]
    }
    else{
      droptest[[p]] <- drop_MIIV_2[[p]]
    }
    names(droptest)[p] <- names(drop_MIIV)[p]
  }

  badmiivs_nextstep <- list()
  for(p in 1:length(droptest)){
    badmiivs_nextstep[[p]] <- list()
    names(badmiivs_nextstep)[p] <- names(droptest)[p]
    for(i in 1:length(droptest[[p]])){
      badmiivs_nextstep[[p]][[i]] <- strsplit(colnames(droptest[[p]])[i], split = "+", fixed = T)[[1]]
    }
  }

  finalobj <- list(#maxdiff = maxdiff,
    #drop_location = drop_location,
    #drop_location_num = drop_location_num,
    drop_MIIV = droptest,
    delMiivsSargan = delMiivsSargan,
    badmiivs = badmiivs,
    badmiivs_nextstep = badmiivs_nextstep
    #badmiivs_all = badmiivs_candidate,
    #sargan_eqnorder = sargan_eqnorder
    #badmiivs = badmiivs_2,
    #delMiivsSargan_nextstep = delMiivsSargan_nextstep
  )
  return(finalobj)

}
