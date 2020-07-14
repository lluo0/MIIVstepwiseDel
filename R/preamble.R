library(lavaan)
library(MIIVsem)

##function that finds the variables with significant sargan
sigSargan <- function(fit){
  v_list <- vector()
  for (p in 1:length(fit$eqn))
    if (fit$eqn[[p]]$sargan.p < .05){
      v_list <- append(v_list,fit$eqn[[p]]$DVobs)
    }
  return(v_list)
}

##get the initial sargian values
getSarganTable <- function(fit){
  sargantable <- matrix(NA, nrow = 3, ncol = length(fit$eqn))
  sargantable <- as.data.frame(sargantable)
  for (p in 1:length(fit$eqn)){
    sargantable[1,p] <- fit$eqn[[p]]$sargan
    sargantable[2,p] <- fit$eqn[[p]]$sargan.p
    sargantable[3,p] <- p
    colnames(sargantable)[p] <- fit$eqn[[p]]$DVobs
  }
  rownames(sargantable) <- c("sargan", "sargan.p","eqnorder")
  return(sargantable)
}

##function that gets the MIIVs for all variables. class is list
getMIIVs <- function(fit){
  miivslist <- list()
  for (p in 1:length(fit$eqn)){
    miivslist[[p]] <- paste0(fit$eqn[[p]]$MIIVs, collapse = "+")
    names(miivslist)[p] <- as.character(fit$eqn[[p]]$DVobs)
  }
  return(miivslist)
}

getMIIVs_noplussign <- function(fit){
  miivslist <- list()
  for (p in 1:length(fit$eqn)){
    miivslist[[p]] <- fit$eqn[[p]]$MIIVs
    names(miivslist)[p] <- as.character(fit$eqn[[p]]$DVobs)
  }
  return(miivslist)
}

##functon that combines the list from "getMIIVs" function to a specification that can miivsem can read
combineMiivs <- function(miivslist){
  printmiivs <- list()
  for (p in 1:length(miivslist)){
    printmiivs[[p]] <- paste0(names(miivslist)[p], sep="~",
                              paste0(miivslist[[p]], collapse = "+"))
  }
  printmiivs_all <- paste0(printmiivs,collapse = "\n")
  return(printmiivs_all)
}

##function that creates the miivs for problems variables that do not need to change
##e.g., for x1-8, if x5 is the prob variabls, this includes all the miivs for other variables besides x5
getgoodmiivs <- function(probVs, miivs_all){
  goodmiivs <- list()
  othermiivs <- list()
  for (p in 1:length(probVs)){
    othermiivs[[p]] <- miivs_all[grepl(probVs[p], unlist(miivs_all))]
    goodmiivs[[p]] <- combineMiivs(othermiivs[[p]])
    names(goodmiivs)[p] <- probVs[p]
  }
  return(goodmiivs)
}


getbadmiivs_combo <- function(single_probVs, fit){
  miivs_all_nosign <- getMIIVs_noplussign(fit)
  # badmiivs <- list()
  single_probV_miivs <- miivs_all_nosign[which(names(miivs_all_nosign)==single_probVs)][[1]]
  # for (p in 1:length(single_probV_miivs)){
  #   badmiivs[[p]] <- combn(single_probV_miivs, m = length(single_probV_miivs)-2, simpflify = T)
  #  # names(badmiivs)[p] <- paste0(setdiff(single_probV_miivs, paste0(combotable[,1])), collapse = "+")
  # }
  badmiivs <- combn(single_probV_miivs, m = 2, simpflify = T)
  return(badmiivs)
}

getbadmiivs_all <- function(probVs, fit){
  miivs_all_nosign <- getMIIVs_noplussign(fit)
  badmiivs <- list()
  single_probV_miivs <- list()
  for (p in 1:length(probVs)){
    single_probV_miivs[[p]] <- miivs_all_nosign[which(names(miivs_all_nosign)==probVs[[p]])][[1]]
    names(single_probV_miivs)[p] <- probVs[p]
  }
  return(single_probV_miivs)
}


#######
ProbMiivs_combo <- function(model, data){
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

# step1del_combo <- function(model, data, setup){
#   probVs <- setup$probVs
#   badmiivs <- setup$badmiivs
#   goodmiivs <- setup$goodmiivs
#   badmiivs_int <- setup$badmiivs_int
#   sig_sargantable_int <- setup$sig_sargantable_int
#
#   fit_del1 <- list()
#   newmiivs <- list()
#   delMiivsSargan <- list()
#   for (p in 1:length(probVs)){
#     ##p is the name of the problematic variable
#     fit_del1[[p]] <- list()
#     newmiivs[[p]] <- list()
#     delMiivsSargan[[p]] <- as.data.frame(matrix(NA, nrow = 2, ncol = ncol(badmiivs[[p]])) )
#     names(delMiivsSargan)[p] <- probVs[p]
#     rownames(delMiivsSargan[[p]]) <- c("sargandef", "newsargan")
#     for (i in 1:ncol(badmiivs[[p]])){
#       ##i for each deleted MIIV
#       ##then get the list of MIIVs
#       newmiivs[[p]][[i]] <- paste0(goodmiivs[[p]], sep = "\n",
#                                    paste0(names(badmiivs)[p], sep = "~",
#                                           paste0(paste0(badmiivs[[p]][,i]), collapse = "+")))
#       ##get the fit for each MIIV deletion
#       fit_del1[[p]][[i]] <- miive(model = model, data = data, var.cov = T, miiv.check = F,
#                                   instruments = newmiivs[[p]][[i]])
#       colnames(delMiivsSargan[[p]])[i] <- paste0(setdiff(badmiivs_int[[p]], paste0(badmiivs[[p]][,i])), collapse = "+")
#       #paste0(paste0(badmiivs[[p]][,i]), collapse = "+")
#       #paste0(setdiff(badmiivs_int[[p]], paste0(badmiivs[[p]][,i])), collapse = "+")
#
#       ##then get the sargan for each miiv deletion
#       delMiivsSargan[[p]][1,i] <- sig_sargantable_int[1,p] - fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
#       delMiivsSargan[[p]][2,i] <- fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
#     }
#   }
#   df1 <- qchisq(.95, df = 1)
#   #store the max chi square drop
#   maxdiff <- sapply(delMiivsSargan, function(i) max(i[1,]))
#   #max needs to be greater than df1
#   maxdiff <- maxdiff[maxdiff>df1]
#   #find list name
#   drop_location <- names(delMiivsSargan)[sapply(seq_along(delMiivsSargan),
#                                                 function(x) {any(delMiivsSargan[[x]]==maxdiff[x])})]
#   drop_location_num <- sapply(drop_location, function(x) match(x, names(delMiivsSargan)))
#
#   drop_MIIV <- sapply(drop_location_num, function(i)
#     colnames(delMiivsSargan[[i]])[delMiivsSargan[[i]][1,] == maxdiff[i]])
#   ##get the equation order for these variables (the fit$eqn[[p]])
#   sargan_eqnorder <- as.numeric(sig_sargantable_int[3,])
#
#   badmiivs_candidate <- list()
#   for (p in 1:length(drop_MIIV)){
#     badmiivs_candidate[[p]] <- strsplit(drop_MIIV[p],split ='+', fixed = T)[[1]]
#     names(badmiivs_candidate)[p] <- probVs[p]
#   }
#
#   finalobj <- list(maxdiff = maxdiff,
#                    drop_location = drop_location,
#                    drop_location_num = drop_location_num,
#                    drop_MIIV = drop_MIIV,
#                    delMiivsSargan = delMiivsSargan,
#                    badmiivs = badmiivs,
#                    badmiivs_all = badmiivs_candidate,
#                    sargan_eqnorder = sargan_eqnorder
#                    #badmiivs = badmiivs_2,
#                    #delMiivsSargan_nextstep = delMiivsSargan_nextstep
#   )
#   return(finalobj)
# }
#
# stepNdel_combo <- function(model, data, setup, stepPrevious_final){
#   probVs <- setup$probVs
#   badmiivs <- setup$badmiivs
#   goodmiivs <- setup$goodmiivs
#   badmiivs_int <- setup$badmiivs_int
#   sig_sargantable_int <- setup$sig_sargantable_int
#
#   badmiivs_new <- stepPrevious_final$badmiivs_all
#   ##first set up the new combos of miivs to be deleted/retained
#   for (p in 1:length(badmiivs_int)){
#     new_colnum <- length(badmiivs_new[[p]])
#     new_rownum <- length(badmiivs_int[[p]])-length(badmiivs_new[[p]])+1
#     badmiivs[[p]] <- matrix(NA, nrow = new_rownum, ncol = new_colnum)
#     for (i in 1:new_rownum){
#       badmiivs[[p]][i,] <- setdiff(badmiivs_int[[p]],badmiivs_new[[p]])[i]
#     }
#     for (q in 1:new_colnum){
#       badmiivs[[p]][new_rownum,q] <- badmiivs_new[[p]][q]
#     }
#   }
#
#   ##then repeat the steps in step1del_combo
#   fit_del1 <- list()
#   newmiivs <- list()
#   delMiivsSargan <- list()
#   for (p in 1:length(probVs)){
#     ##p is the name of the problematic variable
#     fit_del1[[p]] <- list()
#     newmiivs[[p]] <- list()
#     delMiivsSargan[[p]] <- as.data.frame(matrix(NA, nrow = 2, ncol = ncol(badmiivs[[p]])) )
#     names(delMiivsSargan)[p] <- probVs[p]
#     rownames(delMiivsSargan[[p]]) <- c("sargandef", "newsargan")
#     for (i in 1:ncol(badmiivs[[p]])){
#       ##i for each deleted MIIV
#       ##then get the list of MIIVs
#       newmiivs[[p]][[i]] <- paste0(goodmiivs[[p]], sep = "\n",
#                                    paste0(names(badmiivs)[p], sep = "~",
#                                           paste0(paste0(badmiivs[[p]][,i]), collapse = "+")))
#       ##get the fit for each MIIV deletion
#       fit_del1[[p]][[i]] <- miive(model = model, data = data, var.cov = T, miiv.check = F,
#                                   instruments = newmiivs[[p]][[i]])
#       colnames(delMiivsSargan[[p]])[i] <- paste0(setdiff(badmiivs_int[[p]], paste0(badmiivs[[p]][,i])), collapse = "+")
#       #paste0(paste0(badmiivs[[p]][,i]), collapse = "+")
#       #paste0(setdiff(badmiivs_int[[p]], paste0(badmiivs[[p]][,i])), collapse = "+")
#
#       ##then get the sargan for each miiv deletion
#       delMiivsSargan[[p]][1,i] <- sig_sargantable_int[1,p] - fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
#       delMiivsSargan[[p]][2,i] <- fit_del1[[p]][[i]]$eqn[[sig_sargantable_int[3,p]]]$sargan
#     }
#   }
#   df1 <- qchisq(.95, df = 1)
#   #store the max chi square drop
#   maxdiff <- sapply(delMiivsSargan, function(i) max(i[1,]))
#   #max needs to be greater than df1
#   maxdiff <- maxdiff[maxdiff>df1]
#   #find list name
#   drop_location <- names(delMiivsSargan)[sapply(seq_along(delMiivsSargan),
#                                                 function(x) {any(delMiivsSargan[[x]]==maxdiff[x])})]
#   drop_location_num <- sapply(drop_location, function(x) match(x, names(delMiivsSargan)))
#
#   drop_MIIV <- sapply(drop_location_num, function(i)
#     colnames(delMiivsSargan[[i]])[delMiivsSargan[[i]][1,] == maxdiff[i]])
#   ##get the equation order for these variables (the fit$eqn[[p]])
#   sargan_eqnorder <- as.numeric(sig_sargantable_int[3,])
#
#   badmiivs_candidate <- list()
#   for (p in 1:length(drop_MIIV)){
#     badmiivs_candidate[[p]] <- strsplit(drop_MIIV[p],split ='+', fixed = T)[[1]]
#     names(badmiivs_candidate)[p] <- probVs[p]
#   }
#
#   finalobj <- list(maxdiff = maxdiff,
#                    drop_location = drop_location,
#                    drop_location_num = drop_location_num,
#                    drop_MIIV = drop_MIIV,
#                    delMiivsSargan = delMiivsSargan,
#                    badmiivs = badmiivs,
#                    badmiivs_all = badmiivs_candidate,
#                    sargan_eqnorder = sargan_eqnorder
#                    #badmiivs = badmiivs_2,
#                    #delMiivsSargan_nextstep = delMiivsSargan_nextstep
#   )
#   return(finalobj)
#
# }
