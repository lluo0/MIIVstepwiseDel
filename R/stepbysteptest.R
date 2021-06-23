miive('f1=~x5+x6+x8+x4+x1+x7
      f2=~x2+x3+x4+x1+x7', sim3[[1]], var.cov = T)

#step by step test on holzinger data
data <- HolzingerSwineford1939[,7:15]
threshold <- .05
priority <- 'sargan+factorloading_R2'

s1 <- step1_E5(data, threshold, priority)
s2 <- step2_E5(s1, data, threshold, priority)

select_scalingind_stepN(data, threshold, priority, s1)

select_scalingind(data[,c(1:3,7:9)], threshold, priority)

select_scalinging(data[,c(1:3,7:9)], threshold, priority)

select_scalingind_stepN(data, threshold, priority, s1)

data <- sim3[[1]]
s1 <- step1_E5(data, threshold, priority)
s2 <- step2_E5(s1, data, threshold, priority)

stepN_E5(s2, data, threshold, priority)

select_scalingind(data, threshold, 'sargan_factorloading')
r2_order(data)


##6.15
#sim6
data <- sim6[[1]] #2 crossloading variables
threshold <- .05
priority <- 'sargan+factorloading_R2'

s1 <- step1_E5(data, threshold, priority)
s2 <- step2_E5(s1, data, threshold, priority)
s3 <- stepN_E5(s2, data, threshold, priority)


step2_E5 <- function(stepPrev, data, threshold, priority){

  badvar <- stepPrev$badvar
  goodmodelpart <- stepPrev$goodmodelpart
  num_factor <- stepPrev$num_factor

  # scalingindicator <- select_scalingind_stepN(data, threshold, priority, goodmodelpart, badvar, num_factor)
  scalingindicator <- select_scalingind_stepN(data, threshold, priority, stepPrev)

  order_scalingind <- which(badvar==scalingindicator)

  model <- paste(paste0(goodmodelpart, collapse = '\n'),
                 paste(paste0("f",num_factor+1), "=~",paste0(badvar[order_scalingind]), '+',
                       paste0(badvar[-order_scalingind], collapse = "+"), sep = ""),
                 sep = "\n")
  fit <- miive(model, data, var.cov = T)
  badvar <- getbadvar(fit, threshold)
  num_badvar <- length(badvar)

  if(num_badvar==0){
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = num_factor,
                     num_badvar = 0)
  }
  if(num_badvar!=0){
    # #save the problematic variables
    # badvar <- getbadvar(fit, threshold)
    #save the good variables - update the ones from the previous output
    goodvar <- stepPrev$goodvar
    #goodvar[[num_factor]] <- setdiff(stepPrev$badvar, badvar)
    goodvar[[num_factor+1]] <- setdiff(stepPrev$badvar, badvar)
    #reorder the goodvar part again so the new scaling indicator appears first on the new list
    goodvar[[num_factor+1]] <- goodvar[[num_factor+1]][c(match(scalingindicator, goodvar[[num_factor+1]]),
                                                     setdiff(order(goodvar[[num_factor+1]]),match(scalingindicator, goodvar[[num_factor+1]])))]

    # goodmodelpart <- paste(stepPrev$goodmodelpart,
    #                        paste(paste0("f", num_factor), "=~",
    #                              paste(goodvar[[num_factor]], collapse = "+"), sep = ""),
    #                        sep = "\n")
    goodmodelpart <- list(stepPrev$goodmodelpart,
                          paste(paste0("f", num_factor+1), "=~",
                                paste(goodvar[[num_factor+1]], collapse = "+"), sep = ""))
    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = stepPrev$num_factor+1,
                     num_badvar = num_badvar,
                     goodvar = goodvar,
                     badvar = badvar,
                     goodmodelpart = goodmodelpart,
                     nextstep = ifelse(length(badvar!=0), 'yes', 'no'))
  }
  return(finalobj)
}

stepN_E5 <- function(stepPrev, data, threshold, priority){
  ##first crossload the bad variables
  crossloadmodel <- lapply(stepPrev$goodmodelpart, function(x)
    paste0(x, '+',paste0(stepPrev$badvar, collapse = '+')))
  crossloadfit <- miive(paste0(crossloadmodel, collapse = '\n'), data, var.cov = T)
  ##then see if any variables actually crossload
  ##update the bad variable list
  newbadvar <- getbadvar_crossload(crossloadfit, threshold, stepPrev$num_factor, stepPrev$badvar)
  ##add the crossloaded variables to the model
  # newgoodvar_addon <- lapply(newbadvar, function(x) setdiff(x, stepPrev$badvar))
  newgoodvar_addon <- list()
  for(p in 1:length(newbadvar)){
    if(length(newbadvar[[p]])==0){
      newgoodvar_addon[[p]] <- stepPrev$badvar
    }
    if(identical(newbadvar[[p]], stepPrev$badvar)){
      newgoodvar_addon[[p]] <- NULL
    }
    if(!length(newbadvar[[p]])==0 & !identical(newbadvar[[p]], stepPrev$badvar) ){
      newgoodvar_addon[[p]] <- setdiff(stepPrev$badvar,newbadvar[[p]])
    }
  }


#move this down into the if function
  # #create new goodvar and goodmodelvar part based on new crossload addition. 060821
  # newgoodmodelpart <-  mapply(function(x,y) paste0(x, '+', y), stepPrev$goodmodelpart,
  #                             lapply(newgoodvar_addon, function(i) paste0(i, collapse = '+')), SIMPLIFY=FALSE)
  # newgoodvar <- mapply(c, stepPrev$goodvar, newgoodvar_addon, SIMPLIFY  = F)
  #create new badvar
  newbadvar <- setdiff(stepPrev$badvar, unique(unlist(newgoodvar_addon)))

  if(length(newbadvar)<=1){ #we do not need to create a new factor
    #aka after crossloading we either no longer have problematic variables or have only one.

    #add badvar back to the whole model, if any
    if(length(newbadvar==1)){
      #create new goodvar and goodmodelvar part based on new crossload addition. 060821
      mapply(function(x,y) paste0(x, '+', y), stepPrev$goodmodelpart,
             lapply(newgoodvar_addon, function(i) paste0(i, collapse = '+')), SIMPLIFY=FALSE)
      newgoodvar <- mapply(c, stepPrev$goodvar, newgoodvar_addon, SIMPLIFY  = F)

      model <- paste0(paste0(newgoodmodelpart[[-length(newgoodmodelpart)]], collapse = '\n'),
                      paste0(newgoodmodelpart[[length(newgoodmodelpart)]], "+", newbadvar))
    }
    if(length(newbadvar==0)){
      model <- paste0(newgoodmodelpart, collapse = '\n')
    }
    fit <- miive(model, data, var.cov = T)
    badvar <- getbadvar(fit, threshold)
    num_badvar <- length(badvar)

    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = stepPrev$num_factor, #NOTE: same number of factors as the previous step
                     num_badvar = num_badvar,
                     goodvar = newgoodvar,
                     badvar = badvar,
                     goodmodelpart = newgoodmodelpart,
                     nextstep = 'no')
  }
  if(length(newbadvar)>1){ #aka we still need to create a new factor

    #create new goodvar and goodmodelvar part based on new crossload addition. 060821
    newgoodmodelpart <-  mapply(function(x,y) paste0(x, '+', y), stepPrev$goodmodelpart,
                                lapply(newgoodvar_addon, function(i) paste0(i, collapse = '+')), SIMPLIFY=FALSE)
    newgoodvar <- mapply(c, stepPrev$goodvar, newgoodvar_addon, SIMPLIFY  = F)

    scalingindicator <- select_scalingind_stepN(data, threshold, priority, stepPrev)
    order_scalingind <- which(newbadvar==scalingindicator)

    model <- paste(paste0(newgoodmodelpart, collapse = '\n'),
                   paste(paste0("f",stepPrev$num_factor+1), "=~",paste0(newbadvar[order_scalingind]), '+',
                         paste0(newbadvar[-order_scalingind], collapse = "+"), sep = ""),
                   sep = "\n")
    fit <- miive(model, data, var.cov = T)
    badvar <- getbadvar(fit, threshold)
    num_badvar <- length(badvar)

    finalobj <- list(model = model,
                     fit  = fit,
                     num_factor = stepPrev$num_factor+1,
                     num_badvar = num_badvar,
                     goodvar = newgoodvar,
                     badvar = badvar,
                     goodmodelpart = newgoodmodelpart,
                     nextstep = ifelse(length(badvar!=0), 'yes', 'no'))
  }
  return(finalobj)
}

#holzinger data
data <- lavaan::HolzingerSwineford1939[,7:15]
threshold <- .05
priority <- 'sargan+factorloading_R2'

s1 <- step1_E5(data, threshold, priority)
s2 <- step2_E5(s1, data, threshold, priority)
s3 <- stepN_E5(s2, data, threshold, priority)
s4 <- stepN_E5(s3, data, threshold, priority)
stepN_E5(s4, data, threshold, priority)


stepPrev <- s2

miive(model[[1]], data, var.cov=T)

miive('f1=~x6+x4+x5+x8+x9+x7+x3
      f2=~x1+x2+x8+x9+x7+x3', data, var.cov = T)

miive('f1=~x6+x4+x5+x3
      f2=~x1+x2+x3
      f3=~x7+x9', data, var.cov = T)

miive('f1=~x1+x2+x3+x4+x5+x6+x7+x8+x9', data, var.cov = T)

miive('f1=~x6+x2+x3+x4+x5+x1+x7+x8+x9', data, var.cov = T)

summary(cfa('f1=~x1+x2+x3+x4+x5+x6+x7+x8+x9',data))

miive('f1=~x6+x4+x5+x3+x7+x8+x9
      f2=~x1+x2+x3+x7+x8+x9', data, var.cov = T)

miive('f1=~x1+x2+x3
      f2=~x4+x5+x6
      f3=~x7+x8+x9', data, var.cov = T)

summary(cfa('f1=~x1+x2+x3
      f2=~x4+x5+x6
      f3=~x7+x8+x9', data))
lavaanfit <- cfa('f1=~x1+x2+x3
      f2=~x4+x5+x6
      f3=~x7+x8+x9', data)
lavaanfit

atent Variables:
  Estimate  Std.Err  z-value  P(>|z|)
f1 =~
  x1                1.000
x2                0.554    0.100    5.554    0.000
x3                0.729    0.109    6.685    0.000
f2 =~
  x4                1.000
x5                1.113    0.065   17.014    0.000
x6                0.926    0.055   16.703    0.000
f3 =~
  x7                1.000
x8                1.180    0.165    7.152    0.000
x9                1.082    0.151    7.155    0.000


miive('f1=~x1+x2+x3+x8+x9
      f2=~x4+x5+x6+x3+x8+x9
      f3=~x7+x8+x9', data, var.cov = T)

miivs('f1=~x1+x2+x3+x8+x9
      f2=~x6+x5+x4+x3
      f3=~x7+x8+x9')

miive('f1=~x1+x2+x3
      f2=~x6+x5+x4+x3
      f3=~x7+x8+x9+x4',data, var.cov = T)

miive()

miive('')

f1=~x6+x4+x5+x3\nf2=~x1+x2+x3\nf3=~x8+x7+x9

miive('f1=~x6+x4+x5+x3
f2=~x1+x2+x3+x9
f3=~x8+x7+x9', data, var.cov = T)


##6.21.21
#sim3
data <- sim3[[1]]
threshold <- .05
priority <- 'sargan+factorloading_R2'
s1 <- step1_E5(data, threshold, priority)
s2 <- step2_E5(s1,data, threshold, priority)
stepN_E5(s2, data, threshold, priority)

#sm17: 3 factor with crossloading, x7 and x4 crossload
sm17 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 +.5*x7
f2 =~ 1*x5 + .8*x6+ .7*x7 + .7*x8
f3 =~ 1*x9 + .8*x10 + .7*x11 + .x7*x12+.6*x4
f1 ~~ .4*f2
f1 ~~ .4*f3
f2 ~~ .4*f3'
sim17<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim17[[p]] <- simulateData(sm17, sample.nobs = 1000)
}
#test
data <- sim18[[1]]
s1 <- step1_E5(data, threshold, priority)
s2 <- step2_E5(s1, data, threshold, priority)
s3 <- stepN_E5(s2, data, threshold, priority)
s3_new <-  stepN_E5(s2, data, threshold, priority)
s4 <- stepN_E5(s3, data, threshold, priority)
stepPrev <- s3
s5 <- stepN_E5(s4, data, threshold, priority)

sim17out <- EFAmiive5(sim2[[1]], threshold, priority)

data <- sim18[[1]]

EFAmiive5(sim17[[1]], threshold, priority)

EFAmiive5(sim4[[1]], threshold, priority)

s1 <- step1_E5(sim13[[1]], threshold, priority)
step2_E5(s1, sim13[[1]], threshold, priority)
step1_E5(sim18[[1]], threshold, priority)

data <- sim18[[1]]



testlist <- list()
testlist[[1]] <- NULL
testlist[[2]] <- NULL

miive('f1=~x1+x3+x6+x7+x8\nf2=~x5+x4+x2', data, var.cov = T)
EFAmiive5( HolzingerSwineford1939[,7:15], threshold, priority)
data <-  HolzingerSwineford1939[,7:15]


##6/23
#m2: 2 factor model
sm2 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
      f1 ~~ .5*f2'
sim2 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim2[[p]] <- simulateData(sm2, sample.nobs = 1000)
}

#model fit
sim2out <-  EFAmiive5(sim2[[1]], threshold = .05,  priority = 'sargan+factorloading_R2')
sim2out
#m3: 2 factor model w/ x7 crossload on f1
sm3 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2'
sim3 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim3[[p]] <- simulateData(sm3, sample.nobs = 1000)
}

#model fit
sim3out <-  EFAmiive5(sim3[[1]], threshold = .05,  priority = 'sargan+factorloading_R2')
sim3out
#m4: 1 factor w/ x4 ~~ x5, and x4 ~~ x2
sm4 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 + .65*x5 + .6*x6 + .6*x7 + .55*x8
          x4 ~~ .5*x5
          x4 ~~ .4*x2'
sim4 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim4[[p]] <- simulateData(sm4, sample.nobs = 1000)
}

#model fit
sim4out <-  EFAmiive5(sim4[[1]], threshold = .05,  priority = 'sargan+factorloading_R2')
sim4out
#sm17: 3 factor with crossloading, x7 and x4 crossload
sm17 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 +.5*x7
f2 =~ 1*x5 + .8*x6+ .7*x7 + .7*x8
f3 =~ 1*x9 + .8*x10 + .7*x11 + .x7*x12+.6*x4
f1 ~~ .4*f2
f1 ~~ .4*f3
f2 ~~ .4*f3'
sim17<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim17[[p]] <- simulateData(sm17, sample.nobs = 1000)
}

sim17out <- EFAmiive5(sim17[[1]], threshold = .05,  priority = 'sargan+factorloading_R2')
sim17out
#sm18: 4 factor
sm18 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 +.5*x7
f2 =~ 1*x5 + .8*x6+ .7*x7 + .7*x8
f3 =~ 1*x9 + .8*x10 + .7*x11 + .x7*x12+.6*x4
f4 =~ 1*x13 + .8*x14 + .7*x15 + .7*x16 + .5*x10
f1 ~~ .4*f2
f1 ~~ .4*f3
f2 ~~ .4*f3
f1 ~~ .5*f4
f2 ~~ .45*f4
f3 ~~ .4*f4'
sim18<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim18[[p]] <- simulateData(sm18, sample.nobs = 1000)
}

sim18out <- EFAmiive5(sim18[[1]], threshold = .05,  priority = 'sargan+factorloading_R2')
sim18out
