library(lavaan)
##############simulate more data and try on algorithm 3############
#m1: one factor with x4 ~~ x5
sm1 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 + .65*x5 + .6*x6 + .6*x7 + .55*x8
          x4 ~~ .5*x5'
sim1 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim1[[p]] <- simulateData(sm1, sample.nobs = 1000)
}
#m2: 2 factor model
sm2 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
      f1 ~~ .5*f2'
sim2 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim2[[p]] <- simulateData(sm2, sample.nobs = 1000)
}
#m3: 2 factor model w/ x7 crossload on f1
sm3 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2'
sim3 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim3[[p]] <- simulateData(sm3, sample.nobs = 1000)
}
#m4: 1 factor w/ x4 ~~ x5, and x4 ~~ x2
sm4 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 + .65*x5 + .6*x6 + .6*x7 + .55*x8
          x4 ~~ .5*x5
          x4 ~~ .4*x2'
sim4 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim4[[p]] <- simulateData(sm4, sample.nobs = 1000)
}
#m5: 1 factor w/ x4 ~~ x5, x4 ~~ x2, and x1 ~~ x6
sm5 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 + .65*x5 + .6*x6 + .6*x7 + .55*x8
          x4 ~~ .5*x5
          x4 ~~ .4*x2
          x1 ~~ .4*x6'
sim5 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim5[[p]] <- simulateData(sm5, sample.nobs = 1000)
}
#m6: 2 factor model w/ x7 crossload on f1 and x2 crossload on f2
sm6 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8 + .5*x2
        f1 ~~ .5*f2'
sim6 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim6[[p]] <- simulateData(sm6, sample.nobs = 1000)
}
#m7: 2 factor model w/ x7 crossload on f1, x3 crossload on f2, and x1 crossload on f2
sm7 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8 + .5*x2 + .5*x1
        f1 ~~ .5*f2'
sim7 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim7[[p]] <- simulateData(sm7, sample.nobs = 1000)
}
#m8: 1 factor w/  x1 ~~ x6
sm8 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 + .65*x5 + .6*x6 + .6*x7 + .55*x8
          x1 ~~ .4*x6'
sim8 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim8[[p]] <- simulateData(sm8, sample.nobs = 1000)
}
#m9: 2 factor, with x4, x2, x5 on the second factor
sm9<- 'f1 =~ 1*x1 + .7*x3 + .6*x6 + .6*x7 + .55*x8
          f2 =~ 1*x4 + .8*x2 + .65*x5
          f1 ~~ .5*f2'
sim9 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim9[[p]] <- simulateData(sm9, sample.nobs = 1000)
}
#m10: 2 factor, with a crossloading and a pair of correlated error
sm10 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .6*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
      f1 ~~ .5*f2
        x2 ~~ .3*x3'
sim10<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim10[[p]] <- simulateData(sm10, sample.nobs = 1000)
}
#m11: 3 factor
sm11 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4
f2 =~ 1*x5 + .8*x6+ .7*x7 + .7*x8
f3 =~ 1*x9 + .8*x10 + .7*x11 + .x7*x12
f1 ~~ .4*f2
f1 ~~ .4*f3
f2 ~~ .4*f3'
sim11<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim11[[p]] <- simulateData(sm11, sample.nobs = 1000)
}
#sm12: 2 factor but only 2 variables on the 2nd factor
sm12 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .6*x5 + .6*x6
        f2=~ 1*x7 + .7*x8
      f1 ~~ .5*f2'
sim12<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim12[[p]] <- simulateData(sm12, sample.nobs = 1000)
}
#sm13: 2 factor but only 2 variables on the 2nd factor and correlated error (different factor)
sm13 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .6*x5 + .6*x6
        f2=~ 1*x7 + .7*x8
      f1 ~~ .5*f2
      x8 ~~ .4*x2'
sim13<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim13[[p]] <- simulateData(sm13, sample.nobs = 1000)
}

#sm14: 2 factor but only 2 variables on the 2nd factor and correlated error (same factor)
sm14 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .6*x5 + .6*x6
        f2=~ 1*x7 + .7*x8
      f1 ~~ .5*f2
      x3 ~~ .4*x2'
sim14<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim14[[p]] <- simulateData(sm14, sample.nobs = 1000)
}

#sm15: 2 factor but only 2 variables with crossloading x2 on f2
sm15 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .6*x5 + .6*x6
        f2=~ 1*x7 + .7*x8 + .5*x2
      f1 ~~ .5*f2'
sim15<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim15[[p]] <- simulateData(sm15, sample.nobs = 1000)
}
#sm15: 2 factor but only 2 variables with crossloading x8 on f1
sm16 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .6*x5 + .6*x6 + .5*x8
        f2=~ 1*x7 + .7*x8
      f1 ~~ .5*f2'
sim16<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim16[[p]] <- simulateData(sm16, sample.nobs = 1000)
}

