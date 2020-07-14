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

