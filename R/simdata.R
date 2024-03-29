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
#sm16: 2 factor but only 2 variables with crossloading x8 on f1
sm16 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .6*x5 + .6*x6 + .5*x8
        f2=~ 1*x7 + .7*x8
      f1 ~~ .5*f2'
sim16<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim16[[p]] <- simulateData(sm16, sample.nobs = 1000)
}

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

#sm19:
#m19: 2 factor model, with x3 having correlated errors with x4 and x5, but not between x4 and x5.
sm19 <- 'f1 =~ 1*x1 + .8 * x2 + .75*x3 + .7*x4 + .6*x5
        f2=~ 1*x6 + .8*x7 + .75*x8 + .7*x9 + .6*x10
      f1 ~~ .5*f2
      x3 ~~ .4*x4
      x3 ~~ .35*x5
      '
sim19 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim19[[p]] <- simulateData(sm19, sample.nobs = 1000)
}

#sm20
#20: 2 factor model, with correlated errors
sm20 <- 'f1 =~ 1*x1 + .8 * x2 + .75*x3 + .7*x4 + .6*x5
        f2=~ 1*x6 + .8*x7 + .75*x8 + .7*x9 + .6*x10
      f1 ~~ .5*f2
      x2 ~~ .4*x4
      x2 ~~ .35*x5
      x4 ~~ .3*x5
      '
sim20 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim20[[p]] <- simulateData(sm20, sample.nobs = 1000)
}

#sm21
#21: 3 factor model, with correlated errors
sm21 <- 'f1 =~ 1*x1 + .8 * x2 + .75*x3 + .7*x4 + .6*x5
        f2=~ 1*x6 + .8*x7 + .75*x8 + .7*x9 + .6*x10
        f3=~ 1*x11 + .8*x12 + .7*x13 + .65*x14 + .6*x15
      f1 ~~ .5*f2
      x2 ~~ .4*x4
      x2 ~~ .35*x5
      x4 ~~ .3*x5
      f1 ~~ .4*f3
      f2 ~~ .35*f3
      '
sim21 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim21[[p]] <- simulateData(sm21, sample.nobs = 1000)
}

#sm22
#22: 4 factor model, with correlated errors
sm22 <- 'f1 =~ 1*x1 + .8 * x2 + .75*x3 + .7*x4 + .6*x5
        f2=~ 1*x6 + .8*x7 + .75*x8 + .7*x9 + .6*x10
        f3=~ 1*x11 + .8*x12 + .7*x13 + .65*x14 + .6*x15
        f4=~ 1*x16+ .8*x17 + .7*x18+ .65*x19 + .6*x20
      f1 ~~ .5*f2
      x2 ~~ .4*x4
      x2 ~~ .35*x5
      x4 ~~ .3*x5
      f1 ~~ .4*f3
      f2 ~~ .35*f3
      f4~~ .5*f1
      f4~~.4*f2
      f4~~.3*f3
      '
sim22 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim22[[p]] <- simulateData(sm22, sample.nobs = 1000)
}
########small sample size sim#######
#m1: one factor with x4 ~~ x5
sm1 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 + .65*x5 + .6*x6 + .6*x7 + .55*x8
          x4 ~~ .5*x5'
sim1_50 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim1_50[[p]] <- simulateData(sm1, sample.nobs = 50)
}
sim1_100 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim1_100[[p]] <- simulateData(sm1, sample.nobs = 100)
}
#m2: 2 factor model
sm2 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
      f1 ~~ .5*f2'
sim2_50 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim2_50[[p]] <- simulateData(sm2, sample.nobs = 50)
}
sim2_100 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim2_100[[p]] <- simulateData(sm2, sample.nobs = 100)
}
#m3: 2 factor model w/ x7 crossload on f1
sm3 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2'
sim3_50 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim3_50[[p]] <- simulateData(sm3, sample.nobs = 50)
}
sim3_100 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim3_100[[p]] <- simulateData(sm3, sample.nobs = 100)
}
#m4: 1 factor w/ x4 ~~ x5, and x4 ~~ x2
sm4 <- 'f1 =~ 1*x1 + .8*x2 + .7*x3 + .7*x4 + .65*x5 + .6*x6 + .6*x7 + .55*x8
          x4 ~~ .5*x5
          x4 ~~ .4*x2'
sim4_50 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim4_50[[p]] <- simulateData(sm4, sample.nobs = 50)
}
sim4_100 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim4_100[[p]] <- simulateData(sm4, sample.nobs = 100)
}
