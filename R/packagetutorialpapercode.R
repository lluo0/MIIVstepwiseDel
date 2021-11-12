##efamiive paper code
library(psych)
library(lavaan)
library(MIIVsem)
###### simulation illustration########
##simple 2 factor model with one crossload

##2 factor model w/ x7 crossload on f1
sim_demo_model <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2'
set.seed(123.4)
sim_demo <- simulateData(sim_demo_model, sample.nobs = 1000)

# step 1: choose scaling indicator
data <- sim_demo
sigLevel <- .05
scalingCrit <- 'sargan+factorloading_R2'

###### sim 1 - 4 factor models with multiple crossloadings ####### 
sim1_model <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8 + .4*x11
        f3=~ 1*x9 + .75*x10 + .65*x11 + .55*x12 + .45*x15
        f4=~ 1*x13 + .85*x14 + .75*x15 + .6*x16 + .35*x3
        f1 ~~ .5*f2
        f1 ~~ .4*f3
        f1 ~~ .45*f4
        f2 ~~ .35*f3
        f2 ~~ .3*f4
        f3 ~~ .5*f3'
set.seed(123.4)
sim1 <- simulateData(sim1_model, sample.nobs = 1000)

###### sim 2 - including user-specified correlated errors #######
sim2_model <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8 
        f1 ~~ .5*f2
        x1 ~~ .5*x5
      x2 ~~ .55*x6
x3~~ .6*x7
x4~~.65*x8'
set.seed(123.4)
sim2 <- simulateData(sim2_model, sample.nobs = 1000)
EFAmiive5(sim2, .05,'sargan+factorloading_R2' )
EFAmiive5(sim2, .05,'sargan+factorloading_R2' , 
          correlatedErrors = 'x1~~x5
          x2~~x6
          x3~~x7
          x4~~x8')

miive('f1=~x1+x2+x3+x4
      f2=~x5+x6+x7+x8
      x1~~x5
      x2~~x6
      x3~~x7
      x4~~x8', sim2, var.cov = T)
###### empirical 1 - HolzingerSwineford1939’s data######
library(lavaan)
data('HolzingerSwineford1939')

miive('f1=~x1+x2+x3
      f2=~x4+x5+x6
      f3=~x7+x8+x9',HolzingerSwineford1939[,-c(1:6)], var.cov = T )
EFAmiive5(HolzingerSwineford1939[,-c(1:6)],.05, 'sargan+factorloading_R2')

miive('f1=~x6+x4+x5+x8\nf2=~x2+x1+x3+x9+x8\nf3=~x7+x9+x8', HolzingerSwineford1939[,-c(1:6)], var.cov = T )
###### empirical 2 - Bergh's data ######
library(MPsychoR)
data("Bergh")

miive(model = 'f1=~ DP+EP+SP+HP
      f2=~ O1+O2+O3
      f3=~ A1+A2+A3', data = Bergh, var.cov = T)



EFAmiive5(Bergh[,1:10], .05, 'sargan+factorloading_R2')


###### empirical 3 - Ken's perceived accessibility data ######
accdata <- read.table("/Users/lanluo/Downloads/access_raw.txt",header=F,sep=",")
colnames(accdata) <- c(sapply(c(1:6), function(x) paste0('access', x)), sapply(c(1:6), function(x) paste0('easy', x)))


miive(paste0(paste0('f1=~', paste0(sapply(c(1:6), function(x) paste0('access', x)), collapse = '+')), '\n',
             paste0('f2=~', paste0(sapply(c(1:6), function(x) paste0('easy', x)), collapse = '+'))),
      accdata, var.cov = T)


EFAmiive5(accdata, .05, 'sargan+factorloading_R2')
EFAmiive5(accdata, .05, 'sargan+factorloading_R2',
          correlatedErrors = 'access1~~easy1
      access2~~easy2
      access3~~easy3
      access4~~easy4
      access5~~easy5
      access6~~easy6')

#3
sigLevel <- .05
scalingCrit <- 'sargan+factorloading_R2'
data <- sim2
correlatedErrors <- 'x1~~x5
          x2~~x6
          x3~~x7
          x4~~x8'

s1 <- step1_E5(data, sigLevel, scalingCrit, correlatedErrors)
s1
s2 <- step2_E5(s1, data, sigLevel, scalingCrit)
s3 <- stepN_E5(s2, data, sigLevel, scalingCrit)
