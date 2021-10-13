

library(lavaan)
library(MIIVsem)


#m3: 2 factor model w/ x7 crossload on f1
sm3 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2'
sim3 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim3[[p]] <- simulateData(sm3, sample.nobs = 1000)
}


#quick check on sim data
data <- sim3[[1]]
sigLevel <- .05
scalingCrit <- 'sargan+factorloading_R2'
s1 <- step1_E5(data, sigLevel, scalingCrit)
s2 <- step2_E5(s1, data, sigLevel, scalingCrit)


#messy empirical data

big5data <- read.csv(file="https://quantdev.ssri.psu.edu/sites/qdev/files/dataBIG5.csv",
                     header=TRUE)
dim(big5data)
big5data[big5data==0] <- NA #make 0s NAs
big5data <- na.omit(big5data)
dim(big5data)

data <- big5data[,8:57]

s1 <- step1_E5(data, sigLevel, scalingCrit)
s2 <- step2_E5(s1, data, sigLevel, scalingCrit)
