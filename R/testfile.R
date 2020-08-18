#
library(MIIVsem)
library(lavaan)

model1 <- 'f1=~x1+x2+x3+x4+x5+x6+x7+x8'
model2 <- 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8'

#
sim2_setup <- ProbMiivs_c(model1, sim2[[1]])
sim2_step1 <- step1del_c(model1, sim2[[1]], sim2_setup)
sim2_step2 <- stepNdel_c(model1, sim2[[1]], sim2_setup, sim2_step1)

obj_sim6_step3 <- stepNdel_combo(model1, sim6[[1]], obj_sim6, obj_sim6_step2)

sim2_step1B <- step1del_combo(model1, sim2[[2]], sim2_setup)
sim2_step2B <- stepNdel_combo(model1, sim2[[2]], sim2_setup, sim2_step1B)
sim2_step3B <- stepNdel_combo(model1, sim2[[2]], sim2_setup, sim2_step2B)
sim2_step4B <- stepNdel_combo(model1, sim2[[2]], sim2_setup, sim2_step3B)

MIIVdel(model1, sim3[[1]])
MIIVdel(model2, sim3[[1]])

MIIVdel(model2, sim10[[1]])
MIIVdel(model2, sim10[[2]])
MIIVdel(model2, sim10[[4]])

model3f <- 'f1=~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12'
MIIVdel(model3f, sim11[[1]])


miive(model1, sim4[[3]], var.cov = T, miiv.check = F,
      instruments = '
      x2 ~ x3+x5+x6+x8
      x3 ~ x2+x4+x5+x6+x8
      x4 ~ x3+x6+x8+x7
      x5 ~ x2+x3+x6+x8+x7
      x6 ~ x2+x3+x4+x5+x8
      x7 ~ x2+x3+x4+x5+x6+x8
      x8 ~ x2+x3+x4+x5+x6+x7')

# x4 < x2+x5+x7 < step1  > sargan for x4: 1.054
#
# > step2
# x2+x5 > sargan for x4:  5.335
# x2+x7
# x5+x7


miive(model1, sim1[[3]], var.cov = T, miiv.check = F,
      instruments = '
      x2 ~ x3+x7+x4+x5+x6+x8
      x3 ~ x2+x4+x5+x6+x8+x7
      x4 ~ x2+x6+x8+x7+x3
      x5 ~ x2+x3+x6+x8+x7
      x6 ~ x2+x3+x4+x5+x8+x7
      x7 ~ x2+x3+x4+x5+x6+x8
      x8 ~ x2+x3+x4+x5+x6+x7')

miive(model1, sim1[[1]], var.cov = T, miiv.check = F,
      instruments = '
      x2 ~ x3+x7+x4+x5+x6+x8
      x3 ~ x2+x4+x5+x6+x8+x7
      x4 ~ x6+x8+x7+x3
      x5 ~ x2+x3+x6+x8+x7
      x6 ~ x2+x3+x4+x5+x8+x7
      x7 ~ x2+x3+x4+x5+x6+x8
      x8 ~ x2+x3+x4+x5+x6+x7')

miive(model2, sim10[[1]], var.cov = T, miiv.check = F,
      instruments = '
      x2 ~ x3+x4+x5+x6+x7+x8
      x3 ~ x2+x4+x5+x6+x7+x8
      x4 ~ x2+x3+x5+x6+x7+x8
      x6 ~ x1+x2+x3+x4+x7+x8
      x7 ~ x1+x2+x3+x4+x6+x8
      x8 ~ x1+x2+x3+x4+x6+x7')


MIIVdel_new(model2, sim10[[1]])

##test sim12-15
model3 <- 'f1=~x1+x2+x3+x4+x5+x6
f2=~x7+x8'
#sim12
MIIVdel_new(model1, sim12[[3]])
miive('f1=~x1+x2+x3+x4+x5+x6+x7+x8 \n x7~~x8', sim12[[1]], var.cov = T)
#sim13
MIIVdel_new(model1, sim13[[1]])
miive('f1=~x1+x2+x3+x4+x5+x6+x7+x8 \n x2~~x8 \n x7~~x8', sim13[[1]], var.cov = T)

#sim14
MIIVdel_new(model1, sim14[[1]])
#sim15
MIIVdel_new(model1, sim15[[1]])
MIIVdel_new(model1, sim15[[2]])
miive('f1=~x1+x3+x4+x5+x6
f2=~x7+x8+x2', sim15[[1]], var.cov = T)
MIIVdel_new('f1=~x1+x3+x4+x5+x6
f2=~x7+x8+x2', sim15[[1]])
#sim16
MIIVdel_new(model1, sim16[[1]])
miive(model3, sim16[[3]], var.cov = T)
miive(model3, sim16[[3]], var.cov = T)


####
#sim1: one factor with x4 ~~ x5
MIIVdel_new(model1, sim1[[1]])
#sim2: 2 factor model
MIIVdel_new(model1, sim2[[1]])
#sim3: 2 factor model w/ x7 crossload on f1
MIIVdel_new(model1, sim3[[1]])
MIIVdel_new(model2, sim3[[1]])
#sim4: 1 factor w/ x4 ~~ x5, and x4 ~~ x2
MIIVdel_new(model1, sim4[[1]])
#sim5: 1 factor w/ x4 ~~ x5, x4 ~~ x2, and x1 ~~ x6
MIIVdel_new(model1, sim5[[2]])
#sim6: 2 factor model w/ x7 crossload on f1 and x2 crossload on f2
MIIVdel_new(model1, sim6[[4]])
MIIVdel_new('f1=~x1+x3+x4 \n f2=~x5+x6+x7+x8+x2', sim6[[5]])
#sim7: 2 factor model w/ x7 crossload on f1, x3 crossload on f2, and x1 crossload on f2
MIIVdel_new(model1, sim7[[1]])
MIIVdel_new(model1, sim7[[2]])
#sim8: 1 factor w/  x1 ~~ x6
MIIVdel_new(model1, sim8[[1]])
#sim9: 2 factor, with x4, x2, x5 on the second factor
MIIVdel_new(model1, sim9[[1]])
miive('f1=~x1+x2+x3+x4+x5+x6+x7+x8 \n x2 ~~ x4 \n x2~~x5 \n x4~~x5', sim9[[1]], var.cov = T)
#sim10: 2 factor, with a crossloading and a pair of correlated error
MIIVdel_new(model1, sim10[[1]])



########negative construct in fisher's data####
load("/Users/LanLuo/Google Drive/19 spring/859 time series/data/Fisher/FisherData.Rdata")

negfisher <- list()
for (p in 1:length(FisherDataInterp)){
   negfisher[[p]] <-FisherDataInterp[[p]][,c(4:7, 9,11,12,16)]
}

#ind1
EFAmiive3(negfisher[[1]])
MIIVdel_new('f1=~guilty+concentrate+worried+hopeless+anhedonia+irritable+restless+down', negfisher[[1]])
miive('f1=~guilty+concentrate+worried+hopeless+anhedonia+irritable+restless+down \n down~~anhedonia', negfisher[[1]], var.cov = T)
#ind2
EFAmiive3(negfisher[[2]])
MIIVdel_new('f1 =~ down+hopeless+guilty+anhedonia+restless + irritable+worried + concentrate', negfisher[[2]])
miive('f1 =~ down+hopeless+guilty+anhedonia   \n f2=~restless + irritable+worried + concentrate + anhedonia
      \n concentrate~~hopeless ', negfisher[[2]], var.cov = T)
#ind3
EFAmiive3(negfisher[[3]])
MIIVdel_new('f1 =~ restless + worried + hopeless+concentrate +guilty +down+irritable+anhedonia',negfisher[[3]])
MIIVdel_new('f1 =~ restless + worried + hopeless \n f2=~concentrate +guilty +down+irritable+anhedonia',negfisher[[3]])
miive('f1 =~ restless + worried + hopeless \n f2=~concentrate +guilty +down+irritable+anhedonia \n anhedonia~~guilty', negfisher[[3]], var.cov = T)
miive('f1=~restless+worried+hopeless \n f2=~concentrate+guilty+down+irritable+anhedonia \n irritable~~anhedonia', negfisher[[3]], var.cov = T)
#ind4
EFAmiive3(negfisher[[4]])
MIIVdel_new('f1=~anhedonia+guilty+down+hopeless+concentrate+worried+restless+irritable', negfisher[[4]])
miive('f1=~anhedonia+guilty+down+hopeless \n f2=~concentrate+worried+restless+irritable+down', negfisher[[4]], var.cov = T)


##efamiive3 on simulations
#sim 7
EFAmiive3(sim7[[1]])
miive('f1=~x5+x6+x8 \n f2=~x1+x2+x7 \n f3=~x3+x4+x7', sim7[[1]], var.cov = T)
#sim 8
EFAmiive3(sim8[[1]])
#sim 9
EFAmiive3(sim9[[1]])
EFAmiive3(sim9[[2]])
EFAmiive3(sim9[[3]])
#sim 10
EFAmiive3(sim10[[1]])
EFAmiive3(sim10[[2]])
MIIVdel_new('f1=~x5+x6+x8 \n f2=~x1+x2+x3+x4+x7', sim10[[1]])
MIIVdel_new('f1=~x5+x6+x8 \n f2=~x2+x3+x4+x7+x1', sim10[[1]])


##MIIVdel on more negfisher
#ind5 the weird data
MIIVdel_new('f1=~ concentrate+irritable +restless+worried +guilty+hopeless+anhedonia+down', negfisher[[5]])
#ind6

#ind7
MIIVdel_new('f1=~ irritable+hopeless+guilty+restless+down+anhedonia+concentrate+worried', negfisher[[7]])
miive('f1=~ irritable+hopeless+guilty+restless+down+anhedonia+concentrate+worried
      \n concentrate ~~ hopeless \n concentrate ~~restless \n worried~~anhedonia',
      negfisher[[7]], var.cov = T)
miive('f1=~ irritable+hopeless+guilty+restless+down+anhedonia+concentrate+worried
      ',
      negfisher[[7]], var.cov = T)
#ind8
MIIVdel_new('f1=~worried+concentrate+anhedonia+hopeless+restless +guilty+down+irritable', negfisher[[8]])
miive('f1=~worried+concentrate+anhedonia+hopeless+restless +guilty+down+irritable
      \n restless~~irritable \n guilty ~~ down ', negfisher[[8]], var.cov = T)
miive('f1=~worried+concentrate+anhedonia+hopeless+restless +guilty+down+irritable
      \n restless~~irritable \n guilty ~~ down \n restless~~down
      \n guilty~~hopeless', negfisher[[8]], var.cov = T)
MIIVdel_new('f1=~worried+concentrate+anhedonia+hopeless \n f2=~restless +guilty+down+irritable', negfisher[[8]])
miive('f1=~worried+concentrate+anhedonia+hopeless \n f2=~restless +guilty+down+irritable \n irritable~~down',
      negfisher[[8]], var.cov = T)
#ind9
MIIVdel_new('f1=~ down + hopeless +irritable+guilty+restless+anhedonia+concentrate+worried', negfisher[[9]])
miive('f1=~ down + hopeless +irritable+guilty+restless+anhedonia+concentrate+worried \n concentrate ~~worried',
      negfisher[[9]], var.cov = T)
miive('f1=~ down + hopeless +irritable+guilty+restless+anhedonia+concentrate+worried
      \n concentrate ~~anhedonia \n worried~~guilty',
      negfisher[[9]], var.cov = T)
#ind10
MIIVdel_new('f1=~hopeless+guilty+worried+down+anhedonia+irritable+restless+concentrate', negfisher[[10]])
miive('f1=~hopeless+guilty+worried+down+anhedonia+irritable+restless+concentrate
      \n down~~concentrate \n anhedonia ~~down \n irritable~~restless \n concentrate ~~worried',
      negfisher[[10]],var.cov = T)
MIIVdel_new('f1=~hopeless+guilty+worried \n f2=~down+anhedonia+irritable+restless+concentrate', negfisher[[10]])
miive('f1=~hopeless+guilty+worried \n f2=~down+anhedonia+irritable+restless+concentrate
      \n irritable~~restless \n concentrate~~anhedonia', negfisher[[10]], var.cov = T)
#ind11
MIIVdel_new('f1=~ hopeless + guilty+down+restless+worried +concentrate +irritable+anhedonia', negfisher[[11]])
miive('f1=~ hopeless + guilty+down+restless +irritable
      \n f2=~worried + concentrate + anhedonia
      \n anhedonia ~~down', negfisher[[11]], var.cov = T)
MIIVdel_new('f1=~ hopeless + guilty+down+restless \n f2=~worried +concentrate+irritable+anhedonia', negfisher[[11]])
miive('f1=~ hopeless + guilty+down+restless \n f2=~worried +concentrate+irritable+anhedonia
      \n anhedonia ~~concentrate \n anhedonia~~guilty',
      negfisher[[11]], var.cov = T)


###8.17######
sim1smolall <- lapply(sim1smol, function(i) MIIVdel_new(model1, i, .05))
sim1all <- lapply(sim1, function(i) MIIVdel_new(model1, i, .05))
sim2smolall <- lapply(sim2smol, function(i) MIIVdel_new(model1, i, .05))
sim2all <- lapply(sim2, function(i) MIIVdel_new(model1, i, .05))
sim3smolall <- lapply(sim3smol, function(i) MIIVdel_new(model1, i, .05))
sim3all <- lapply(sim3, function(i) MIIVdel_new(model1, i, .05))
sim4smolall <- lapply(sim4smol, function(i) MIIVdel_new(model1, i, .05))
sim4all <- lapply(sim4, function(i) MIIVdel_new(model1, i, .05))

EFAmiive3(sim2smol[[1]])
