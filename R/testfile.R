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
