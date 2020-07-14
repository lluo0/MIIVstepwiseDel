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
