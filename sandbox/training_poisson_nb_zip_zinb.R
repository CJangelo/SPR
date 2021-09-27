
rm(list = ls())
gc()
#source('./sandbox/sim_data_types.R')
library(SPR)

N = 1e3
Beta <- matrix(c(.2, 1), nrow = 2)

#
sim.out <- SPR::sim_dat_types(N=N,
                              data.type = 'Poisson',
                              reg.formula =  formula(~Group),
                              Beta = Beta,
                              subject.var = 0,
                              number.timepoints = 1)

df <- sim.out$dat

sim.out <- SPR::sim_dat_types(N=N,
                              data.type = 'NegBinom',
                              reg.formula =  formula(~Group),
                              Beta = Beta,
                              subject.var = 0,
                              number.timepoints = 1)

df$Y_nb <- sim.out$dat$Y_nb


sim.out <- SPR::sim_dat_types(N=N,
                              data.type = 'ZIP',
                              reg.formula =  formula(~Group),
                              Beta = Beta,
                              subject.var = 0,
                              number.timepoints = 1)

df$Y_zip <- sim.out$dat$Y_zip


sim.out <- SPR::sim_dat_types(N=N,
                              data.type = 'ZINB',
                              reg.formula =  formula(~Group),
                              Beta = Beta,
                              subject.var = 0,
                              number.timepoints = 1)

df$Y_zinb <- sim.out$dat$Y_zinb


library(glmmTMB)
mod1 <- glmmTMB::glmmTMB(Y_pois ~ Group, family = poisson, data = df)
summary(mod1)
mod2 <- glmmTMB::glmmTMB(Y_nb ~ Group, family = nbinom2, data = df)
summary(mod2)
mod3 <- glmmTMB(Y_zip ~ Group,
                    ziformula = ~ 1,
                    data = df,
                    family = poisson)
summary(mod3)

mod4 <- glmmTMB(Y_zinb ~ Group,
                    ziformula = ~ 1,
                    data = df,
                    family = nbinom2(link = "log"))
summary(mod4)


df <- df[, c('USUBJID', 'Group', 'Y_pois', 'Y_nb', 'Y_zip', 'Y_zinb')]

colnames(df) <- c('USUBJID', 'Group', 'Y1', 'Y2', 'Y3', 'Y4')

write.table(df, na = '.', quote = F,
            sep = ', ', col.names = T,
            row.names = F, file = 'dataset4.txt')
