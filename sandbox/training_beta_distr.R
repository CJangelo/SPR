
rm(list = ls())
gc()
#source('./sandbox/sim_data_types.R')
library(SPR)

N = 200
Beta <- matrix(c(.2, .5), nrow = 2)
set.seed(8122021)
#
sim.out <- SPR::sim_dat_types(N=N,
                              data.type = 'Beta',
                              reg.formula =  formula(~Group),
                              Beta = Beta,
                              subject.var = 0,
                              number.timepoints = 1)

df <- sim.out$dat

library(glmmTMB)

mod1 <- glmmTMB::glmmTMB(Y_beta ~ Group, 
                         family = beta_family(link = "logit"), 
                         data = df)
summary(mod1)

library(betareg)
mod2 <- betareg::betareg(Y_beta ~ Group, data = df)
summary(mod2)

df <- df[, c('USUBJID', 'Group', 'Y_beta')]

write.table(df, na = '.', quote = F,
            sep = ', ', col.names = T,
            row.names = F, file = 'dataset5.txt')
