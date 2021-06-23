
rm(list = ls())
gc()

library(MASS)
library(polycor)

###
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat_ord.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/dropout.R")
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/VarCov_gls.R")

 set.seed(321)

  out <- sim_dat_ord(N = 10000,
                 number.groups = 2 ,
                 number.timepoints = 4,
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 0.5,
                 thresholds = c(0.25, 0.50, 0.75),
                 corr = 'ar1',
                 cor.value = 0.4)


  dat <- out$dat
  str(dat)

  dat <- dropout(dat = dat,
                 type_dropout  = c('mcar'),
                 prop.miss = 0.3)


cor(x = dat$Y_comp[dat$Time == 'Time_1'], y = dat$Y_comp[dat$Time == 'Time_2'])
cor(x = dat$Y_comp[dat$Time == 'Time_2'], y = dat$Y_comp[dat$Time == 'Time_3'])
cor(x = dat$Y_comp[dat$Time == 'Time_3'], y = dat$Y_comp[dat$Time == 'Time_4'])

aggregate(Y_comp ~ Group + Time, FUN = function(x) mean(x, na.rm = T), dat = dat, na.action = na.pass)
aggregate(Y_comp ~ Time + Group, FUN = function(x) prop.table(table(x)), dat = dat, na.action = na.pass)


polychor(x = dat$Y_comp[dat$Time == 'Time_1'], y = dat$Y_comp[dat$Time == 'Time_2'])
#polychor(x = as.factor(dat$Y_comp[dat$Time == 'Time_1']), y = as.factor(dat$Y_comp[dat$Time == 'Time_2']))
polychor(x = dat$Y_comp[dat$Time == 'Time_2'], y = dat$Y_comp[dat$Time == 'Time_3'])
polychor(x = dat$Y_comp[dat$Time == 'Time_3'], y = dat$Y_comp[dat$Time == 'Time_4'])



# Ordinal
mod.ord1 <- MASS:::polr(as.factor(Y_comp) ~ Group + Time + Group*Time , data= dat, method = 'probit', Hess = T)
mod.ord2 <- MASS:::polr(as.factor(Y_mcar) ~ Group + Time + Group*Time , data= dat, method = 'probit', Hess = T)

out$Beta
matrix(mod.ord1$coefficients, ncol = 1, dimnames = list(names(mod.ord1$coefficients)))
matrix(mod.ord2$coefficients, ncol = 1, dimnames = list(names(mod.ord2$coefficients)))
out$thresholds
mod.ord1$zeta
mod.ord2$zeta

library(ordinal)
mod.clmm <- clmm(as.factor(Y_comp) ~ Group + Time + Group*Time + (1|USUBJID),
                 data= dat)
summary(mod.clmm)

