
# 6.16.21: we seem to be able to recapture the generating parameters using
# the conditional CLMM model, without any correction for "population average"
# like we did with the GLMM sim study. Why is that?
# Data generation is almost EXACTLY the same as the binomial!
# wtf!
# nothing to do with the drop-out
# I guess do the simulation study thing with an ordinal gee?
# see if this works to recapture the estimates?


rm(list = ls())
gc()

library(MASS)
library(polycor)
library(SPR)
###
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat.R")
#source("C:/Users/ciaconangelo/OneDrive - OPEN Health/Documents/SPR/sandbox/sim_dat_ord_logistic.R")
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/dropout.R")
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/VarCov_gls.R")

 set.seed(6162021)

  sim.out <- SPR::sim_dat_ord_logistic(N = 500,
                 number.groups = 2 ,
                 number.timepoints = 4,
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 0.5,
                 thresholds = c(0),
                 corr = 'ar1',
                 cor.value = 0.4)


  dat <- sim.out$dat
  str(dat)

# Check the correlations:
polychor(x = dat$Y_comp[dat$Time == 'Time_1'], y = dat$Y_comp[dat$Time == 'Time_2'])
polychor(x = dat$Y_comp[dat$Time == 'Time_2'], y = dat$Y_comp[dat$Time == 'Time_3'])
polychor(x = dat$Y_comp[dat$Time == 'Time_3'], y = dat$Y_comp[dat$Time == 'Time_4'])


library(ordinal)
mod.clm <- clm(as.factor(Y_comp) ~ Group + Time + Group*Time, data= dat)
summary(mod.clm)
out$Beta
out$thresholds
library(emmeans)
emmeans::emmeans(mod.clm, pairwise ~ Group | Time)

st <- Sys.time()
mod.clmm <- clmm(as.factor(Y_comp) ~ Group + Time + Group*Time + (1|USUBJID),
                 data= dat)
et <- Sys.time()
et - st
summary(mod.clmm)
mod.clm$coefficients
mod.clmm$coefficients

mod.polr <- MASS:::polr(as.factor(Y_comp) ~ Group + Time + Group*Time ,
                        data= dat, method = 'logistic', Hess = T)
summary(mod.polr)
# same estimates, good!

# maybe not correct?
sigma2 <- ordinal::VarCorr(mod.clmm)
sigma2 <- as.numeric(as.data.frame(sigma2))
denom.logit <- sqrt( (sigma2 + (pi^2/3))/(pi^2/3) )

#------
# Generalized Estimating Equations:
library(multgee)
dat.gee <- dat
tmp <- order(dat.gee$id.geepack)
dat.gee <- dat.gee[tmp, ]
mod.gee <- ordLORgee(formula = Y_comp ~ Time + Group + Time*Group,
                     data = dat.gee,
                     id = id.geepack,
                     repeated = Time, LORstr = "uniform")
summary(mod.gee)
coef(mod.gee)

# library(geepack)
# dat.gee <- dat
# tmp <- order(dat.gee$id.geepack)
# dat.gee <- dat.gee[tmp, ]
# dat.gee$Y_comp <- dat.gee$Y_comp + 1
# mod.gee1 <- ordgee(ordered(Y_comp) ~ Time + Group + Time*Group,
#                    id = id.geepack,
#                    data = dat.gee,
#                    corstr = "exchangeable")
# summary(mod.gee1)

cbind('Gen_par' = c(out$thresholds, as.vector(out$Beta[-1,])),
      'CLM' = coef(mod.clm),
     # 'CLM_polr' = coef(mod.polr),
      'ord_GEE' = coef(mod.gee),
      'CLMM' = coef(mod.clmm),
      'CLMM_pop_avg_logit' = coef(mod.clmm)/denom.logit)
