
# Conditional Model
# Generated data using conditional model
# Other CLMM Vignettes were using marginal data generation

# TODO: run a simulation study to check model parameter recovery;
# there looks to be something suspect going on with the Beta, all the
# itneractions ahve the same beta value...

rm(list = ls())
gc()

#library(MASS)
library(polycor)

###
source("C:/Users/ciaconangelo/OneDrive - OPEN Health/Documents/SPR/sandbox/sim_dat_ord_logistic_conditional.R")

 set.seed(6212021)
 # N = 500
 # number.groups = 2
 # number.timepoints = 4
 # reg.formula <-  formula(~ Group + Time + Time*Group)
 # Beta = 0.5
 # thresholds = c(-2, -1, 0, 1)
 # subject.var = 1
 # residual.var = 1
 # cond.mcar = F
 # Covariate = F

 out <- sim_dat_ord_logistic_conditional(N = 500,
                                         number.groups = 2 ,
                                         number.timepoints = 4,
                                         reg.formula = formula( ~ Time + Group + Time*Group),
                                         Beta = 0.5,
                                         thresholds = c( 0),
                                         subject.var = 2,
                                         cond.mcar = F,
                                         Covariate = F)


  dat <- out$dat
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

library(lme4)
mod.glmm <- glmer(Y_comp ~ Time + Group + Time*Group + (1|USUBJID),
                  family = 'binomial', data = dat)
summary(mod.glmm)
fixef(mod.glmm)

s2 <- lme4::VarCorr(mod.glmm)
as.data.frame(s2)
s2 <- ordinal::VarCorr(mod.clmm)
as.data.frame(s2)

