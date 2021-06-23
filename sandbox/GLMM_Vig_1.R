

# 5.28.21
# First step is to recover estimates with full data
# use GEE
# use GLMM and then convert back to population average estimate
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650580/#R11
# https://projecteuclid.org/download/pdf_1/euclid.ss/1009212671

# Then figure out a way to compute the standard errors of the population average
# unclear how to do that - that's only thing preventing you from abandoning
# the GEE!

rm(list = ls())
gc()

library(MASS)
library(SPR)
library(geepack)

###
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat.R")
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat_binom.R")
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/dropout.R")
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/VarCov_gls.R")

 set.seed(321)

  out <- SPR::sim_dat_binom(N = 1000,
                 number.groups = 2 ,
                 number.timepoints = 4,
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 1,
                 corr = 'ar1',
                 cor.value = 0.4)


  dat <- out$dat
  str(dat)

  # dat <- SPR::dropout(dat = dat,
  #                type_dropout  = c('mcar'),
  #                prop.miss = 0.3)
  #
# Binomial
#aggregate(Y_comp ~ Group + Time, FUN = function(x) mean(x, na.rm = T), dat = dat, na.action = na.pass)
mod.glm1 <- glm(Y_comp ~ Time + Group + Time*Group, family = 'binomial', data = dat)
  matrix(coef(mod.glm1), ncol = 1, dimnames = list(names(coef(mod.glm1))))
#mod.glm2 <- glm(Y_mcar ~ Time + Group + Time*Group, family = 'binomial', data = dat)
 # matrix(coef(mod.glm2), ncol = 1, dimnames = list(names(coef(mod.glm1))))
#out$Beta


# Generalized Estimating Equations:
#library(geepack)
###
mod.gee1 <- geeglm(Y_comp ~ Time + Group + Time*Group, id = id.geepack, data = dat, family = binomial, corstr = "ar1")
#summary(mod.gee1)
#matrix(coef(mod.gee1), ncol = 1, dimnames = list(names(coef(mod.gee1))))
#matrix(coef(mod.glm1), ncol = 1, dimnames = list(names(coef(mod.glm1))))

library(lme4)
mod.glmm <- glmer(Y_comp ~ Time + Group + Time*Group + (1|USUBJID),
                  family = 'binomial', data = dat)
summary(mod.glmm)
sigma2 <- lme4::VarCorr(mod.glmm)
sigma2 <- as.data.frame(sigma2)$vcov
denom.logit <- sqrt( (sigma2 + (pi^2/3))/(pi^2/3) )
denom.probit <- sqrt( 1 + sigma2)
#fixef(mod.glmm)/denom

cbind('Gen_par' = out$Beta,
      'GLM' = coef(mod.glm1),
      'GEE' = coef(mod.gee1),
      'GLMM' = fixef(mod.glmm),
      'GLMM_pop_avg_logit' = fixef(mod.glmm)/denom.logit,
      'GLMM_pop_avg_probit' = fixef(mod.glmm)/denom.probit)


#
# This will hang R
mod.gee2 <- geeglm(Y_mcar ~ Time + Group + Time*Group, id = USUBJID, data = dat, family = binomial, corstr = "ar1")
summary(mod.gee2)
# # As soon as there's missing data, these don't line up - pretty close though
matrix(coef(mod.gee2), ncol = 1, dimnames = list(names(coef(mod.gee2))))
matrix(coef(mod.glm2), ncol = 1, dimnames = list(names(coef(mod.glm1))))


# ------------------------------------------------------------------------------
# Conditional MCAR:

  out <- sim_dat_binom(N = 10000,
                 number.groups = 2 ,
                 number.timepoints = 4,
                 cond.mcar = T,
                 #reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 1,
                 corr = 'ar1',
                 cor.value = 0.4)


 dat <- out$dat
  str(dat)

  dat <- dropout(dat = dat,
                 type_dropout  = c('mcar'),
                 prop.miss = 0.3)

  out$reg.formula
# Binomial
mod.glm1 <- glm(as.factor(Y_comp) ~ Group + Time + Covariate + Time * Group + Covariate * Time, family = 'binomial', data = dat)
  matrix(coef(mod.glm1), ncol = 1, dimnames = list(names(coef(mod.glm1))))
mod.glm2 <- glm(Y_mcar ~ Group + Time + Covariate + Time * Group + Covariate * Time, family = 'binomial', data = dat)
  matrix(coef(mod.glm2), ncol = 1, dimnames = list(names(coef(mod.glm2))))
  out$Beta

# Misspecified Model:
mod.glm3 <- glm(Y_comp ~ Group + Time + Time * Group, family = 'binomial', data = dat)
  matrix(coef(mod.glm3), ncol = 1, dimnames = list(names(coef(mod.glm3))))
  out$Beta
# Way off, as expected
