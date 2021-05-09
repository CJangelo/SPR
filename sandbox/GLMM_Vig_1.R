
rm(list = ls())
gc()

library(MASS)

###
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat_binom.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/dropout.R")
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/VarCov_gls.R")

 set.seed(321)
  
  out <- sim_dat_binom(N = 100, 
                 number.groups = 2 , 
                 number.timepoints = 4, 
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 1,
                 corr = 'ar1', 
                 cor.value = 0.4) 

  
  dat <- out$dat
  str(dat)
  
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mcar'), 
                 prop.miss = 0.3)
  
# Binomial
aggregate(Y_comp ~ Group + Time, FUN = function(x) mean(x, na.rm = T), dat = dat, na.action = na.pass)
mod.glm1 <- glm(as.factor(Y_comp) ~ Time + Group + Time*Group, family = 'binomial', data = dat)
  matrix(coef(mod.glm1), ncol = 1, dimnames = list(names(coef(mod.glm1))))
mod.glm2 <- glm(Y_mcar ~ Time + Group + Time*Group, family = 'binomial', data = dat)
  matrix(coef(mod.glm2), ncol = 1, dimnames = list(names(coef(mod.glm1))))
out$Beta


# Generalized Estimating Equations:
library(geepack)
###
mod.gee1 <- geeglm(Y_comp ~ Time + Group + Time*Group, id = id.geepack, data = dat, family = binomial, corstr = "ar1")
summary(mod.gee1)
matrix(coef(mod.gee1), ncol = 1, dimnames = list(names(coef(mod.gee1))))
matrix(coef(mod.glm1), ncol = 1, dimnames = list(names(coef(mod.glm1))))

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
  