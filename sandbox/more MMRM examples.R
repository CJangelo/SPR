

rm(list = ls())
gc()

library(SPR)
library(MASS)
library(glmmTMB) # https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
library(emmeans)
library(nlme)
library(lme4)

  set.seed(03212021)

  sim.out <- sim_dat(N = 100,
                 number.groups = 2 ,
                 number.timepoints = 4,
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 0,
                 corr = 'ar1',
                 cor.value = 0.8,
                 var.values = 2)


  dat <- sim.out$dat
  dat <- dropout(dat = dat,
                 type_dropout  = c('mar'),
                 prop.miss = 0.3)


    mod.gls1 <- nlme::gls(Y_comp ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)

     mod.us1 <- glmmTMB::glmmTMB(Y_comp ~ Group + Time + Group*Time + us(Time + 0 | USUBJID),
                  data=dat,
                  REML = T,
                  dispformula=~0)

     mod.lme <- nlme::lme(Y_comp ~ Group + Time + Group*Time,
                           random =  ~ 0 +  Time | USUBJID, data = dat)

     summary(mod.lme)
     # Fixed effects:
     fixed.effects(mod.lme)
     glmmTMB::fixef(mod.us1)
     coef(mod.gls1)
     sim.out$Beta
     # Variance-covariance matrix:
     nlme::getVarCov(mod.lme, type = 'marginal')
     SPR::VarCov_gls(mod.gls1)
     out <- glmmTMB::VarCorr(mod.us1)
     out$cond$USUBJID[1:4, 1:4]
     #attributes(out$cond$USUBJID)
     sim.out$sigma

  emmeans::emmeans(mod.gls1, pairwise ~ Group | Time, mode = 'satt')
  emmeans::emmeans(mod.lme, pairwise ~ Group | Time, mode = 'satt')
  emmeans::emmeans(mod.lme, pairwise ~ Group | Time, mode = 'df.error')
