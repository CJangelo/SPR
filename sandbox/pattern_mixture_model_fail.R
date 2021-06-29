

rm(list = ls())
gc()

library(SPR)
library(MASS)
library(glmmTMB) # https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
library(emmeans)
library(nlme)
library(lme4)



  sim.out <- SPR::sim_dat(N = 2000,
                 number.groups = 2 ,
                 number.timepoints = 4,
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 1,
                 corr = 'ar1',
                 cor.value = 0.8,
                 var.values = 1)

  sim.out$Beta
  dat <- sim.out$dat
  dat <- SPR::dropout(dat = dat,
                 type_dropout  = c('mcar', 'mar', 'mnar'),
                 prop.miss = 0.3)


  score <- 'Y_mnar'

  f1 <- as.formula(paste0(score, '~ Group + Time + Group*Time'))

#   MMRM
 # library(nlme)
  mod.gls1 <- gls( f1,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)

# pattern mixture model
  for (i in unique(dat$USUBJID)){

    dat[dat$USUBJID == i,'cohort'] <-
      paste0(1*is.na(dat[dat$USUBJID == i, score ]), collapse = '')

    }

    f2 <- as.formula(paste0(score,
                            '~ Group + Time + Group*Time + cohort + Group*cohort'))

  mod.gls2 <- gls( f2,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)



  # Compute marginal means
  #library(emmeans)
  mod.emms.gls1 <- emmeans(mod.gls1,
                           pairwise ~ Group | Time,
                           adjust = 'none',
                           mode = 'df.error')
  mod.emms.gls2 <- emmeans(mod.gls2,
                           pairwise ~ Group + cohort| Time,
                           adjust = 'none',
                           mode = 'df.error')

  mod.emms.gls1$contrasts
  mod.emms.gls2$contrasts

 Group_1 0000 - Group_2 0000 -0.69118
 Group_1 0001 - Group_2 0001 -0.64531
  Group_1 0011 - Group_2 0011 -0.71364
 Group_1 0111 - Group_2 0111 -1.11789

aggregate(Y_mnar ~ Time, FUN = function(x) mean(is.na(x)), na.action = na.pass,data = dat)
xtabs(~ cohort + Time, data = dat)
  -0.69118*.7 +
  -0.64531*.1 +
  -0.71364*.1 +
  -1.11789*.1
