

# Simulate longitudinal data
# Use complete data only (for now)
# Compare the following:
# means, conditional model, marginal model

# 9.29.21 - this should be next vignette
# unclear how useful this is - just says that the random intercept makes it
# so that the model marginal means do not align with the descriptive means

# Review slides on Marginal vs Conditional
# slide 184 here - http://www.drizopoulos.com/courses/EMC/CE08.pdf
# Run through this R code to establish the pattern
# want everyone to be aware of this issue
# Tie it back to "How to walk sponsor through statistical evidence?"
# You are just fundamentally estimating something different w a conditional model
# The odds ratios are just different, the rates are different
# You are estimating a time trend conditional on the subject
# Frank Harrell comments are basically that this is personalized medicine
# https://www.fharrell.com/post/robcov/

# Review slides on marg vs conditional
# R code
# maybe pass dataset for a count distribution to show this further
# Review the F Harrell article or commentary - was it his blog?

rm(list = ls())
gc()
library(SPR)
library(glmmTMB)
set.seed(9292021)



#---
# continuous data with gaussian/normal error distribution
#
sim.out <- SPR::sim_dat_types(N=1000,
                              data.type = 'Gaussian',
                              reg.formula =  formula(~ Time),
                              subject.var = 1,
                              number.timepoints = 4)
dat <- sim.out$dat
aggregate(Y_gaussian ~ Time, FUN = mean, data = dat)
modc <- glmmTMB::glmmTMB(Y_gaussian ~ Time + (1|USUBJID), data = dat)
as.data.frame(emmeans::emmeans(modc, ~ Time))

modm <- glmmTMB::glmmTMB(Y_gaussian ~ Time + us(Time + 0|USUBJID),
                         REML = T,
                         dispformula=~0,
                         data = dat)
as.data.frame(emmeans::emmeans(modm, ~ Time))
# means align with the model values, both conditional model and marginal model


# Simulate MMRM marginal data
  gen.beta <- sim.out$Beta
sim.out <- sim_dat(N = 1000,
                 number.groups = 2 ,
                 number.timepoints = 4,
                 reg.formula = formula( ~ Time),
                 Beta = gen.beta,
                 corr = 'ar1',
                 cor.value = 0.8,
                 var.values = 2)

sim.out$Beta
sim.out$cor.mat
dat <- sim.out$dat
str(dat)
aggregate(Y_comp ~ Time, FUN = mean, data = dat)
modc <- glmmTMB::glmmTMB(Y_comp ~ Time + (1|USUBJID), data = dat)
as.data.frame(emmeans::emmeans(modc, ~ Time))

modm <- glmmTMB::glmmTMB(Y_comp ~ Time + us(Time + 0|USUBJID),
                         REML = T,
                         dispformula=~0,
                         data = dat)
as.data.frame(emmeans::emmeans(modm, ~ Time))

# okay at least with the AR autocorrelation structure and complete data
# everything aligns. Can probably not include.

#-----------------------------
# Binary data
# conditional and marginal models do not align
# this is due to the logistic link

sim.out <- SPR::sim_dat_types(N=1000,
                              data.type = 'Binomial',
                              reg.formula =  formula(~Time),
                              subject.var = 3,
                              number.timepoints = 4)
dat <- sim.out$dat
aggregate(Y_binom ~ Time, FUN = mean, data = dat)
modc <- glmmTMB::glmmTMB(Y_binom ~ Time + (1|USUBJID), data = dat, family = binomial)
as.data.frame(emmeans::emmeans(modc, ~ Time, type = 'response'))
exp(as.data.frame(emmeans::emmeans(modc, ~ Time))$emmean)
ptab <- prop.table(xtabs(~ Y_binom + Time, data = dat), margin = 2)
  ptab[2,]/(1 - ptab[2,])


#-----------
# GEE - general estimating equation
library(geepack)
mod.gee <- geepack::geeglm(Y_binom ~ Time,
                              id = id.geepack,
                              data = dat,
                              family = binomial, corstr = "ar1")

as.data.frame(emmeans::emmeans(mod.gee, ~ Time, type = 'response'))
exp(as.data.frame(emmeans::emmeans(mod.gee, ~ Time))$emmean)
# Conclusion: marginal model (GEE) aligns exactly with the means
# conditional model does not
# (i think the conditional model is unusually close here, typically it's not)


#----------
# Ordinal data
# skip this - marginal means too complicated
# https://github.com/paul-buerkner/brms/issues/980
# https://github.com/rvlenth/emmeans/blob/master/R/ordinal-support.R

sim.out <- SPR::sim_dat_types(N=1000,
                              data.type = 'Ordinal',
                              reg.formula =  formula(~Time),
                              subject.var = 1,
                              number.timepoints = 4)
dat <- sim.out$dat
ptab <- prop.table(xtabs(~ Y_ord + Time, data = dat), margin = 2)
ptab['3', ]/(1 - ptab['3', ])

# Conditional Model
library(ordinal)
modc <- ordinal::clmm(as.factor(Y_ord) ~ Time + (1|USUBJID), nAGQ = 5, data= dat)
exp(coef(modc))
exp(as.data.frame(emmeans::emmeans(modc, ~ Time))$emmean)

# Marginal model
# Generalized Estimating Equations:
library(multgee)
dat.gee <- dat
tmp <- order(dat.gee$id.geepack)
  dat.gee <- dat.gee[tmp, ]
mod.gee <- multgee::ordLORgee(formula = Y_ord ~ Time ,
                                data = dat.gee,
                                id = id.geepack,
                                repeated = Time, LORstr = "uniform")
as.data.frame(emmeans::emmeans(mod.gee, ~ Time))
# not supported
exp(coef(mod.gee))

# ignore the covariance structure
mod.ig <- ordinal::clm(as.factor(Y_ord) ~ Time, data= dat)
exp(as.data.frame(emmeans::emmeans(mod.ig, ~ Time))$emmean)
exp(coef(mod.ig))


#---
# Poisson
sim.out <- SPR::sim_dat_types(N=1000,
                              data.type = 'Poisson',
                              reg.formula =  formula(~Time),
                              subject.var = 1,
                              number.timepoints = 4)
dat <- sim.out$dat
aggregate(Y_pois ~ Time, FUN = mean, data = dat)
modc <- glmmTMB(Y_pois ~ Time + (1|USUBJID), data = dat, family = poisson)
as.data.frame(emmeans::emmeans(modc, ~ Time, type = 'response'))
# means and conditional model do not align
# does not matter if you're using the 'marginal means"
# you're fitting a conditional model

mod.ig <- glmmTMB(Y_pois ~ Time, data = dat, family = poisson)
as.data.frame(emmeans::emmeans(mod.ig, ~ Time, type = 'response'))

# the difference comes from the random effect
# simulate it with zero variance
sim.out2 <- SPR::sim_dat_types(N=1000,
                              data.type = 'Poisson',
                              reg.formula =  formula(~Time),
                              subject.var = 0,
                              number.timepoints = 4)

dat2 <- sim.out2$dat
aggregate(Y_pois ~ Time, FUN = mean, data = dat2)
modc2 <- glmmTMB(Y_pois ~ Time + (1|USUBJID), data = dat2, family = poisson)
as.data.frame(emmeans::emmeans(modc2, ~ Time, type = 'response'))

# notice how the estimates from the first conditional model align
# pretty closely with the second set of averages
# with subject variance = 1:
as.data.frame(emmeans::emmeans(modc, ~ Time, type = 'response'))
# with subject variance = 0:
aggregate(Y_pois ~ Time, FUN = mean, data = dat2)



#-----------------
# Use to simulate a dataset - ask them to fit it
# This is basically same as Poisson
# Compute two negative binomial variables
# one with var = 3, another with var = 0
# ask them to fit means and conditional models
# describe the subject-level effects observed in each

sim.out <- SPR::sim_dat_types(N=1000,
                              data.type = 'NegBinom',
                              reg.formula =  formula(~Time),
                              subject.var = 1,
                              number.timepoints = 4)

dat <- sim.out$dat
aggregate(Y_nb ~ Time, FUN = mean, data = dat)
modc <- glmmTMB(Y_nb ~ Time + (1|USUBJID),
                data = dat,
                family = nbinom2(link = "log"))

as.data.frame(emmeans::emmeans(modc, ~ Time, type = 'response'))

