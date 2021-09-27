


# c('Ordinal', 'Beta', 'Gaussian', 'NegBinom', 'ZIP', 'ZINB', 'Multinom')

rm(list = ls())
gc()
#source('./sandbox/sim_data_types.R')
library(SPR)
library(glmmTMB)

# Check if we can create cross-sectional data:
library(ordinal)
sim.out <- SPR::sim_dat_types(N=5000,
                              data.type = 'Ordinal',
                              reg.formula =  formula(~Group),
                              subject.var = 0,
                              number.timepoints = 1)

dat <- sim.out$dat
modt1 <- ordinal::clm(ordered(Y_ord) ~ Group, data = dat)
summary(modt1)
# looks like we can create cross-sectional data just fine


#---
# ZINB
sim.out <- SPR::sim_dat_types(N=2000,
                         data.type = 'ZINB',
                         reg.formula = formula(~ Time))
  dat <- sim.out$dat
mod <- glmmTMB::glmmTMB(Y_zinb ~ Time + (1|USUBJID),
                ziformula = ~ 1,
                data = dat,
                family = nbinom2(link = "log"))
summary(mod)

#--------
# ZIP
sim.out <- SPR::sim_dat_types(N=2000,
                         data.type = 'ZIP',
                         reg.formula = formula(~ Time))
  dat <- sim.out$dat
mod <- glmmTMB(Y_zip ~ Time + (1|USUBJID),
                ziformula = ~ 1,
                data = dat,
                family = poisson)
summary(mod)

#--------------
library(lme4)
sim.out <- SPR::sim_dat_types(data.type = 'Gaussian')
  dat <- sim.out$dat
mod1 <- lme4::lmer(Y_gaussian ~ Group + Time + Group*Time + (1|USUBJID),
             data = dat)
summary(mod1)

#---
# Gaussian - aligns with lme4
sim.out <- SPR::sim_dat_types(data.type = 'Gaussian')
  dat <- sim.out$dat
mod0 <- glmmTMB(Y_gaussian ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat,
                family = gaussian)
summary(mod0)


#---
# Poisson
sim.out <- SPR::sim_dat_types(data.type = 'Poisson')
  dat <- sim.out$dat
mod1 <- glmmTMB(Y_pois ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat,
                family = poisson)
summary(mod1)

#---
# Negative Binomial
sim.out <- SPR::sim_dat_types(data.type = 'NegBinom')
  dat <- sim.out$dat
mod2 <- glmmTMB(Y_nb ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat,
                family = nbinom2(link = "log"))
summary(mod2)



#---
# Beta distribution
sim.out <- SPR::sim_dat_types(data.type = 'Beta', subject.var = 1,
                              reg.formula = formula(~ Group),
                              shape2 = 500, phi = 5)
  dat <- sim.out$dat
mod3 <- glmmTMB(Y_beta ~ Group  + (1|USUBJID),
                data = dat,
                family = beta_family(link = "logit"))
summary(mod3)
mu <- exp(1)/(1 + exp(1))
phi <- 5
vv <- mu*(1-mu)/(phi + 1)
vv


library(ordinal)
#---
# Ordinal
sim.out <- SPR::sim_dat_types(data.type = 'Ordinal')
  dat <- sim.out$dat
mod4 <- clmm(ordered(Y_ord) ~ Group + Time + Group*Time + (1|USUBJID),
             data = dat)
summary(mod4)
sim.out$thresholds
sim.out$Beta


#---
# Nominal data
sim.out <- sim_dat_types(N = 5000,
                         data.type = 'Multinom',
                         reg.formula =  formula(~ Time) )
  dat <- sim.out$dat

# Fit Models:
library(nnet)
dat$Y_nom_factor <- as.factor(dat$Y_nom)
mod <- nnet:::multinom(Y_nom_factor ~ Time, data = dat)
summary(mod)
t(coef(mod))
sim.out$Beta

mod5 <- clmm(ordered(Y_nom) ~ Group + Time + Group*Time + (1|USUBJID),
            # weights = Freq,
             data = dat)
mod6 <- clm2(Y_nom_factor ~ 1, nominal = ~ Time, data = dat)
summary(mod6)
