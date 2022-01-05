



rm(list = ls())
gc()
library(SPR)
library(glmmTMB)
set.seed(9282021)

# zero out the random intercept
# Beta distribution
sim.out1 <- SPR::sim_dat_types(N=1000,
                              data.type = 'Beta',
                              number.timepoints = 3,
                              subject.var = 0,
                              reg.formula = formula(~ Group))
  dat1 <- sim.out1$dat

mod1 <- glmmTMB(Y_beta ~ Group + (1|USUBJID),
                data = dat1,
                family = beta_family(link = "logit"))
summary(mod1)
unlist(fixef(mod1))
xb <- c('Group_1' = -0.01497992, 'Group_2' = 1.05174065)
mu <- exp(xb)/(1 + exp(xb))
sim.out1$phi
phi <- 10
vv <- mu*(1-mu)/(phi + 1)
vv
aggregate(Y_beta ~ Group, var, data = dat1)
# model variance and observed variance line up

N = 1000
sim.out2 <- SPR::sim_dat_types(N=N,
                              data.type = 'Beta',
                              number.timepoints = 3,
                              subject.var = 1,
                              reg.formula = formula(~ Group))
dat2 <- sim.out2$dat
sim.out2$Beta
# dat2$Y_beta <- ifelse(dat2$Y_beta > 0.999, 0.999,
#                       ifelse(dat2$Y_beta < 0.001, 0.001, dat2$Y_beta))
mod2 <- glmmTMB(Y_beta ~ Group + (1|USUBJID),
                data = dat2,
                family = beta_family(link = "logit"))
summary(mod2)
aggregate(Y_beta ~ Group, var, data = dat2)
unlist(fixef(mod2))
xb <- c('Group_1' = 0, 'Group_2' = 1)
mu <- exp(xb)/(1 + exp(xb))
sim.out2$phi
phi <- 10
vv <- mu*(1-mu)/(phi + 1) # variance (call it residual variance)
vv # use estimates, get approximately this variance for each group
aggregate(Y_beta ~ Group, var, data = dat2) # observed total variance
# add vv and the variance of the predicted values from the model
# sums to the total observed variance
# var(mod.mu.i0)+mean(mod.mu.i0*(1-mod.mu.i0)/(1+mod.phi0))
# "mod.mu.i0" are the predicted values
# mod.mu.i0 <- exp(predict(mod0))/(1+exp(predict(mod0)))


# so the ICC (or rho) would be:
denom <- aggregate(Y_beta ~ Group, var, data = dat2)[,'Y_beta']
1 - vv/denom

