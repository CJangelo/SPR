


# 9.9.21:
# Courage data project

# This is a check - can we just take the inverse of the binary data?
# Shouldn't the point estimates be 1-Beta?
# shoudln't the null hypothesis be P=0 instead of P=1?
# This suggests that yes it should
# Unclear why it differed with empirical data


rm(list = ls())
gc()
#source('./sandbox/sim_data_types.R')
library(SPR)
library(glmmTMB)


#---
# Poisson
sim.out <- SPR::sim_dat_types(data.type = 'Ordinal', threshold = c(0) )
dat <- sim.out$dat

# Check the data:
str(dat)
table(dat$Y_ord)
hist(dat$Y_ord)
aggregate(Y_ord ~ Group*Time, mean, data = dat)

# Fit the models:
mod1 <- glmmTMB(Y_ord ~ Group + Time + Group*Time + us(1|USUBJID),
                data = dat,
                family = poisson)
summary(mod1)

# Now invert the data:
dat$Y_inv <- 1 - dat$Y_ord
xtabs(~ Y_inv + Y_ord, data = dat)
# okay now model this:

mod2 <- glmmTMB(Y_inv ~ Group + Time + Group*Time + us(1|USUBJID),
                data = dat,
                family = poisson)
summary(mod2)

# easier to look at the marginal means:
library(emmeans)
emm1 <- emmeans(mod1, ~ Group | Time, type = 'response')
emm2 <- emmeans(mod2, ~ Group | Time, type = 'response')

as.data.frame(emm1)[7:8, ]
as.data.frame(emm2)[7:8, ]

# Was it because there was drop-out?
dat$Y_comp <- dat$Y_ord
dat2 <- SPR::dropout(dat, type_dropout = c('mcar', 'mar', 'mnar'), prop.miss = 0.5)
str(dat2)


# MCAR drop-out:
mod3 <- glmmTMB(Y_mcar ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat2,
                family = poisson)
emm3 <- emmeans(mod3, ~ Group | Time, type = 'response')
as.data.frame(emm3)[7:8, ]
# Compare to the complete data:
as.data.frame(emm1)[7:8, ]

# Invert data:
dat2$Y_mcar_inv <- 1 - dat2$Y_mcar
mod4 <- glmmTMB(Y_mcar_inv ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat2,
                family = poisson)
emm4 <- emmeans(mod4, ~ Group | Time, type = 'response')
as.data.frame(emm4)[7:8, ]
1- as.data.frame(emm4)[7:8,'rate' ]
as.data.frame(emm3)[7:8, ]
# Works fine! No issue here!

# MAR drop-out
mod5 <- glmmTMB(Y_mar ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat2,
                family = poisson)
emm5 <- emmeans(mod5, ~ Group | Time, type = 'response')
as.data.frame(emm5)[7:8, ]
# Compare to the complete data:
as.data.frame(emm1)[7:8, ]

# Invert data:
dat2$Y_mar_inv <- 1 - dat2$Y_mar
mod6 <- glmmTMB(Y_mar_inv ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat2,
                family = poisson)
emm6 <- emmeans(mod6, ~ Group | Time, type = 'response')
as.data.frame(emm6)[7:8, ]
1- as.data.frame(emm5)[7:8,'rate' ]


# Okay, one more time but with MNAR drop-out?
mod7 <- glmmTMB(Y_mnar ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat2,
                family = poisson)
emm7 <- emmeans(mod7, ~ Group | Time, type = 'response')
as.data.frame(emm7)[7:8, ]
# Compare to the complete data:
as.data.frame(emm1)[7:8, ]

# Invert data:
dat2$Y_mnar_inv <- 1 - dat2$Y_mnar
mod8 <- glmmTMB(Y_mnar_inv ~ Group + Time + Group*Time + (1|USUBJID),
                data = dat2,
                family = poisson)
emm8 <- emmeans(mod8, ~ Group | Time, type = 'response')
as.data.frame(emm8)[7:8, ]
1- as.data.frame(emm7)[7:8,'rate' ]
