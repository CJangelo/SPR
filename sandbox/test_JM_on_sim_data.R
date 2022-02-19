

# 11.24.2021
# Simulate data for the Joint Model
# This is the PRO data part
# output the mt

# Not able to break the LMM estimation of PRO scores across groups
# the drop-out across Groups is not very different
# the LMM can recover the gen parameters just fine even though the survival
# process relies on the LMM process
# TODO Break it then figure out how to model it correctly



#---
rm(list = ls())
gc()


#source("C:/Users/ciaconangelo/OneDrive - OPEN Health/Documents/SPR/sandbox/sim_JM_PRO.R")
source("C:/Users/ciaconangelo/OneDrive - OPEN Health/Documents/SPR/sandbox/sim_JM_PRO_random_slopes.R")
source("C:/Users/ciaconangelo/OneDrive - OPEN Health/Documents/SPR/sandbox/sim_JM_surv.R")

# Generate Data
set.seed(11302021)
out <- sim_dat_PRO_random_slopes(N = 1e3,
                   number.groups = 2,
                   number.timepoints = 4,
                   reg.formula = formula( ~ Group*Time),
                   sigma2 = 1)
dat <- out$dat

dat[dat$USUBJID == 'Subject_0001', ]

# TEST/CHECK: Fit Model to PRO data
library(nlme)
test <- lme(Y ~ Group*Time, random =  ~ Time | USUBJID, data = dat)
summary(test)
fixed.effects(test)
out$Beta
# g2g

#---------
# Generate survival data
out <- sim_dat_surv(dat.PRO = dat, reg.formula = formula(~ Group))
  dat.surv <- out$dat.surv
  dat.surv <- dat.surv[, c('USUBJID', 'Group', 'time.surv', 'status.surv')]


#------------------------------
# Drop PRO data after survival process event:
  dat$PRO <- dat$Y
  id <- unique(dat$USUBJID)
for (i in id) {
  max.time <- dat.surv[dat.surv$USUBJID == i, 'time.surv']
  drop <- dat[dat$USUBJID == i, 'Time'] > max.time
  dat[dat$USUBJID == i, 'PRO'][drop] <- NA
}


#-------------------
# Analyze - confirm can recover gen param
library(survival)
mod.surv <- coxph(Surv(time.surv, status.surv) ~ Group, data = dat.surv, x = T)
summary(mod.surv)
coef(mod.surv)


library(nlme)
aggregate(Y ~ Group + Time, FUN = mean, data = dat) # complete data
aggregate(PRO ~ Group + Time, FUN = mean, data = dat) # with missing
aggregate(PRO ~ Group + Time, FUN = function(x) mean(is.na(x)), data = dat, na.action = na.pass)

dat.PRO <- dat[complete.cases(dat), ]
dat.PRO <- dat.PRO[order(dat.PRO$USUBJID), ]
mod.lme <- lme(PRO ~ Group*Time, random =  ~ Time | USUBJID, data = dat.PRO)


library(JM)
st <- Sys.time()
mod.JM <- JM::jointModel(mod.lme, mod.surv,
                            timeVar = 'Time',
                            method = 'Cox-PH-GH')
et <- Sys.time()
print(et - st)

fixed.effects(test)
fixed.effects(mod.lme)
coef(summary(mod.JM))$Long[,'Value']
