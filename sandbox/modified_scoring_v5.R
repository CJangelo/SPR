
# Next step is to pull the empirical data. If that actually works,
# you're in business. It needs to look like this.
# Then you can figure out how to fit the longitudinal ordinal model
#   https://hbiostat.org/proj/covid19/

#  https://hbiostat.org/R/rms/blrm.html

# 2.8.2022 -
# 1. Re-run the CLMM with the baseline value as a predictor
# 2. Adjust power accordingly - can bump down power of MMRM
# Goal is to have a simulation study where the MMRM sucks and the rank sum
# and the PO model are strong.
# 3. Fit the ordinal GEE instead of the CLMM, see if the OR better aligns
# with the empirical C-index.
# 4. alternative is to see if Hedeker has the ordinal thing available somewhere
# sim study: rank sum test, PO model, both with and without missing data
# Ideally we would have the empirical data and show an example, then get
# into the expected outcomes part
# Methods: explain it,
# sim study to evaluate power & type I error,
# empirical example


# bump up the power of the survival analysis to 0.8
# then see why the MMRM is working so well
# see if transofmred conditional to marginal is closer to empirical C-index




# 2.7.2022 - Modified Scoring version 3
# Assign a score of "4" at time of drop-out AND ALL SUBSEQUENT TIMEPOINTS!
# TO DO: simulate survival data so drop-out happens earlier than timepoint 1!

#-------------------------------------------------------------------
rm(list = ls())
gc()

source('./sandbox/sim_surv_data.R')
source('./sandbox/sim_PRO_data.R')
library(ordinal)
library(emmeans)
library(nlme)

c.index <- vector()
pval <- vector()
pval2 <- vector()
pval3 <- vector()
repl <- 7 # 11, 14, 15, 17, 18, 24


for (repl in 1:30){

  set.seed(02082022 + repl)
#----------------------------------
# Simulate Data:
sim.out <- sim_PRO_data()
dat.PRO <- sim.out$dat.PRO
dat.surv <- sim.out$dat.surv

#------------------------------
# Drop PRO data after survival process event:
  dat.PRO$PRO <- dat.PRO$Y_comp
  id <- unique(dat.PRO$USUBJID)
for (i in id) {
  max.time <- dat.surv[dat.surv$USUBJID == i, 'time.surv']
  drop <- dat.PRO[dat.PRO$USUBJID == i, 'Time'] > max.time
  dat.PRO[dat.PRO$USUBJID == i, 'PRO'][drop] <- NA
}


#-------------------
  # Analyze - confirm can recover gen param
  library(survival)
  mod.surv <- coxph(Surv(time.surv, status.surv) ~ Group, data = dat.surv, x = T)
  km_fit <- survfit(Surv(time.surv, status.surv) ~ Group, data=dat.surv)
  library(survminer)
  #
  #ggsurvplot(km_fit)
  pval3 <- rbind(pval3, 'Cox PH Reg' = summary(mod.surv)$waldtest['pvalue'])
#----------------------------------------------

  # Modified Scoring Algorithm: make drop-out the highest/worst score:
    #unique(as.factor(dat.PRO$PRO))
    dat.PRO$PRO_modified <- dat.PRO$PRO
    dat.PRO$PRO_modified[is.na(dat.PRO$PRO_modified)] <- 4

#-------------------------------------
  # MMRM - typical approach!
  dat.mmrm <- dat.PRO
  dat.mmrm$trans <-
    100*(dat.mmrm$PRO)/(diff(range(dat.mmrm$PRO, na.rm = T)))
  tmp <- dat.mmrm[dat.mmrm$Visit == 'Visit_1', c('USUBJID', 'trans')]
    colnames(tmp) <- c('USUBJID', 'trans_bl')
  dat.mmrm <- merge(x = dat.mmrm, y = tmp, by = 'USUBJID')
  dat.mmrm$delta <- dat.mmrm$trans - dat.mmrm$trans_bl

  mod.mmrm <- glmmTMB::glmmTMB(
    delta ~ Group + Visit + Group*Visit + us(Time + 0 | USUBJID),
    data = dat.mmrm[dat.mmrm$Visit != 'Visit_1', ],
    REML = T, dispformula=~0)

  emm.mmrm <- emmeans(mod.mmrm, pairwise ~ Group|Visit, adjust = 'none',  mode = 'df.error')

#------------------------------------------------------------------------

#------------
  # CLM
  mod.clm <- ordinal::clm(as.factor(Y_comp) ~ Group + Visit + Group*Visit,
                          nAGQ = 5, data= dat.PRO)

  mod.clm.missing <- ordinal::clm(as.factor(PRO) ~ Group + Visit + Group*Visit,
                          nAGQ = 5, data= dat.PRO)

  mod.clm.modified <- ordinal::clm(as.factor(PRO_modified) ~ Group + Visit + Group*Visit,
                          nAGQ = 5, data= dat.PRO)

  mod.mixed.modified <-
    ordinal::clmm(as.factor(PRO_modified) ~ Group + Visit + Group*Visit + (1|USUBJID),
                          nAGQ = 5, data= dat.PRO)

  # Marginal Means
  emm.clm <- emmeans::emmeans(mod.clm, pairwise ~ Group | Visit)
  emm.clm.missing <- emmeans::emmeans(mod.clm.missing, pairwise ~ Group | Visit)
  emm.clm.modified <- emmeans::emmeans(mod.clm.modified, pairwise ~ Group | Visit)
  emm.mixed.modified <- emmeans::emmeans(mod.mixed.modified, pairwise ~ Group | Visit)



  pval <- rbind(pval, c(
    'mmrm' = as.data.frame(emm.mmrm$contrasts)['3', 'p.value'], # Timepoint 3 bc no baseline!
    'PO_complete' = as.data.frame(emm.clm$contrasts)['4', 'p.value'],
    'PO_missing' = as.data.frame(emm.clm.missing$contrasts)['4', 'p.value'],
    'PO_modified' = as.data.frame(emm.clm.modified$contrasts)['4', 'p.value'],
    'PO_mixed_modified' = as.data.frame(emm.mixed.modified$contrasts)['4', 'p.value']))
 #
 #
  emm.clm <- as.data.frame(emm.clm$contrasts)['4', 'estimate']
  emm.clm.missing <- as.data.frame(emm.clm.missing$contrasts)['4', 'estimate']
  emm.clm.modified <- as.data.frame(emm.clm.modified$contrasts)['4', 'estimate']
  #emm.mixed.modified <- as.data.frame(emm.mixed.modified$contrasts)['4', 'estimate']

  # Transform from conditional to marginal, see if this is closer to C-index
   sigma2 <- unlist(VarCorr.clmm(mod.mixed.modified))
   denom.logit <- sqrt( (sigma2 + (pi^2/3))/(pi^2/3) )

  emm.mixed.modified2 <-
    as.data.frame(emm.mixed.modified$contrasts)['4', 'estimate']/denom.logit

  emm.mixed.modified <- as.data.frame(emm.mixed.modified$contrasts)['4', 'estimate']

#------------------------------------------------------------
  # Concordance Index
  g1 <- dat.PRO[dat.PRO$Visit == 'Visit_4' & dat.PRO$Group == 'Group_1', 'Y_comp']
  g2 <- dat.PRO[dat.PRO$Visit == 'Visit_4' & dat.PRO$Group == 'Group_2', 'Y_comp']
    wt.comp <- wilcox.test(x = g1, y = g2, alternative = "two.sided")
  g1.re <- sample(x = g1, size = 1e4, replace = T)
  g2.re <- sample(x = g2, size = 1e4, replace = T)

  C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
  C_emp_complete <- mean(C_emp, na.rm = T)

#------------------------------------------------------------
  # Concordance Index
  g1 <- dat.PRO[dat.PRO$Visit == 'Visit_4' & dat.PRO$Group == 'Group_1', 'PRO']
  g2 <- dat.PRO[dat.PRO$Visit == 'Visit_4' & dat.PRO$Group == 'Group_2', 'PRO']
      wt.missing <- wilcox.test(x = g1, y = g2, alternative = "two.sided")

  g1.re <- sample(x = g1, size = 1e4, replace = T)
  g2.re <- sample(x = g2, size = 1e4, replace = T)

  C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
  C_emp_missing <- mean(C_emp, na.rm = T)

#----------------------------------------------------------------
  # Concordance Index
  g11 <- dat.PRO[dat.PRO$Visit == 'Visit_4' & dat.PRO$Group == 'Group_1', 'PRO_modified']
  g22 <- dat.PRO[dat.PRO$Visit == 'Visit_4' & dat.PRO$Group == 'Group_2', 'PRO_modified']
      wt.modified <- wilcox.test(x = g11, y = g22, alternative = "two.sided")

  g11.re <- sample(x = g11, size = 1e4, replace = T)
  g22.re <- sample(x = g22, size = 1e4, replace = T)

  C_emp_modified <- ifelse(g22.re > g11.re, 1,
                             ifelse(g22.re < g11.re, 0, 0.5))
  C_emp_modified <- mean(C_emp_modified)
  # This is the empiircal C index



#----------------------------------------------------------------------
  C_true <- 1/(1 + exp(-1*sim.out$Beta['GroupGroup_2:VisitVisit_4', ]/1.52))
  #OR <- exp(sim.out$Beta['GroupGroup_2:Time_factorTime_4', ])
  #(OR^0.66)/(1 + OR^0.66)
  C_clm <- 1/(1 + exp(1*emm.clm/1.52))
  C_clm_missing <- 1/(1 + exp(1*emm.clm.missing/1.52))
  C_clm_modified <- 1/(1 + exp(1*emm.clm.modified/1.52))
  C_mixed_modified <- 1/(1 + exp(1*emm.mixed.modified/1.52))
  C_mixed_modified2 <- 1/(1 + exp(1*emm.mixed.modified2/1.52))


  pval2 <- rbind(pval2, c('Rank sum - complete' = wt.comp$p.value,
                          'Rank sum - missing' = wt.missing$p.value,
                          'Rank sum - modified' = wt.modified$p.value))

  # Output
  c.index <- rbind(c.index, c(
    'Empirical_C_index_complete' = C_emp_complete,
    'Empirical_C_index_missing' = C_emp_missing,
    'Empirical_C_index_modified' = C_emp_modified,
    'CLM_True_Beta' = C_true,
    'CLM_complete' = C_clm,
    'CLM_missing' = C_clm_missing,
    'CLM_modified' = C_clm_modified,
    'Mixed_modified' = C_mixed_modified,
    'Mixed_modified2' = C_mixed_modified2))


  cat('Replication: ', repl, '\n')

}# end replications


mod <- ordinal::clmm2(as.factor(PRO_modified) ~ Group + Visit + Group*Visit,
                      random = USUBJID, data = dat.PRO)

predict(mod.mixed.modified)
fitted(mod.mixed.modified)
#----------------------------------
# Results

# Power
apply(pval < 0.05, 2, mean)
apply(pval2 < 0.05, 2, mean)
apply(pval3 < 0.05, 2, mean) # perfect, power is 0.80 lol so perfect

# Concordance Probability:
#apply(c.index, 2, mean)
apply(1-c.index, 2, mean)[ 1:3]
apply(1-c.index, 2, mean)[ 4:ncol(c.index)]


# let's pull a replication we like:
use.repl <- pval[, 1] > 0.05 & # MMRM approach is not stat sig
  pval[,3] > 0.05 &  # prop odds model not stat sig
  pval[,4] < 0.05 &  # modified approach is stat sig
  pval2[,2] > 0.05 & # Wilcoxon test is not stat sig
  pval2[,3] < 0.05 & # rank sum test is stat sig  w modified approach
  pval3[,1] < 0.05  # Cox PH model is stat sig

which(use.repl)
