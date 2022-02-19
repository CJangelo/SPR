
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

for (repl in 1:100){

  set.seed(02072022 + repl)
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
  # aggregate(PRO ~ Group + Time,
  #           function(x) mean(is.na(x)),
  #           data = dat.PRO,
  #           na.action = na.pass)

  # MMRM - typical approach!
  dat.mmrm <- dat.PRO
  dat.mmrm$trans <-
    100*(dat.mmrm$PRO)/(diff(range(dat.mmrm$PRO, na.rm = T)))
  tmp <- dat.mmrm[dat.mmrm$Time_factor == 'Time_1', c('USUBJID', 'trans')]
    colnames(tmp) <- c('USUBJID', 'trans_bl')
  dat.mmrm <- merge(x = dat.mmrm, y = tmp, by = 'USUBJID')
  dat.mmrm$delta <- dat.mmrm$trans - dat.mmrm$trans_bl

  # barplot(table(dat.mmrm$delta))
  # mod.ols <- lm(delta ~ Group + Time_factor + Group*Time_factor,
  #               data = dat.mmrm[dat.mmrm$Time_factor != 'Time_1', ])
  # emmeans(mod.ols, pairwise ~ Group|Time_factor)
  mod.mmrm <- nlme::gls(delta ~ Group + Time_factor + Group*Time_factor,
                        data = dat.mmrm[dat.mmrm$Time_factor != 'Time_1', ],
                        correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                        weights = varIdent(form = ~ 1 | Time_factor),          #  freely estimate variance at subsequent timepoints
                        na.action = na.exclude)
  emm.mmrm <- emmeans(mod.mmrm, pairwise ~ Group|Time_factor, adjust = 'none',  mode = 'df.error')

  #
  # aggregate(cbind(Y_comp, PRO, PRO_modified) ~ Group + Time,
  #           function(x) mean(x, na.rm = T),
  #           data = dat.PRO,
  #           na.action = na.pass)
  # Initial look is that it over-corrects, but we are more interested
  # in the Concordance Probability, which comes from the Beta estimate

#------------------------------------------------------------------------

#------------
  # CLM
  mod.clm <- ordinal::clm(as.factor(Y_comp) ~ Group + Time_factor + Group*Time_factor,
                          nAGQ = 5, data= dat.PRO)

  mod.clm.missing <- ordinal::clm(as.factor(PRO) ~ Group + Time_factor + Group*Time_factor,
                          nAGQ = 5, data= dat.PRO)

  mod.clm.modified <- ordinal::clm(as.factor(PRO_modified) ~ Group + Time_factor + Group*Time_factor,
                          nAGQ = 5, data= dat.PRO)

  # Marginal Means
  emm.clm <- emmeans::emmeans(mod.clm, pairwise ~ Group | Time_factor)
  emm.clm.missing <- emmeans::emmeans(mod.clm.missing, pairwise ~ Group | Time_factor)
  emm.clm.modified <- emmeans::emmeans(mod.clm.modified, pairwise ~ Group | Time_factor)

  pval <- rbind(pval, c(
    'mmrm' = as.data.frame(emm.mmrm$contrasts)['3', 'p.value'], # Timepoint 3 bc no baseline!
    'PO_complete' = as.data.frame(emm.clm$contrasts)['4', 'p.value'],
    'PO_missing' = as.data.frame(emm.clm.missing$contrasts)['4', 'p.value'],
    'PO_modified' = as.data.frame(emm.clm.modified$contrasts)['4', 'p.value']))
 #
 #
  emm.clm <- as.data.frame(emm.clm$contrasts)['4', 'estimate']
  emm.clm.missing <- as.data.frame(emm.clm.missing$contrasts)['4', 'estimate']
  emm.clm.modified <- as.data.frame(emm.clm.modified$contrasts)['4', 'estimate']

#------------------------------------------------------------
  # Concordance Index
  g1 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_1', 'Y_comp']
  g2 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_2', 'Y_comp']
    wt.comp <- wilcox.test(x = g1, y = g2, alternative = "two.sided")
  g1.re <- sample(x = g1, size = 1e4, replace = T)
  g2.re <- sample(x = g2, size = 1e4, replace = T)

  C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
  C_emp_complete <- mean(C_emp, na.rm = T)

#------------------------------------------------------------
  # Concordance Index
  g1 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_1', 'PRO']
  g2 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_2', 'PRO']
      wt.missing <- wilcox.test(x = g1, y = g2, alternative = "two.sided")

  g1.re <- sample(x = g1, size = 1e4, replace = T)
  g2.re <- sample(x = g2, size = 1e4, replace = T)

  C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
  C_emp_missing <- mean(C_emp, na.rm = T)

#----------------------------------------------------------------
  # Concordance Index
  g11 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_1', 'PRO_modified']
  g22 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_2', 'PRO_modified']
      wt.modified <- wilcox.test(x = g11, y = g22, alternative = "two.sided")

  g11.re <- sample(x = g11, size = 1e4, replace = T)
  g22.re <- sample(x = g22, size = 1e4, replace = T)

  C_emp_modified <- ifelse(g22.re > g11.re, 1,
                             ifelse(g22.re < g11.re, 0, 0.5))
  C_emp_modified <- mean(C_emp_modified)
  # This is the empiircal C index

  C_true <- 1/(1 + exp(-1*sim.out$Beta['GroupGroup_2:Time_factorTime_4', ]/1.52))
  #OR <- exp(sim.out$Beta['GroupGroup_2:Time_factorTime_4', ])
  #(OR^0.66)/(1 + OR^0.66)
  C_clm <- 1/(1 + exp(1*emm.clm/1.52))
  C_clm_missing <- 1/(1 + exp(1*emm.clm.missing/1.52))
  C_clm_modified <- 1/(1 + exp(1*emm.clm.modified/1.52))


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
    'CLM_modified' = C_clm_modified))


  cat('Replication: ', repl, '\n')

}# end replications


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

# #
# > apply(1-c.index, 2, mean)[ 1:3]
# Empirical_C_index_complete  Empirical_C_index_missing Empirical_C_index_modified
#                  0.6411083                  0.6092965                  0.6565900
# > apply(1-c.index, 2, mean)[ 4:ncol(c.index)]
# CLM_True_Beta  CLM_complete   CLM_missing  CLM_modified
#     0.6587873     0.6441456     0.6082749     0.7082663

# let's pull a replication we like:
use.repl <- pval[, 1] > 0.05 & # MMRM approach is not stat sig
  pval[,3] > 0.05 &  # prop odds model not stat sig
  pval[,4] < 0.05 &  # modified approach is stat sig
  pval2[,2] > 0.05 & # Wilcoxon test is not stat sig
  pval2[,3] < 0.05 & # rank sum test is stat sig  w modified approach
  pval3[,1] < 0.05  # Cox PH model is stat sig

# Pull a replication (one of the 33), and plot the KM curves, add the p-value
# Plot the MMRM change from baseline curves, add not stat sig!
# Table:
# PH Cox reg model p-value
# MMRM p-value
# Wilcoxon p-value
# Prop odds model
# Modified approach
# Wilcoxon
# prop odds
# Compute the concordance probability,
# show that PO model yields easily interpreted output, same as Wilcoxon test

# Expected outcomes:
# Show that the PO model can also acommodate demo and clinical characteristics
# Can yield predictions about expected outcomes

# Then in Appendix or somewhere, show power study
# Basically all the same thing, with and without
# power much higher, seems to make sense given what we see here
# not clear that power even makes sense here, is there a true difference?
# I guess between the two models, both have a stat sig txa - both the PRO and PH
# a little strange though. Sort of a JM model?
table(use.repl)
# 38% meet this criteria, you're golden

#
# > # Power
# > apply(pval < 0.05, 2, mean)
#        mmrm PO_complete  PO_missing PO_modified
#        0.35        0.82        0.30        0.96
# > apply(pval2 < 0.05, 2, mean)
# Rank sum - complete  Rank sum - missing Rank sum - modified
#                0.81                0.31                0.94
# > apply(pval3 < 0.05, 2, mean) # perfect, power is 0.80 lol so perfect
# pvalue
#   0.93
#
# # Concordance Probability:
# > #apply(c.index, 2, mean)
# > apply(1-c.index, 2, mean)[ 1:3]
# Empirical_C_index_complete  Empirical_C_index_missing Empirical_C_index_modified
#                  0.5890300                  0.5632694                  0.6167785
# > apply(1-c.index, 2, mean)[ 4:ncol(c.index)]
# CLM_True_Beta  CLM_complete   CLM_missing  CLM_modified
#     0.6131390     0.5956188     0.5675418     0.6600486
# >
