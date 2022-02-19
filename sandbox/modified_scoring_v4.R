

# 2.7.2022
# Modified Scoring version 4 - this assigns a score of "4" at the visit
# closest to the time of event. Visits after that timepoint are NA.
# So this is closer to a proper discrete time to event survival model.


# Does not seem to help a ton, although i guess you do get a bump up
# Use the other approach
#
#----------------------------------------------------------------------

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

# Set that as PRO missing
colnames(dat.PRO)[ colnames(dat.PRO) == 'PRO' ] <- 'PRO_missing'
#dat.PRO$PRO <- dat.PRO$Y_comp

# Now include the alternative modified scoring approach:
#------------------------------
# Add "4" at event, then missing afterwards:
# i <- 'Subject_0002'
  dat.PRO$PRO_modified <- dat.PRO$Y_comp
  id <- unique(dat.PRO$USUBJID)
for (i in id) {
      max.time <- dat.surv[dat.surv$USUBJID == i , 'time.surv']
      missing <- dat.PRO[dat.PRO$USUBJID == i, 'Time'] > max.time
      dat.PRO[dat.PRO$USUBJID == i, 'PRO_modified'][which(missing)] <- NA
      event <- min( dat.PRO[dat.PRO$USUBJID == i, 'Time'][which(missing)])
      dat.PRO[dat.PRO$USUBJID == i &
                dat.PRO$Time == event, 'PRO_modified'] <- 4

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


  # MMRM - typical approach!
  dat.mmrm <- dat.PRO
  dat.mmrm$trans <-
    100*(dat.mmrm$Y_comp)/(diff(range(dat.mmrm$Y_comp, na.rm = T)))
  tmp <- dat.mmrm[dat.mmrm$Time_factor == 'Time_1', c('USUBJID', 'trans')]
    colnames(tmp) <- c('USUBJID', 'trans_bl')
  dat.mmrm <- merge(x = dat.mmrm, y = tmp, by = 'USUBJID')
  dat.mmrm$delta <- dat.mmrm$trans - dat.mmrm$trans_bl

  #barplot(table(dat.mmrm$delta))
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

  mod.clm.missing <- ordinal::clm(as.factor(PRO_missing) ~ Group + Time_factor + Group*Time_factor,
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
  g1 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_1', 'PRO_missing']
  g2 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_2', 'PRO_missing']
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
  C_emp_modified <- mean(C_emp_modified, na.rm = T)
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
#
# > #----------------------------------
# > # Results
# >
# > # Power
# > apply(pval < 0.05, 2, mean)
#        mmrm PO_complete  PO_missing PO_modified
#        0.58        0.84        0.26        0.53
# > apply(pval2 < 0.05, 2, mean)
# Rank sum - complete  Rank sum - missing Rank sum - modified
#                0.85                0.29                0.50
# > apply(pval3 < 0.05, 2, mean) # perfect, power is 0.80 lol so perfect
# pvalue
#   0.96
# >
# > # Concordance Probability:
# > #apply(c.index, 2, mean)
# > apply(1-c.index, 2, mean)[ 1:3]
# Empirical_C_index_complete  Empirical_C_index_missing Empirical_C_index_modified
#                  0.5903540                  0.5619345                  0.5796407
# > apply(1-c.index, 2, mean)[ 4:ncol(c.index)]
# CLM_True_Beta  CLM_complete   CLM_missing  CLM_modified
#     0.6131390     0.5978384     0.5658118     0.5909688
# >
