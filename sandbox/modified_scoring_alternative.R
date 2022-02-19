

# Test out the other way of adding a missing value!


rm(list = ls())
gc()

source('./sandbox/sim_surv_data.R')
source('./sandbox/sim_PRO_data.R')
library(ordinal)
library(emmeans)
library(nlme)

repl <- 1
set.seed(02042022 + repl)
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
      max.time <- dat.surv[dat.surv$USUBJID == i , 'time.surv']
      missing <- dat.PRO[dat.PRO$USUBJID == i, 'Time'] > max.time
      dat.PRO[dat.PRO$USUBJID == i, 'PRO'][which(missing)] <- NA
      event <- which.min( dat.PRO[dat.PRO$USUBJID == i, 'Time'][which(missing)])
      dat.PRO[dat.PRO$USUBJID == i, 'PRO'][event] <- 4

}

  dat.PRO[dat.PRO$USUBJID == i, ]

#-------------------
  # Modified Scoring Algorithm: make drop-out the highest/worst score:
    #unique(as.factor(dat.PRO$PRO))
    #dat.PRO$PRO_modified <- dat.PRO$PRO
    #dat.PRO$PRO_modified[is.na(dat.PRO$PRO_modified)] <- 4
#------------------------------------------------------------------------
  # CLM
  mod.clm <- ordinal::clm(as.factor(Y_comp) ~ Group + Time_factor + Group*Time_factor,
                          nAGQ = 5, data= dat.PRO)

  mod.clm.missing <- ordinal::clm(as.factor(PRO) ~ Group + Time_factor + Group*Time_factor,
                          nAGQ = 5, data= dat.PRO)

  # Marginal Means
  emm.clm <- as.data.frame(emmeans::emmeans(mod.clm, pairwise ~ Group | Time_factor)$con)
  emm.clm.missing <- as.data.frame(emmeans::emmeans(mod.clm.missing, pairwise ~ Group | Time_factor)$con)

  # Compare them
  emm.clm
  emm.clm.missing
  # okay so average 1.11 (this should be 0.75) and 1.19
  sim.out$Beta_PRO

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
  C_emp_complete

#----C_emp_complete--------------------------------------------------------
  # Concordance Index
  g1 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_1', 'PRO']
  g2 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_2', 'PRO']
      #wt.missing <- wilcox.test(x = g1, y = g2, alternative = "two.sided")

  g1.re <- sample(x = g1, size = 1e4, replace = T)
  g2.re <- sample(x = g2, size = 1e4, replace = T)

  C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
  C_emp_missing <- mean(C_emp, na.rm = T)

  C_emp_complete
  C_emp_missing

    emm.clm <- emm.clm['4', 'estimate']
  emm.clm.missing <- emm.clm.missing['4', 'estimate']
 C_clm <- 1/(1 + exp(1*emm.clm/1.52))
  C_clm_missing <- 1/(1 + exp(1*emm.clm.missing/1.52))

  C_clm
C_clm_missing
 C_emp_complete
  C_emp_missing
