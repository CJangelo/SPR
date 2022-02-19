

# 1.28.22 - looks like you can really break this thing this way
# next step is to see if you can create an ordered scale that works better!

# Once you know this works, put it all together into a cleaner function
# pass the correct variables


rm(list = ls())
gc()

source('./sandbox/sim_surv_data.R')
source('./sandbox/sim_PRO_data.R')
library(ordinal)
library(emmeans)

out <- vector()

for (repl in 1:30){

  set.seed(02032022 + repl)
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

  # Modified Scoring Algorithm: make drop-out the highest/worst score:
    unique(as.factor(dat.PRO$PRO))
    dat.PRO$PRO_modified <- dat.PRO$PRO
    dat.PRO$PRO_modified[is.na(dat.PRO$PRO_modified)] <- 5


   #emmeans fucks up the contrasts
  dat.PRO$Group <- factor(dat.PRO$Group, levels = c('Group_2', 'Group_1'))

#-------------------------------------
  # cor(dat.PRO$Y_comp, dat.PRO$PRO, use = 'pairwise')
  # aggregate(PRO ~ Group + Time_factor,
  #           function(x) mean(is.na(x)),
  #           data = dat.PRO,
  #           na.action = na.pass)
  #
  #   aggregate(cbind(Y_comp, PRO) ~ Group + Time_factor,
  #           function(x) mean(x, na.rm = T),
  #           data = dat.PRO,
  #           na.action = na.pass)
  #------------------------------------------------------------------------


  #-------------------------------------------------
  # Model Fitting
  # Model 1 - full data
  mod.clmm <- ordinal::clmm(as.factor(Y_comp) ~ Group + Time_factor +
                              Group*Time_factor + (1|USUBJID),
                            nAGQ = 5, data= dat.PRO)

  emm.clmm <- emmeans::emmeans(mod.clmm, pairwise ~ Group | Time_factor)


  # Model 2 - Missing data
  mod.clmm2 <- ordinal::clmm(as.factor(PRO) ~ Group + Time_factor +
                              Group*Time_factor + (1|USUBJID),
                            nAGQ = 5, data= dat.PRO)

  emm.clmm2 <- emmeans::emmeans(mod.clmm2, pairwise ~ Group | Time_factor)


  # Model 3 - modified scoring algorithm
  mod.clmm3 <- ordinal::clmm(as.factor(PRO_modified) ~ Group + Time_factor +
                              Group*Time_factor + (1|USUBJID),
                            nAGQ = 5, data= dat.PRO)

  emm.clmm3 <- emmeans::emmeans(mod.clmm3, pairwise ~ Group | Time_factor)


# Concordance Probability:
  # ----------------------
    # Empirical Concordance Probability
    g1 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_1', 'Y_comp']
    g2 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_2', 'Y_comp']

    g1.re <- sample(x = g1, size = 100, replace = T)
    g2.re <- sample(x = g2, size = 100, replace = T)

    C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
    C_emp_complete <- mean(C_emp)
#----------------------------------------------------------

 # ----------------------
    # Empirical Concordance Probability
    g1 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_1', 'PRO_modified']
    g2 <- dat.PRO[dat.PRO$Time_factor == 'Time_4' & dat.PRO$Group == 'Group_2', 'PRO_modified']

    g1.re <- sample(x = g1, size = 1000, replace = T)
    g2.re <- sample(x = g2, size = 1000, replace = T)

    C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
    C_emp_modified <- mean(C_emp)

#----------------------------------------------------------
  # Use generating Beta
    true.beta <- tail(sim.out$Beta, 1)
    C_true <- 1/(1 + exp(-1*true.beta/1.52))

  # Full data estimate
    beta.hat <- tail(as.data.frame(emm.clmm$con)[, 'estimate'], 1)
    C_full <- 1/(1 + exp(-1*beta.hat/1.52))

  # Missing data estimate
    beta.hat2 <- tail(as.data.frame(emm.clmm2$con)[, 'estimate'], 1)
    C_miss <- 1/(1 + exp(-1*beta.hat2/1.52))

  # Proposed Modified Scoring Algorithm
    beta.hat3 <- tail(as.data.frame(emm.clmm3$con)[, 'estimate'], 1)
    C_proposed <- 1/(1 + exp(-1*beta.hat3/1.52))


    out <- rbind(out, c(
      'True' = C_true,
      'Empirical_complete' = C_emp_complete,
      'Empirical_modified' = C_emp_modified,
      'Full_data' = C_full,
      'Missing_data' = C_miss,
      'Proposed_approach' = C_proposed))

    cat('Replication: ', repl, '\n')

}

# ----------------------
# Results
    apply(out, 2, mean)

# apply(out, 2, mean)
#               True Empirical_complete Empirical_modified          Full_data       Missing_data
#          0.6587873          0.6198333          0.6462000          0.6530969          0.6953801
#  Proposed_approach
#          0.7391053
