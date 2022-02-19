

# Cindex_long_ordinal
# now with drop-out!

# Marginal logistic simulation - no subject-level random effect
# just use the OLS estimation
# next step is to use GEE estimation

# I understand that this is MNAR drop-out rather than a survival process
# But I am concerned about the discrepancy between the
# empirical C-index using the modified scoring approach
# versus
# C-index estimated using the beta hat from CLM using modified scoring approach
# There is a larger difference, 6% versus the usual 2%

rm(list = ls())
gc()

library(SPR)
library(ordinal)
library(multgee)
library(emmeans)

# Initialize output:
out.clm <- vector()
out.clmm <- vector()
out.gee <- vector()
out.emms.clm <- vector()
out.emms.clmm <- vector()
number.repl <- 100
repl <- 1
est <- vector()
c.index <- vector()

#--------------------------
for (repl in 1:number.repl) {

  set.seed(02032022 + repl)

  sim.out <- SPR::sim_dat_ord_logistic(
    N = 2000,
    number.groups = 2 ,
    number.timepoints = 4,
    reg.formula = formula( ~ Time + Group + Time*Group),
    Beta = 2,
    thresholds = c(-2, 0, 2, 3),
    corr = 'ar1',
    cor.value = 0.4)

  dat <- sim.out$dat

  dat <- SPR::dropout(dat = dat,
                      type_dropout  = c('mcar', 'mar', 'mnar'),
                      prop.miss = 0.5)

#------------------------------------------------------------
  # Modified Scoring Algorithm: make drop-out the highest/worst score:
    unique(as.factor(dat$Y_comp))
    dat$Y_modified <- dat$Y_mnar
    dat$Y_modified[is.na(dat$Y_modified)] <- 5

#-------------------------------------
  # aggregate(Y_mnar ~ Group + Time,
  #           function(x) mean(is.na(x)),
  #           data = dat,
  #           na.action = na.pass)
  #
  # aggregate(cbind(Y_comp, Y_mnar, Y_modified) ~ Group + Time,
  #           function(x) mean(x, na.rm = T),
  #           data = dat,
  #           na.action = na.pass)
  # Initial look is that it over-corrects, but we are more interested
  # in the Concordance Probability, which comes from the Beta estimate

#------------------------------------------------------------------------


#------------
  # CLM
  mod.clm <- ordinal::clm(as.factor(Y_comp) ~ Group + Time + Group*Time,
                          nAGQ = 5, data= dat)

  mod.clm.mnar <- ordinal::clm(as.factor(Y_mnar) ~ Group + Time + Group*Time,
                          nAGQ = 5, data= dat)

  mod.clm.modified <- ordinal::clm(as.factor(Y_modified) ~ Group + Time + Group*Time,
                          nAGQ = 5, data= dat)


  # Marginal Means
  emm.clm <- emmeans::emmeans(mod.clm, pairwise ~ Group | Time)
  emm.clm.mnar <- emmeans::emmeans(mod.clm.mnar, pairwise ~ Group | Time)
  emm.clm.modified <- emmeans::emmeans(mod.clm.modified, pairwise ~ Group | Time)

 #
  emm.clm <- as.data.frame(emm.clm$contrasts)['4', 'estimate']
  emm.clm.mnar <- as.data.frame(emm.clm.mnar$contrasts)['4', 'estimate']
  emm.clm.modified <- as.data.frame(emm.clm.modified$contrasts)['4', 'estimate']

#------------------------------------------------------------
  # Concordance Index
  g1 <- dat[dat$Time == 'Time_4' & dat$Group == 'Group_1', 'Y_comp']
  g2 <- dat[dat$Time == 'Time_4' & dat$Group == 'Group_2', 'Y_comp']
  g1.re <- sample(x = g1, size = 1000, replace = T)
  g2.re <- sample(x = g2, size = 1000, replace = T)

  C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
  C_emp_Y_comp <- mean(C_emp)

#----------------------------------------------------------------
  # Concordance Index
  g11 <- dat[dat$Time == 'Time_4' & dat$Group == 'Group_1', 'Y_modified']
  g22 <- dat[dat$Time == 'Time_4' & dat$Group == 'Group_2', 'Y_modified']
  g11.re <- sample(x = g11, size = 1000, replace = T)
  g22.re <- sample(x = g22, size = 1000, replace = T)

  C_emp_Y_modified <- ifelse(g22.re > g11.re, 1,
                             ifelse(g22.re < g11.re, 0, 0.5))
  C_emp_Y_modified <- mean(C_emp_Y_modified)
  # This is the empiircal C index

  C_true <- 1/(1 + exp(-1*sim.out$Beta['TimeTime_4:GroupGroup_2', ]/1.52))
  C_clm <- 1/(1 + exp(1*emm.clm/1.52))
  C_clm.mnar <- 1/(1 + exp(1*emm.clm.mnar/1.52))
  C_clm.modified <- 1/(1 + exp(1*emm.clm.modified/1.52))

  # Output
  c.index <- rbind(c.index, c(
    'True' = C_true,
    'Empirical_Y_complete' = C_emp_Y_comp,
    'Empirical_Y_modified' = C_emp_Y_modified,
    'CLM_complete' = C_clm,
    'CLM_mnar' = C_clm.mnar,
    'CLM_modified' = C_clm.modified))


  cat('Replication: ', repl, '\n')

}# end replications


#----------------------------------
# Results
apply(c.index, 2, mean)

# > apply(c.index, 2, mean)
#                 True Empirical_Y_complete Empirical_Y_modified         CLM_complete
#            0.7884803            0.7720200            0.7044400            0.7881752
#             CLM_mnar         CLM_modified
#            0.6855609            0.7619360
