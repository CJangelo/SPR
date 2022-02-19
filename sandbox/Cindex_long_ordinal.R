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

  sim.out <- SPR::sim_dat_ord_logistic_conditional(N = 500,
                                                   number.groups = 2 ,
                                                   number.timepoints = 3,
                                                   reg.formula = formula( ~ Time + Group + Time*Group),
                                                   Beta = 1,
                                                   thresholds = c(-2, 0, 2),
                                                   subject.var = 1,
                                                   cond.mcar = F,
                                                   Covariate = F)


  dat <- sim.out$dat


#------------
  # CLM
  mod.clm <- ordinal::clm(as.factor(Y_comp) ~ Group + Time + Group*Time,
                          nAGQ = 5, data= dat)


  # CLMM - Random Effect
  mod.clmm <- ordinal::clmm(as.factor(Y_comp) ~ Group + Time + Group*Time + (1|USUBJID),
                            nAGQ = 5, data= dat)


  # Generalized Estimating Equations:
  dat.gee <- dat
  tmp <- order(dat.gee$id.geepack)
  dat.gee <- dat.gee[tmp, ]
  mod.gee <- multgee::ordLORgee(formula = Y_comp ~ Time + Group + Time*Group,
                                data = dat.gee,
                                id = id.geepack,
                                repeated = Time, LORstr = "uniform")


  # Marginal Means
  #emmeans::emmeans(mod.gee, pairwise ~ Group | Time) NOT AVAILABLE
  emm.clm <- emmeans::emmeans(mod.clm, pairwise ~ Group | Time)
  emm.clmm <- emmeans::emmeans(mod.clmm, pairwise ~ Group | Time)
  emm.gee <- mod.gee$coeff['TimeTime_3'] +
    mod.gee$coeff['GroupGroup_2'] +
    mod.gee$coeff['TimeTime_3:GroupGroup_2']

  emm.clm <- as.data.frame(emm.clm$contrasts)['3', 'estimate']
  emm.clmm <- as.data.frame(emm.clmm$contrasts)['3', 'estimate']

  # Concordance Index
  g1 <- dat[dat$Time == 'Time_3' & dat$Group == 'Group_1', 'Y_comp']
  g2 <- dat[dat$Time == 'Time_3' & dat$Group == 'Group_2', 'Y_comp']
  g1.re <- sample(x = g1, size = 1000, replace = T)
  g2.re <- sample(x = g2, size = 1000, replace = T)

  C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
  C_emp <- mean(C_emp)
  # This is the empiircal C index

  C_true <- 1/(1 + exp(-1*sim.out$Beta['TimeTime_3:GroupGroup_2', ]/1.52))
  C_gee <- 1/(1 + exp(1*emm.gee/1.52))
  C_clm <- 1/(1 + exp(1*emm.clm/1.52))
  C_clmm <- 1/(1 + exp(1*emm.clmm/1.52))

  # Output
  c.index <- rbind(c.index, c(
    'True' = C_true,
    'Empirical' = C_emp,
    'GEE' = C_gee,
    'CLM' = C_clm,
    'mixed_mod' = C_clmm))


  est <- rbind(est, c(
               'True' = sim.out$Beta['TimeTime_3:GroupGroup_2', ],
               'GEE' = emm.gee,
               'clm' = emm.clm,
               'mixed_mod' = emm.clmm))


  cat('Replication: ', repl, '\n')

}# end replications



apply(c.index, 2, mean)
# apply(c.index, 2, mean)
#           True      Empirical GEE.TimeTime_3            CLM      mixed_mod
#      0.6587873      0.6293500      0.6353439      0.6367801      0.6588329
# >
