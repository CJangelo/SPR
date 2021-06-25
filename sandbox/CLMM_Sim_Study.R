
# 6.25.21:
# it looks like the conditional data generation does
# NOT allow you to just fit a model without random intercepts
# this is distinct from the marginal models!

# marginal models only the SE is affected, i think
# with the conditional models, you can't ignore the subject effects
# or else you bias the estimates!

# So of course the fixed effects only model (clm) won't handle
# missing data correctly, because it can't handle complete data correctly!

#

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


#--------------------------
for (repl in 1:number.repl) {

  set.seed(6252021 + repl)

  sim.out <- SPR::sim_dat_ord_logistic_conditional(N = 300,
                                                   number.groups = 2 ,
                                                   number.timepoints = 4,
                                                   reg.formula = formula( ~ Time + Group + Time*Group),
                                                   Beta = 1,
                                                   thresholds = c(-2, 0, 2),
                                                   subject.var = 2,
                                                   cond.mcar = F,
                                                   Covariate = F)


  dat <- sim.out$dat



# Visualization:
  # Quick barplot, base R plot functions
  # barplot(100*table(dat$Y_comp)/sum(table(dat$Y_comp)),
  #       ylim = c(0, 100), ylab = 'Percentage',
  #       col = 'grey', main = 'Ordinal')



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
  # emmeans::emmeans(mod.gee, pairwise ~ Group | Time) NOT AVAILABLE
  emm.clm <- emmeans::emmeans(mod.clm, pairwise ~ Group | Time)
  emm.clmm <- emmeans::emmeans(mod.clmm, pairwise ~ Group | Time)


  # Output
  out.gee <- rbind(out.gee, coef(mod.gee))
  out.clm <- rbind(out.clm, coef(mod.clm))
  out.clmm <- rbind(out.clmm, coef(mod.clmm))
  out.emms.clm  <- rbind(out.emms.clm,
                         as.data.frame(emm.clm$contrast)$estimate)
  out.emms.clmm <- rbind(out.emms.clmm,
                         as.data.frame(emm.clmm$contrast)$estimate)



  cat(paste0('Replication: ', repl, '\n'))

}# end replications

sim.out$Beta
sim.out$thresholds
colMeans(out.clm)
colMeans(out.clmm)
colMeans(out.gee)
colMeans(out.emms.clm)
colMeans(out.emms.clmm)

# All parameters
cbind('clm' = colMeans(out.clm),
      'gee' = colMeans(out.gee),
      'clmm' = colMeans(out.clmm))

# Marginal means:
cbind('clm' = colMeans(out.emms.clm),
      'clmm' = colMeans(out.emms.clmm))
