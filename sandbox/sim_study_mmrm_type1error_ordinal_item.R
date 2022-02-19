


#------------------------------------
# 2.8.22 - I think the typical approach to evaluating a single ordinal item
# using MMRM inflates the Type I error. This would be useful to note in
# the manuscript.

# Result: no evidence of Type I error inflation
             #  CLM              CLMM      MMRM_glmmTMB MMRM_lmer_Kenward
             # 0.03              0.03              0.07              0.06

rm(list = ls())
gc()

library(SPR)
library(ordinal)
#library(multgee)

# Initialize output:
#score.names <- c('Y_comp', 'Y_mcar', 'Y_mar', 'Y_mnar')
score.names <- c('Y_comp')
pval <- vector()
number.repl <- 100
repl <- 1


#--------------------------
  st <- Sys.time()
for (repl in 1:number.repl) {

  set.seed(02082022 + repl)


  sim.out <- SPR::sim_dat_ord_logistic_conditional(N = 200 ,
                                             number.groups = 2,
                                             number.timepoints = 4,
                                             reg.formula = formula( ~ Time + Group + Time*Group),
                                             Beta = 0,
                                             thresholds = c( -1, 0, 1),
                                             subject.var = 1,
                                             cond.mcar = F,
                                             Covariate = F)

  dat <- sim.out$dat

    # formula:
    mf <- as.formula(ordered(Y_comp) ~ Group + Time + Group*Time)

    # CLM - no random effect
    mod.ols1 <- ordinal::clm(mf, nAGQ = 5, data= dat)
    mod.emms.ols1 <- emmeans(mod.ols1, pairwise ~ Group | Time,
                             adjust = 'none')


    # # Generalized Estimating Equations:
    # dat.gee <- dat
    #   tmp <- order(dat.gee$id.geepack)
    # dat.gee <- dat.gee[tmp, ]
    # mod.gee <- ordLORgee(formula = mf,
    #                  data = dat.gee,
    #                  id = id.geepack,
    #                  repeated = Time, LORstr = "uniform")

  # Current Practice: MMRM
  dat.mmrm <- dat
  dat.mmrm$trans <-
    100*(dat.mmrm$Y_comp)/(diff(range(dat.mmrm$Y_comp, na.rm = T)))
  tmp <- dat.mmrm[dat.mmrm$Time == 'Time_1', c('USUBJID', 'trans')]
    colnames(tmp) <- c('USUBJID', 'trans_bl')
  dat.mmrm <- merge(x = dat.mmrm, y = tmp, by = 'USUBJID')
  dat.mmrm$delta <- dat.mmrm$trans - dat.mmrm$trans_bl


 # MMRM -
  mod.lmer1 <- lme4::lmer(delta~ Group + Time + Group*Time + ( -1 + Time | USUBJID),
                  data=dat.mmrm[dat.mmrm$Time != 'Time_1', ],
                    control = lmerControl(check.nobs.vs.nRE = 'ignore'))
  #fixef(mod.lmer1)
  mod.emms.lmer1.kr <- emmeans(mod.lmer1, pairwise ~ Group | Time, adjust = 'none', mode = "kenward" )

  #emmeans::emmeans(mod.lmer1, pairwise ~ Group|Time, adjust = 'none',  mode = 'df.error')

  mod.mmrm <- glmmTMB::glmmTMB(delta ~ Group + Time + Group*Time + us(Time + 0 | USUBJID),
                  data=dat.mmrm[dat.mmrm$Time != 'Time_1', ],
                  REML = T,
                  dispformula=~0)
  #glmmTMB::fixef(mod.mmrm)

  emm.mmrm <- emmeans::emmeans(mod.mmrm, pairwise ~ Group|Time, adjust = 'none',  mode = 'df.error')

    # CLMM - Random Effect - True model:
  mod.clmm <- ordinal::clmm(as.factor(Y_comp) ~ Group + Time + Group*Time + (1|USUBJID),
                            nAGQ = 5, data= dat)

  emm.clmm <- emmeans::emmeans(mod.clmm, pairwise ~ Group | Time)



  pval <- rbind(pval, c(
    'CLM' = as.data.frame(mod.emms.ols1$contrasts)['4', 'p.value'],
    'CLMM' = as.data.frame(emm.clmm$contrasts)['4', 'p.value'],
    'MMRM_glmmTMB' = as.data.frame(emm.mmrm$contrasts)['3', 'p.value'],
    'MMRM_lmer_Kenward' = as.data.frame(mod.emms.lmer1.kr$contrasts)['3', 'p.value']))



  cat('Replication: ', repl, '\n')

}# end replications


#---------------------------------------------
# Results
apply(pval < 0.05, 2, mean)
# > apply(pval < 0.05, 2, mean)
#               CLM              CLMM      MMRM_glmmTMB MMRM_lmer_Kenward
#              0.03              0.03              0.07              0.06
