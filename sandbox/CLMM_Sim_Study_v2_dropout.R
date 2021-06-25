
# 6.25.21:
# Conditional Model
# How does it handle drop-out?


rm(list = ls())
gc()

library(SPR)
library(ordinal)
library(emmeans)

# Initialize output:
score.names <- c('Y_comp', 'Y_mcar', 'Y_mar', 'Y_mnar')

out <- data.frame('Param_type' = NA,
                  'Data' = NA,
                  'Estimate' = NA,
                  'Gen_param' = NA)

number.repl <- 100


#--------------------------
  st <- Sys.time()
for (repl in 1:number.repl) {

  set.seed(6252021 + repl)

  sim.out <- SPR::sim_dat_ord_logistic_conditional(
    N = 300,
    number.groups = 2 ,
    number.timepoints = 4,
    reg.formula = formula( ~ Time + Group + Time*Group),
    Beta = 1,
    thresholds = c(-2, 0, 2),
    subject.var = 2,
    cond.mcar = F,
    Covariate = F)


  dat <- sim.out$dat

  dat <- SPR::dropout(dat = dat,
                      type_dropout  = c('mcar', 'mar', 'mnar'),
                      prop.miss = 0.3)


  aggregate(Y_mar ~ Time, FUN = function(x) sum(!is.na(x)),
            na.action = na.pass, data = dat)

  for (score in score.names) {

    # formula:
    mf <- as.formula(paste0('ordered(', score,') ~ Group + Time + Group*Time + (1|USUBJID)'))

    # CLMM - Random Effect
    mod.clmm <- ordinal::clmm(mf, nAGQ = 5, data= dat)

    tmp <- data.frame(
      'Param_type' = names(coef(mod.clmm)),
      'Data' = score,
      'Estimate' = coef(mod.clmm),
      'Gen_param' = c(sim.out$thresholds, as.vector(sim.out$Beta)[-1]))
    # drop the intercept from Beta!

    out <- rbind.data.frame(out, tmp)

    # Marginal Means
    emm.clmm <-
      suppressMessages(
        emmeans::emmeans(mod.clmm, pairwise ~ Group | Time)
      )

    emm <- as.data.frame(emm.clmm$contrasts)

    # Output
    tmp <- data.frame(
      'Param_type' = paste0('emm_', emm$Time),
      'Data' = score,
      'Estimate' = emm$estimate,
      'Gen_param' = as.vector(-1*sim.out$Beta[5:8,]))

    out <- rbind.data.frame(out, tmp)


  }# score loop

  cat(paste0('Replication: ', repl, '\n'))

}# end replications

# Runtime:
et <- Sys.time()
et - st

#---------------------------------------------------------------
# Data Mgmt
out$Data <- ifelse(out$Data == 'Y_comp', 'Complete',
       ifelse(out$Data == 'Y_mcar', 'MCAR',
              ifelse(out$Data == 'Y_mar', 'MAR',
                    ifelse(out$Data == 'Y_mnar', 'MNAR', NA))))

out$Data <- factor(out$Data, levels = c('Complete', 'MCAR', 'MAR', 'MNAR'))
#-------------------------------------------------------------
# Output
tab.out <- aggregate(
  cbind(Estimate,'Bias' = Estimate - Gen_param) ~ Data + Param_type,
  FUN = function(x)  mean(x, na.rm = T),
  na.action = na.pass,
  data = out)

tab.out
tab.out[tab.out$Param_type == 'emm_Time_4', ]
tab.out[tab.out$Param_type == 'GroupGroup_2:TimeTime_4', ]
# g2g
