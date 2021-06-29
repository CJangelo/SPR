
# Mixed Effect Model recovering MMRM parameter estimates
rm(list = ls())
gc()

library(SPR)
library(MASS)
library(emmeans)
library(nlme)
library(lme4)

# Initialize output:
score.names <- c('Y_comp', 'Y_mcar', 'Y_mar', 'Y_mnar')

out <- data.frame('Approach' = NA,
                  'Param_type' = NA,
                  'Data' = NA,
                  'Estimate' = NA,
                  'Gen_param' = NA)

number.repl <- 1000


#--------------------------
st <- Sys.time()
for (repl in 1:number.repl) {

  set.seed(6302021 + repl)

  sim.out <- SPR::sim_dat(N = 300,
                          number.groups = 2 ,
                          number.timepoints = 4,
                          reg.formula = formula( ~ Time + Group + Time*Group),
                          Beta = 1,
                          corr = 'ar1',
                          cor.value = 0.8,
                          var.values = 1)



  dat <- sim.out$dat

  dat <- SPR::dropout(dat = dat,
                      type_dropout  = c('mcar', 'mar', 'mnar'),
                      prop.miss = 0.3)


  for (score in score.names) {

    # formula:
    f1 <- as.formula(paste0(score, '~ Group + Time + Group*Time'))

    # MMRM
    mod.gls1 <- gls( f1,
                     data = dat,
                     correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                     weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                     na.action = na.exclude)


    f2 <- as.formula(paste0(score, '~ Group + Time + Group*Time + (1|USUBJID)'))

    # Mixed Effects
    mod.lmer1 <- lmer(f2, data = dat)

    # Compute marginal means
    #library(emmeans)
    emms.gls1 <- emmeans(mod.gls1,
                         pairwise ~ Group | Time,
                         adjust = 'none',
                         mode = 'df.error')

    emms.lmer1 <- emmeans(mod.lmer1,
                          pairwise ~ Group | Time,
                          adjust = 'none',
                          mode = 'asymptotic')

    emms.gls1 <- as.data.frame(emms.gls1$contrasts)
    emms.lmer1 <- as.data.frame(emms.lmer1$contrasts)


    # Output:
    tmp <- data.frame(
      'Approach' = 'gls',
      'Param_type' = emms.gls1$Time,
      'Data' = score,
      'Estimate' = emms.gls1$estimate,
      'Gen_param' = c(0, 0.25, 0.625, 1) )

    out <- rbind.data.frame(out, tmp)

    # Output
    tmp <- data.frame(
      'Approach' = 'lmer',
      'Param_type' = emms.lmer1$Time,
      'Data' = score,
      'Estimate' = emms.lmer1$estimate,
      'Gen_param' = c(0, 0.25, 0.625, 1) )

    out <- rbind.data.frame(out, tmp)

    # Reset, avoid accidental duplicates
    emms.lmer1 <- emms.gls1 <- NA


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
  cbind(Estimate,'Bias' = Estimate - Gen_param) ~ Approach + Data + Param_type,
  FUN = function(x)  mean(x, na.rm = T),
  na.action = na.pass,
  data = out)

tab.out
tab.out[tab.out$Param_type == 'Time_4', ]
# ???
