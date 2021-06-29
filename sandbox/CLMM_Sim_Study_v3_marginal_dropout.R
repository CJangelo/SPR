
# Marginal Model
# How does it handle drop-out?

rm(list = ls())
gc()

library(SPR)
library(ordinal)
library(multgee)
library(polycor)

# Initialize output:
score.names <- c('Y_comp', 'Y_mcar', 'Y_mar', 'Y_mnar')

out <- data.frame('Approach' = NA,
                  'Param_type' = NA,
                  'Data' = NA,
                  'Estimate' = NA,
                  'Gen_param' = NA)

number.repl <- 100


#--------------------------
  st <- Sys.time()
for (repl in 1:number.repl) {

  set.seed(6292021 + repl)

  sim.out <- SPR::sim_dat_ord_logistic(
    N = 300,
    number.groups = 2 ,
    number.timepoints = 4,
    reg.formula = formula( ~ Time + Group + Time*Group),
    Beta = 1,
    thresholds = c(-1, 0, 1),
    corr = 'ar1',
    cor.value = 0.8)


  dat <- sim.out$dat

  dat <- SPR::dropout(dat = dat,
                      type_dropout  = c('mcar', 'mar', 'mnar'),
                      prop.miss = 0.6)

# # Check your data:
#   aggregate(Y_mnar ~ Time, FUN = function(x) sum(!is.na(x)),
#             na.action = na.pass, data = dat)
#
#   aggregate(cbind(Y_comp, Y_mar, Y_mnar) ~ Time + Group, FUN = function(x) mean(x, na.rm = T),
#             na.action = na.pass, data = dat)
#
# # Check the correlations:
# tmp <- tidyr::pivot_wider(data = dat,
#                    id_cols = c('USUBJID'),
#                    names_from = c('Time'),
#                    values_from = c('Y_comp'))
# polycor::hetcor(tmp[, -1])$correlations



  for (score in score.names) {

    # formula:
    mf <- as.formula(paste0('ordered(', score,') ~ Group + Time + Group*Time'))

    # CLM - no random effect
    mod.clm <- ordinal::clm(mf, nAGQ = 5, data= dat)

    # Generalized Estimating Equations:
    dat.gee <- dat
      tmp <- order(dat.gee$id.geepack)
    dat.gee <- dat.gee[tmp, ]
    mod.gee <- ordLORgee(formula = mf,
                     data = dat.gee,
                     id = id.geepack,
                     repeated = Time, LORstr = "uniform")


  # Output:
    tmp <- data.frame(
      'Approach' = 'clm',
      'Param_type' = names(coef(mod.clm)),
      'Data' = score,
      'Estimate' = coef(mod.clm),
      'Gen_param' = c(sim.out$thresholds, as.vector(sim.out$Beta)[-1]))
    # drop the intercept from Beta!

    out <- rbind.data.frame(out, tmp)

    # Output
    tmp <- data.frame(
      'Approach' = 'gee',
      'Param_type' = names(coef(mod.gee)),
      'Data' = score,
      'Estimate' = coef(mod.gee),
      'Gen_param' = c(sim.out$thresholds, as.vector(sim.out$Beta)[-1]))

    out <- rbind.data.frame(out, tmp)

    # Test:
    mod.clm <- mod.gee <- NA

      # Conditional Model fit to Marginal Data is SLOW!
    # skip for now
    # mff <- as.formula(paste0('ordered(', score,') ~ Group + Time + Group*Time + (1|USUBJID)'))
    # mod.clmm <- ordinal::clmm(mff, nAGQ = 5, data= dat)
    # sigma2 <- ordinal::VarCorr(mod.clmm)
    #   sigma2 <- as.numeric(as.data.frame(sigma2))
    #   denom.logit <- sqrt( (sigma2 + (pi^2/3))/(pi^2/3) )


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
tab.out[tab.out$Param_type == 'GroupGroup_2:TimeTime_4', ]
# ???

sum(coef(mod.gee)[c('GroupGroup_2', 'GroupGroup_2:TimeTime_4')])
    sum(coef(mod.clm)[c('GroupGroup_2', 'GroupGroup_2:TimeTime_4')])
