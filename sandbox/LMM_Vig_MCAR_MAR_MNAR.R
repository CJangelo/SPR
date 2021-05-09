
# 4.7.21 - somewhere in here you've got to show the alignment between 
# means, ANOVA, t-test, linear regression, and the MMRM when you've got complete data. FDA reviewers are idiots. 

## 3.21.21 - Vignette 1
# Go through a quick blurb on MCAR, MAR, and MNAR

# Quick example

# Vignette: MCAR, MAR and MNAR
# Examples:
# 1. MCAR (standard interpretation of MCAR - descriptive statistics sufficient!)
# 2. MAR with correlation = 0.8  (standard interpretation of MAR dropout - descriptive statistics, OLS not sufficient!)
# 3. MNAR with correlation = 0.4 (standard interprtation of MNAR - nothing works!)

# However, thinking about it more:
# 4. MAR with correlation = 0.1  (Approaches MCAR!)
# 5. MNAR with correlation = 0.95  (approaches MAR!)
# 6. MNAR with correlation = 0.0 (you're fucked!)

# Take-away: MNAR is ugly, but if correlations are high, you're good!

# Run everythign through, save the output, then do the knitr, don't want to re-run it every damn time!

# Next: 
# Type I error 
# Use FWE correction as well
# Power
# Then go through the actual training?? Seems like a lot of work
# Idea would be to show all the different types of drop-out
# For example, show that as correlation increases, MNAR becomes MAR
# can show with just one replication
# Show the relationship between MCAR, MAR, and MNAR is adjusted via those correlations
# Where is the code to use previous timepoint to predict future drop-out? 
# MCAR test, and also there was a conditional MCAR test version of that
# Where is that?
# Make Vig 2 a Type I error simulation
# Make Vig 3 a Power simulation
# bootstrap example as well!
# Then get into the illustrations of various concepts 
# Not really suited to Vignettes, where can I post that?


#---------------------------------------------------------------------------------------
#
#
#
#
#
#   VIGNETTE: HOW DO I HANDLE MISSING DATA?

rm(list = ls())
gc()

library(MASS)
# library(glmmTMB) # https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
library(emmeans)
library(nlme)
# library(lme4)
# library(lmerTest)

###
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/dropout.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/VarCov_gls.R")


  set.seed(03232021)

#----------------------------------------------------------------------------------------
# MCAR Example
# 1. MCAR (standard interpretation of MCAR - descriptive statistics sufficient!)
  
  out <- sim_dat(N = 2000, 
                 number.groups = 2 , 
                 number.timepoints = 4, 
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 2,
                 corr = 'ar1', 
                 cor.value = 0.8, 
                 var.values = 2)
                
  
  dat <- out$dat
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mcar'), 
                 prop.miss = 0.3)
  
 
#  OLS Regression Model 
  mod.ols1 <- lm(Y_comp ~ Time + Group + Time*Group, 
                 data = dat)

  mod.ols2 <- lm(Y_mcar ~ Group + Time + Group*Time, 
                            data = dat)

#   MMRM 
  library(nlme)
  mod.gls1 <- gls(Y_comp ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)


  mod.gls2 <- gls(Y_mcar ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)
  
  
  # Compute marginal mean contrasts
  library(emmeans)
  mod.emms.ols1 <- emmeans(mod.ols1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

 # mod.emms.gls1 <- emmeans(mod.gls1, pairwise ~ Group | Time, adjust = 'none',  mode = 'satterthwaite')
  mod.emms.gls1 <- emmeans(mod.gls1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
 # mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'satterthwaite')
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
# Analytical Satterthwaite method not available; using appx-satterthwaite
# Error: Can't estimate Satterthwaite parameters.
#   Try adding the argument 'mode = "df.error"'



  # Examine output:
  out$Beta
  # Look at drop-out:
  aggregate(Y_mcar ~ Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  aggregate(Y_mcar ~ Group + Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  # Approximately equal across both treatment arms - already know there's no way for the drop-out to bias
  # the treatment arm comparisons!
  
  # Descriptive Statistics 
  # Sufficient with MCAR drop out - both complete and missing align, yield unbiased estimates:
  aggregate(Y_comp ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),5), dat = dat, na.action = na.pass)
  aggregate(Y_mcar ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),5), dat = dat, na.action = na.pass)

  # OLS likewise recovers unbiased parameters using complete data AND with MCAR drop-out
  mod.emms.ols1$contrasts
  mod.emms.ols2$contrasts
  # Here we are only running a single replication to illustrate the unbiased estimates
  # if we ran a full simulation study, we would see that RMSE increases under MCAR drop-out,
  # resulting from the reduction in sample size - no getting away from that!
  mod.emms.gls1$contrasts
  mod.emms.gls2$contrasts
  # Can see that the correctly specified model - MMRM - yields the same estimates
#
# Note: Can test for MCAR drop out!
# Test for MCAR
# Does your score at a previous timepoint predict your score at the next timepoint?
# Let's make our life easy and just do the following:

dat1 <- dat[dat$Time == 'Time_3', 'Y_mcar']
dat2 <- dat[dat$Time == 'Time_4', 'Y_mcar']
dropped.out <- !is.na(dat1) & is.na(dat2) # these subjects dropped out at timepoint 4
table(dropped.out)

# Does the score at timepoint 3 predict drop-out at timepoint 4?
mod <- glm(dropped.out ~ dat1, family = 'binomial')
summary(mod)
# This is evidence suggesting that the dropout may be MCAR
# Later we will re-compute this test on MAR dropout

# ----------------------------------------------------------------------------------------------------
# MAR example
# 2. MAR with correlation = 0.8  (standard interpretation of MAR dropout - descriptive statistics, OLS not sufficient!)

  dat <- out$dat
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mar'), 
                 prop.miss = 0.3)
 
# Test for MCAR
# Does your score at a previous timepoint predict your score at the next timepoint?
dat1 <- dat[dat$Time == 'Time_3', 'Y_mar']
dat2 <- dat[dat$Time == 'Time_4', 'Y_mar']
dropped.out <- !is.na(dat1) & is.na(dat2) # these subjects dropped out at timepoint 4
table(dropped.out)

# Does the score at timepoint 3 predict drop-out at timepoint 4?
mod <- glm(dropped.out ~ dat1, family = 'binomial')
summary(mod)  
# Reject null hypothesis; timepoint 3 predicts dropout at next timepoint
# Evidence against MCAR drop-out



#  OLS Regression Model 
  mod.ols2 <- lm(Y_mar ~ Group + Time + Group*Time, data = dat)

# MMRM
  mod.gls2 <- gls(Y_mar ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)

# Marginal Contrasts:
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  
  # Examine output:
  out$Beta
  # Look at drop-out:
  aggregate(Y_mar ~ Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  aggregate(Y_mar ~ Group + Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  # Substantial differential rates of drop-out across the two treatment arms

  # Descriptive Statistics 
  # Substantial bias in the descriptive statistics when MAR dropout
  aggregate(Y_comp ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
  aggregate(Y_mar ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
  # Do the arithmetic with the decsriptive statistics, compare to OLS:
  
  # OLS likewise has substantial bias with MAR dropout:
  mod.emms.ols1$contrasts
  mod.emms.ols2$contrasts
  # Notice that the OLS contrasts are equal to those using descriptive statistics

  # MMRM is correctly specified model:
  mod.emms.gls1$contrasts
  mod.emms.gls2$contrasts
  # MMRM unbiased with MAR drop-out
  # In this single replication, we see very little bias
  # Not bad for 30% drop-out!



# ----------------------------------------------------------------------------------------------------
# MNAR example  
# MNAR with correlation = 0.8 (standard interpretation of MNAR - nothing works!)

  dat <- out$dat
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mnar'), 
                 prop.miss = 0.3)

# Test for MCAR 
# Does your score at a previous timepoint predict your score at the next timepoint?
dat1 <- dat[dat$Time == 'Time_3', 'Y_mnar']
dat2 <- dat[dat$Time == 'Time_4', 'Y_mnar']
dropped.out <- !is.na(dat1) & is.na(dat2) # these subjects dropped out at timepoint 4
table(dropped.out)

# Does the score at timepoint 3 predict drop-out at timepoint 4?
mod <- glm(dropped.out ~ dat1, family = 'binomial')
summary(mod)  
# Reject null hypothesis; timepoint 3 predicts dropout at next timepoint
# Evidence against MCAR drop-out

 
#  OLS Regression Model 
  mod.ols2 <- lm(Y_mnar ~ Group + Time + Group*Time, 
                            data = dat)

# MMRM
  mod.gls2 <- gls(Y_mnar ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)  
  
# Marginal Contrasts:
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
 
  # Examine output:
  out$Beta
  # Look at drop-out:
  aggregate(Y_mnar ~ Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  aggregate(Y_mnar ~ Group + Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  # Substantial differential rates of drop-out across the two treatment arms

  # Descriptive Statistics 
  # Substantial bias in the descriptive statistics when MAR dropout
  aggregate(Y_comp ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
  aggregate(Y_mnar ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)

  # OLS likewise has substantial bias with MAR dropout:
  mod.emms.ols1$contrasts
  mod.emms.ols2$contrasts
  # 25% bias!

  # MMRM is correctly specified model:
  mod.emms.gls1$contrasts
  mod.emms.gls2$contrasts
  # MMRM still biased with MNAR drop-out
  # However, it's not AS BAS as OLS and descriptive statistics
  

