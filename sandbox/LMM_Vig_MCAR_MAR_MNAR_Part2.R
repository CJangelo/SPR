# Continuation

# Let's examine the relationship between MCAR, MAR, and MNAR more in-depth
# These aren't clear delineations, even with simulations, which by definition are clear and neat and tidy
# In the real world, the distinctions between these types of drop-out are even murkier

# However, thinking about it more:
# 4. MAR with correlation = 0.1  (Approaches MCAR!)
# 5. MNAR with correlation = 0.95  (approaches MAR!)
# 6. MNAR with correlation = 0.0 (you're fucked!)

# Take-away: MNAR is ugly, but if correlations are high, you're good!

# Run everythign through, save the output, then do the knitr, don't want to re-run it every damn time!

#---------------------------------------------------------------------------------------
#
#
#
#
#
#   VIGNETTE: HOW DO I HANDLE MISSING DATA? Part 2 - What is the relationship between MCAR/MAR/MNAR?



# --------------------------------------------------------------------------------
# MNAR with correlation = 0.95  (approaches MAR!)
# We just saw that the MMRM helps the estimates of MNAR drop-out
# Let's explore how much

# INcrease the correlation to 0.95, AKA, it's the same across timepoints
 set.seed(03232021)

  out <- sim_dat(N = 2000, 
                 number.groups = 2 , 
                 number.timepoints = 4, 
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 2,
                 corr = 'ar1', 
                 cor.value = 0.95, 
                 var.values = 2)
                
  
  dat <- out$dat
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mnar'), 
                 prop.miss = 0.3)

#  OLS Regression Model 
  mod.ols1 <- lm(Y_comp ~ Time + Group + Time*Group, 
                 data = dat)

  mod.ols2 <- lm(Y_mnar ~ Group + Time + Group*Time, 
                            data = dat)

#   MMRM 
  library(nlme)
  mod.gls1 <- gls(Y_comp ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)


  mod.gls2 <- gls(Y_mnar ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)

  
  # Compute marginal mean contrasts
  # Don't worry about the degrees of freedom, we aren't focused on the Type I error or Power
  # Just look at the parameter estimates to check for bias
  library(emmeans)
  mod.emms.ols1 <- emmeans(mod.ols1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

  mod.emms.gls1 <- emmeans(mod.gls1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

 # Examine output:
  out$Beta
  # Look at drop-out:
  aggregate(Y_mnar ~ Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  aggregate(Y_mnar ~ Group + Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  # Substntially different drop out rates across the treatment arms 
  
  # Descriptive Statistics 
  # Sufficient with MCAR drop out - both complete and missing align, yield unbiased estimates:
  aggregate(Y_comp ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
  aggregate(Y_mnar ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
# > 1.15+.29
# [1] 1.44
# > 1.44/2
# [1] 0.72
  
  # OLS likewise has substantial bias in the parameter estimates with MNAR drop out, as we expect
  mod.emms.ols1$contrasts
  mod.emms.ols2$contrasts
 
  mod.emms.gls1$contrasts
  mod.emms.gls2$contrasts
  # MMRM yields decent estimates - still some bias remaining, but much less
  # This suggests that the MMRM can still help you 
  out$cor.mat
  # especially if your residuals are highly correlated 
  # Hypothesis: if your correlations are very high, then your response at timepoint t is very well predicted by 
  # response at timepoint t-1;
  # Thus, whether or not you drop out at timepoint t is very well predicted by your response at timepoint t-1
  # The latter is the definition of MAR
  # We can see that the distinction between MAR and MNAR becomes murkier as the correlation across timepoints increases
  
# --------------------------------------------------------------------------------
# MNAR with correlation = 0.05 (you're fucked)
# We can further explore this hypothesis:
 set.seed(03232021)

  out <- sim_dat(N = 2000, 
                 number.groups = 2 , 
                 number.timepoints = 4, 
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 2,
                 corr = 'ar1', 
                 cor.value = 0.05, 
                 var.values = 2)
                
  
  dat <- out$dat
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mnar'), 
                 prop.miss = 0.3)  
 
#  OLS Regression Model 
  mod.ols1 <- lm(Y_comp ~ Time + Group + Time*Group, 
                 data = dat)

  mod.ols2 <- lm(Y_mnar ~ Group + Time + Group*Time, 
                            data = dat)

#   MMRM 
  library(nlme)
  mod.gls1 <- gls(Y_comp ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)


  mod.gls2 <- gls(Y_mnar ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)

  
  # Compute marginal mean contrasts
  # Don't worry about the degrees of freedom, we aren't focused on the Type I error or Power
  # Just look at the parameter estimates to check for bias
  library(emmeans)
  mod.emms.ols1 <- emmeans(mod.ols1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

  mod.emms.gls1 <- emmeans(mod.gls1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

  
# Examine output:
  out$Beta
  # Look at drop-out:
  aggregate(Y_mnar ~ Group + Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  # Substntially different drop out rates across the treatment arms 
  
  # Descriptive Statistics 
  # Sufficient with MCAR drop out - both complete and missing align, yield unbiased estimates:
  aggregate(Y_comp ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
  aggregate(Y_mnar ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
# Again, not good
  
  # OLS likewise has substantial bias in the parameter estimates with MNAR drop out, as we expect
  mod.emms.ols1$contrasts
  mod.emms.ols2$contrasts
  # bad again
 
  mod.emms.gls1$contrasts
  mod.emms.gls2$contrasts
  # MMRM doesn't help you at all. Estimates just as bad as the OLS
  out$cor.mat  
  # This correlation matrix explains why - you get NO help in what your score would've been at timepoint t
  # from your score at timepoint t-1. 
  # This is the most extreme version of MNAR. 
  # There's not much you can do here. 

    
#--------------------------------------------------------------------------------------------------------
# MAR with correlation = 0.1  (Approaches MCAR!)
# 
# if our hypothesis about the correlations is accurate, then we should be able to predict what happens here as well
# MAR dropout with a correlation approaching zero (0.05, let's say) will be approximately MCAR drop out
# In other words, the bias in descriptive statistics with MAR dropout will go away 
  
  dat <- out$dat
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mar'), 
                 prop.miss = 0.3)  
 
#  OLS Regression Model 
  mod.ols2 <- lm(Y_mar ~ Group + Time + Group*Time, data = dat)

#   MMRM 
  library(nlme)
  mod.gls2 <- gls(Y_mar ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)

  # Compute marginal mean contrasts
  # Don't worry about the degrees of freedom, we aren't focused on the Type I error or Power
  # Just look at the parameter estimates to check for bias
  library(emmeans)
  mod.emms.ols1 <- emmeans(mod.ols1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

  mod.emms.gls1 <- emmeans(mod.gls1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

  
# Examine output:
  out$Beta
  # Look at drop-out:
  aggregate(Y_mar ~ Group + Time, FUN = function(x) mean(is.na(x)), dat = dat, na.action = na.pass)
  # Substantially different drop out rates across the treatment arms 
  
  # Descriptive Statistics 
  # Sufficient with MCAR drop out - both complete and missing align, yield unbiased estimates:
  aggregate(Y_comp ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
  aggregate(Y_mar ~ Group + Time, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
# No bias - looks good to go
  
  # OLS likewise has substantial bias in the parameter estimates with MNAR drop out, as we expect
  mod.emms.ols1$contrasts
  mod.emms.ols2$contrasts
  # No bias! Good to go - similar to MCAR, you have a little noise due to 30% dropout 
 
  mod.emms.gls1$contrasts
  mod.emms.gls2$contrasts
  # Again, we see this turning out like MCAR dropout
  out$cor.mat  

# Our hypothesis seems validated - it looks like as correlations go to zero under MAR dropout, 
# the dropout type approaches MCAR

  