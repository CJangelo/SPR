
## 3.21.21 - Vignette 1
# Quick example

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
#   VIGNETTE: WHAT TO DO IF YOU'RE ASKED TO FIT A MMRM IN R?


rm(list = ls())
gc()

library(MASS)
library(glmmTMB) # https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
library(emmeans)
library(nlme)
library(lme4)
# library(lmerTest)

###
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/dropout.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/VarCov_gls.R")


  set.seed(03212021)
  
  out <- sim_dat(N = 100, 
                 number.groups = 2 , 
                 number.timepoints = 4, 
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 0,
                 corr = 'ar1', 
                 cor.value = 0.8, 
                 var.values = 2)
                
  
  dat <- out$dat
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mar'), 
                 prop.miss = 0.3)
  

#  OLS Regression Model 
  mod.ols1 <- lm(Y_comp ~ Time + Group + Time*Group, 
                 data = dat)

  mod.ols2 <- lm(Y_mar ~ Group + Time + Group*Time, 
                            data = dat)


#   MMRM 
  library(nlme)
  mod.gls1 <- gls(Y_comp ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)


  mod.gls2 <- gls(Y_mar ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)
  
  
# MMRM - 
  library(glmmTMB)
  mod.us1 <- glmmTMB(Y_comp ~ Group + Time + Group*Time + us(Time + 0 | USUBJID), 
                  data=dat, 
                  REML = T,
                  dispformula=~0)
  
  mod.us2 <- glmmTMB(Y_mar ~ Group + Time + Group*Time + us(Time + 0 | USUBJID),
                  data=dat,
                  REML = T,
                  dispformula=~0)
  
  
  # MMRM - 
  library(lme4)
  #  https://bit.ly/2ZSxBRG
  #  https://drive.google.com/file/d/1sOZUAFOc004H4jO8vuUc_4HyYHEgu45b/view?usp=sharing&usp=embed_facebook
  mod.lmer1 <- lmer(Y_comp ~ Group + Time + Group*Time + ( -1 + Time | USUBJID), 
                    data = dat, 
                    control = lmerControl(check.nobs.vs.nRE = 'ignore')) 
  
  
  mod.lmer2 <- lmer(Y_mar ~ Group + Time + Group*Time + (0 + Time | USUBJID), 
                    data = dat, 
                    control = lmerControl(check.nobs.vs.nRE = 'ignore')) 
  



  # Compute marginal means
  library(emmeans)
  mod.emms.ols1 <- emmeans(mod.ols1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

  mod.emms.gls1 <- emmeans(mod.gls1, pairwise ~ Group | Time, adjust = 'none',  mode = 'satterthwaite')
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'satterthwaite')

  mod.emms.us1 <- emmeans(mod.us1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.us2 <- emmeans(mod.us2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')


  # Kenward Rogers Degrees of Freedom
  mod.emms.lmer1.kr <- emmeans(mod.lmer1, pairwise ~ Group | Time, adjust = 'none', mode = "kenward" )
  mod.emms.lmer2.kr <- emmeans(mod.lmer2, pairwise ~ Group | Time, adjust = 'none', mode = "kenward" )
# Does not align with SAS output
  
  # Satterthwaite Degrees of Freedom:
  mod.emms.lmer1.sat <- emmeans(mod.lmer1, pairwise ~ Group | Time, adjust = 'none', mode = 'satterthwaite' )
  mod.emms.lmer2.sat <- emmeans(mod.lmer2, pairwise ~ Group | Time, adjust = 'none', mode = 'satterthwaite' )


  # Examine output:
  out$Beta
  mod.emms.ols1$contrasts
  mod.emms.gls1$contrasts
  mod.emms.us1$contrasts
  mod.emms.lmer1.kr$contrasts
  mod.emms.lmer1.sat$contrasts
 
  # Compare the variance-covariance matrices:
  out$sigma # Generating
  VarCov_gls(mod.gls1) # gls()
  summary(mod.ols1)$sigma
  glmmTMB::VarCorr(mod.us1)$cond
  
  # The lmer() hack yields an unidentified variance-covariance matrix:
  VarCorr(mod.lmer1)
  glmmTMB::VarCorr(mod.us1)
  # https://stat.ethz.ch/pipermail/r-sig-mixed-models/2015q2/023520.html
  
  
  
# Print to console the different model estimates:
  
for(nn in c('ols1', 'ols2', 'gls1', 'gls2', 'us1', 'us2', 'lmer1.kr', 'lmer2.kr', 'lmer1.sat', 'lmer2.sat')){
    
    mod.out <- get(paste0('mod.emms.', nn))
    mod.out <- as.data.frame(mod.out$contrasts)
    cat(paste0('Model: ', nn, '\n'))
    print(mod.out)
    
  }
    
 