

## Simulation study to evaluate Power of MMRM 
# 2.25.21 - TODO resolve issues with KR df, comparison to OLS df
# TODO incorporate the bootstrap as another comparison


  #   load(file = 'out_Sim_Study_MMRM_Power_v2.RData')

rm(list = ls())
gc()

library(MASS)
library(glmmTMB) # https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
library(emmeans)
library(nlme)
library(lme4)
# library(foreach)
# library(doParallel)

# Initialize parallelization
# cl<-makeCluster(8, outfile="useless.txt")
# print(cl)
# registerDoParallel(cl)
# getDoParWorkers()
# getDoParName()

###
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/Function_Generate_Long_Data.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/Function_Implement_Missingness_v2.R")
#source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/Function_bootstrap_MMRM_lmer.R")
source('C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/Fit_MMRM_custom_est_function.R')


  N = 500 #this should be divisible by however many groups you use!
  number.groups <- 2
  number.timepoints <- 3
  score = 'Y_mar'   #dat$Y_mcar <- dat$Y_cmcar <- dat$Y_mar <- dat$Y_mnar <- dat$Y_comp
  #number.boot.repl <- 1000
  all.p <- vector()
  all.est <- vector()
  start.time <- Sys.time()
  repl <- 4
  
  
#for(repl in 1:100){
  
  set.seed(03172021 + repl)
  
  dat <- Generate_Long_Data(N = N, number.groups = number.groups , number.timepoints = number.timepoints)
  dat <- Implement_Missingness(dat = dat, score = score, N = N, number.timepoints = number.timepoints)

#  OLS Regression Model 
  mod.ols1 <- lm(Y_comp ~ Group + Time + Group*Time, data = dat)
  mod.ols2 <- lm(as.formula(paste0(score, '~ Group + Time + Group*Time')), data = dat)


# MMRM - library(glmmTMB)
  mod.us1 <- glmmTMB(Y_comp ~ Group + Time + Group*Time + us(Time + 0 | USUBJID), 
                  data=dat, 
                  REML = T,
                  dispformula=~0)
  
  mod.us2 <- glmmTMB(as.formula(paste0(score, '~ Group + Time + Group*Time + us(Time + 0 | USUBJID)')),
                  data=dat,
                  REML = T,
                  dispformula=~0)
  
  
  # MMRM - library(lme4)
  #  https://bit.ly/2ZSxBRG
  #  https://drive.google.com/file/d/1sOZUAFOc004H4jO8vuUc_4HyYHEgu45b/view?usp=sharing&usp=embed_facebook
  mod.lmer1 <- lmer(Y_comp ~ Group + Time + Group*Time + ( 0 + Time | USUBJID), data = dat, 
                          control = lmerControl(check.nobs.vs.nRE = 'ignore')) 
  
  
  mod.lmer2 <- lmer(as.formula(paste0(score, '~ Group + Time + Group*Time + (0 + Time | USUBJID)')), data = dat, 
                          control = lmerControl(check.nobs.vs.nRE = 'ignore')) 
  

#   MMRM - 
  #library(nlme)
  mod.mmrm1 <- gls(Y_comp ~ Group + Time + Group*Time,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)


  mod.mmrm2 <- gls(as.formula(paste0(score, '~ Group + Time + Group*Time')),
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)

  
  mod.custom1 <- custom_MMRM_estimation(Y_comp ~ Group + Time + Time*Group, 
                                       subject_id = 'USUBJID',  
                                       time_var = 'Time',  
                                       var_matrix = 'unstructured',
                                       dat = dat)
  
  # mod.formula <- as.formula(paste0(score, '~ Group + Time + Group*Time'))
  # var_matrix = 'independent'
  
  
  mod.custom2 <- custom_MMRM_estimation(as.formula(paste0(score, '~ Group + Time + Group*Time')), 
                                        subject_id = 'USUBJID', 
                                        time_var = 'Time', 
                                        var_matrix = 'unstructured', 
                                        dat = dat)

  # mod.custom3 <- custom_MMRM_estimation(as.formula(paste0(score, '~ Group + Time + Group*Time')), 
  #                                       subject_id = 'USUBJID', 
  #                                       time_var = 'Time',  
  #                                       var_matrix = 'independent',
  #                                       dat = dat)

 
  cbind(
        'OLS' = coef(mod.ols1),
        'LMM_lme4' = fixef(mod.lmer1),
        'MMRM_nlme' = coef(mod.mmrm1),
        'MMRM_glmmTMB' = fixef(mod.us1)$cond,
        'Custom' = as.vector(mod.custom1$Beta.hat)
  )

  cbind(
        'OLS' = coef(mod.ols2),
        'LMM_lme4' = fixef(mod.lmer2),
        'MMRM_nlme' = coef(mod.mmrm2),
        'MMRM_glmmTMB' = fixef(mod.us2)$cond,
        'Custom' = as.vector(mod.custom2$Beta.hat))
       # 'Custom_ind' = as.vector(mod.custom3$Beta.hat)
       # 'Custom_comp' = mod.custom4$Beta.hat
    
  
  G1 <- matrix(c(1, 0, 1, 1, 0, 0), ncol = nrow(mod.custom2$Beta.hat))
  G2 <- matrix(1, ncol = nrow(mod.custom2$Beta.hat))
  C  <- rbind(G1, G2)
  sum(C %*% coef(mod.ols1))
  sum(C %*% fixef(mod.lmer1))
  
  
  sum(C %*% coef(mod.ols2))
  sum(C %*% coef(mod.mmrm2))
  sum(C %*% fixef(mod.lmer2))
  sum(C %*% fixef(mod.us2)$cond)
  sum(C %*% mod.custom2$Beta.hat)
  #sum(C %*% mod.custom3$Beta.hat)

mod.custom2$sigma.hat
