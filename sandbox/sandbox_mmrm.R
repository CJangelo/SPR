

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
#library(lmerTest)
library(lme4)


###
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/Function_Generate_Long_Data.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/Function_Implement_Missingness_v2.R")


  N = 100 #this should be divisible by however many groups you use!
  number.groups <- 2
  number.timepoints <- 4
  score = 'Y_mar'   #dat$Y_mcar <- dat$Y_cmcar <- dat$Y_mar <- dat$Y_mnar <- dat$Y_comp
  reg <- formula(~ Group + Time + Time*Group)

  all.p <- vector()
  all.est <- vector()
  start.time <- Sys.time()
  repl <- 1
  
for(repl in 1:100){
  
  set.seed(03182021 + repl)
  
  dat <- Generate_Long_Data(N = N, number.groups = number.groups , number.timepoints = number.timepoints)
  #dat <- Generate_Long_Data(N = N, number.groups = number.groups)

  dat <- Implement_Missingness(dat = dat, score = score, N = N, number.timepoints = number.timepoints)

  #write.table(dat, na = '.', quote = F, sep = ', ', row.names = F, file = paste0('dat_sim_study_repl_', repl, '_18Mar.txt'))

#  OLS Regression Model 
  mod.ols1 <- lm(Y_comp ~ Time + Group + Time*Group, data = dat)

  mod.ols2 <- lm(as.formula(paste0(score, '~ Group + Time + Group*Time')), data = dat)


# MMRM - library(glmmTMB)
  mod.us1 <- glmmTMB(Y_comp ~ Group + us(Time + 0 | USUBJID), 
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
                          REML = T,
                          control = lmerControl(check.nobs.vs.nRE = 'ignore')) 
  
  
  mod.lmer2 <- lmer(as.formula(paste0(score, '~ Group + Time + Group*Time + (0 + Time | USUBJID)')), data = dat, 
                          REML = T,
                          control = lmerControl(check.nobs.vs.nRE = 'ignore')) 
  

  summary(mod.lmer2, ddf = c("Kenward-Roger"))
  summary(mod.lmer2, ddf = c("Satt"))
  summary(mod.lmer2, ddf = c("lme4"))
  
  
#   MMRM - 
  # library(nlme)  
  mod.gls1 <- gls(Y_comp ~ Group,
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)


  mod.gls2 <- gls(as.formula(paste0(score, '~ Group + Time + Group*Time')),
                  data = dat,
                  correlation = corSymm(form = ~ 1 | USUBJID),    #  unstructured correlation
                  weights = varIdent(form = ~ 1 | Time),          #  freely estimate variance at subsequent timepoints
                  na.action = na.exclude)
  

  # Compute marginal means
  mod.emms.ols1 <- emmeans(mod.ols1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

  mod.emms.us1 <- emmeans(mod.us1, pairwise ~ Group, adjust = 'none',  mode = 'df.error')
  mod.emms.us2 <- emmeans(mod.us2, pairwise ~ Group, adjust = 'none',  mode = 'df.error')


  # Kenward Rogers Degrees of Freedom
  mod.emms.lmer1 <- emmeans(mod.lmer1, pairwise ~ Group | Time, adjust = 'none',   # mode = 'df.error')
                                mode = "kenward" )
  mod.emms.lmer2 <- emmeans(mod.lmer2, pairwise ~ Group | Time, adjust = 'none', #mode = 'asymptotic')
                                mode = "satterthwaite" )
  
  pbkrtest::vcovAdj(mod.lmer2)
  L <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0), nrow = 1)
  L <- matrix(c(1, 0, 1, 0, 0, 0, 0, 0), nrow = 1)
  L <- matrix(c(1, 0, 0, 1, 0, 0, 0, 0), nrow = 1)
  L <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0), nrow = 1)
  L <- matrix(c(1, 1, 0, 0, 1, 0, 0, 1, 
                1, 0, 0, 0, 1, 0, 0, 0), nrow = 2, ncol = 8)
  out <- pbkrtest::KRmodcomp(mod.lmer2, L)
  out$test$ddf
  out$test
  

  mod.emms.gls1 <- emmeans(mod.gls1, pairwise ~ Group, adjust = 'none', #mode = 'df.error')
                               mode = "satterthwaite" )
  
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none', # mode = 'df.error')
                                mode = "satterthwaite"  )
  
  # Can't use this with missing values?
  mod.emms.lmer2.satt <- emmeans(mod.lmer2, pairwise ~ Group | Time, adjust = 'none',  #  mode = 'df.error')
                                mode = "satterthwaite" , data = dat)

  
  # Post Processing
    est <- vector()
    p <- vector()
    la <- vector()
  for(nn in c('ols1', 'ols2', 'us1', 'us2', 'lmer1', 'lmer2')){
    
    mod.out <- get(paste0('mod.emms.', nn))
    mod.out <- as.data.frame(mod.out$contrasts)
    est <- c(est, mod.out[ , 'estimate'] )
    p <- c(p, mod.out[ , 'p.value'] )
    la <- c(la, paste0(nn, '_timepoint_', 1:number.timepoints) )
    
  }
    
  names(est) <- names(p) <- la
  all.est <- rbind(all.est, est)
  all.p <- rbind(all.p, p)
  cat(paste0('Replication: ', repl, '\n'))

}# end repl
  
  end.time <- Sys.time()
  end.time - start.time
  
  # Save the Simulation Output
  ## save.image(file = 'out_Sim_Study_MMRM_Power_2.RData')
  ## load(file = 'out_Sim_Study_MMRM_Power_2.RData')

  
  # Rejection Rate of the individual tests/predictors
  colMeans(all.p < 0.05)
  round(colMeans(all.est), 2)

  # Family-Wise - The FWE is key part of evaluating Type I error control
  for(nn in c('ols1', 'ols2', 'us1', 'us2', 'lmer1', 'lmer2')){
      pcomp <- all.p[ , grep(nn, colnames(all.p))] 
        #print( table(apply(pcomp < 0.05, 1, function(x) sum(x) > 0)) )
        cat(nn ,table(apply(pcomp < 0.05, 1, function(x) sum(x) > 0)), '\n' )
    }

  # Common multiplicity correction is Benjamini & Hochberg’s FDR
  # Evaluate how FDR affects the power
    
### 
  for(nn in c('ols1', 'ols2', 'us1', 'us2', 'lmer1', 'lmer2')){
      out <- all.p[ , grep(nn, colnames(all.p))] 
      out.FDR <- matrix(0, nrow = nrow(out), ncol = ncol(out))
        colnames(out.FDR) <- colnames(out)
        alpha0 <- 0.05
        
        for(repl in 1:nrow(out)){
            # Order the p-values
            tmp <- sort(out[repl,], decreasing=T) # step-up test proceeds from least to most significant p-value
            #eq15 <- alpha0 * (1/(1 + number.endpoints - number.endpoints:1)) # Eq 14 in the paper (Hochberg step-up test)
            eq15 <- alpha0 * (number.timepoints:1/number.timepoints) # Eq 15 in the paper (Benjamini & Hochberg’s FDR)
            # Note: Benjamini & Hochberg’s FDR - this is flipped from the way it's written in their 1995 paper, this way I can use match() function
              if(any(tmp <= eq15)){ #if something is stat sig, otherwise, next repl
                  stat.sig.endpoints <- names(tmp)[ match(TRUE, tmp <= eq15):number.timepoints] # first sig p-value and everything smaller, match() will pull the first "TRUE" - "match returns a vector of the positions of (first) matches of its first argument in its second."
                  out.FDR[repl, stat.sig.endpoints ] <- 1
                }# end if statement
  
          }# end loop over replications

        # Summarize
      print(colMeans(out.FDR))
      table( apply(out.FDR, 1, function(x) sum(x) > 0) )

  }# end loop over model types

  
# FDR works well - next step is to compare to bootstrap approach
