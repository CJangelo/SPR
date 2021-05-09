
## 3.23.21 - Vignette 1
# TYpe I error and Power
# Two different Vignettes
# Upload the RData to github as backup
# Use Rdata to create the Vignette file
# TODO: make the Output look clean, you need the tables with the rejection rate to look nice
# use the R2Word package with it
# Then dump the bootstrap stuff up there as a TODO vignette - need to decide which is the best MMRM approach!
# I think it's gls(), not lmer()
# Can change that later, not crucial for right now - R packages typically have functionality for a bootstrap anyway
# DO all this then switch to Alcon


rm(list = ls())
gc()

library(MASS)
library(glmmTMB) # https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
library(emmeans)
library(nlme)
library(lme4)


###
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/sim_dat.R")
source("C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/dropout.R")


  # N = 100 #this should be divisible by however many groups you use!
  # number.groups <- 2
  number.timepoints <- 4
  # type_dropout = c('mar')  
  # reg.formula <- formula( ~ Time + Group + Time*Group)

  test.sim <- F
  all.p <- vector()
  all.est <- vector()
  start.time <- Sys.time()
  repl <- 1
  
  
for(repl in 1:100){
  
  set.seed(03232021 + repl)
  set.seed(17)
  out <- sim_dat(N = 40, 
                 number.groups = 2 , 
                 number.timepoints = 4, 
                 reg.formula = formula( ~ Time + Group + Time*Group),
                 Beta = 0,
                 corr = 'ar1', 
                 cor.value = 0.8, 
                 var.values = 1)
                
  
  dat <- out$dat
  dat <- dropout(dat = dat, 
                 type_dropout  = c('mar'), 
                 prop.miss = 0.3)
  

#  OLS Regression Model 
  mod.ols1 <- lm(Y_comp ~ Time + Group + Time*Group, 
                 data = dat)

  mod.ols2 <- lm(Y_mar ~ Group + Time + Group*Time, 
                            data = dat)


if (test.sim == FALSE){
  
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
  mod.lmer1 <- lmer(Y_comp ~ Group + Time + Group*Time + ( 0 + Time | USUBJID), 
                    data = dat, 
                    control = lmerControl(check.nobs.vs.nRE = 'ignore')) 
  
  
  mod.lmer2 <- lmer(Y_mar ~ Group + Time + Group*Time + (0 + Time | USUBJID), 
                    data = dat, 
                    control = lmerControl(check.nobs.vs.nRE = 'ignore')) 
  
}

  # Compute marginal means
  library(emmeans)
  mod.emms.ols1 <- emmeans(mod.ols1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.ols2 <- emmeans(mod.ols2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

if (test.sim == FALSE) {

  mod.emms.gls1 <- emmeans(mod.gls1, pairwise ~ Group | Time, adjust = 'none',  mode = 'satterthwaite')
  mod.emms.gls2 <- emmeans(mod.gls2, pairwise ~ Group | Time, adjust = 'none',  mode = 'satterthwaite')

  mod.emms.us1 <- emmeans(mod.us1, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')
  mod.emms.us2 <- emmeans(mod.us2, pairwise ~ Group | Time, adjust = 'none',  mode = 'df.error')

  # Kenward Rogers Degrees of Freedom
  mod.emms.lmer1 <- emmeans(mod.lmer1, pairwise ~ Group | Time, adjust = 'none', mode = "kenward" )
  mod.emms.lmer2 <- emmeans(mod.lmer2, pairwise ~ Group | Time, adjust = 'none', mode = "kenward" )

}
  
  # Post Processing
    est <- vector()
    p <- vector()
    la <- vector()

    
if (test.sim == TRUE) { 
    mod.list <- c('ols1', 'ols2') 
} else {
    mod.list <- c('ols1', 'ols2', 'gls1', 'gls2', 'us1', 'us2', 'lmer1', 'lmer2')
}
  
    
  for(nn in mod.list){
   
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
  ## save.image(file = 'out_Sim_Study_MMRM_TypeI_23Mar2021.RData')
  ## save.image(file = 'out_Sim_Study_MMRM_Power_23Mar2021.RData')

  
  # Rejection Rate of the individual tests/predictors
  colMeans(all.p < 0.05)
  round(colMeans(all.est), 2)

  # Family-Wise - The FWE is key part of evaluating Type I error control
  ##for(nn in c('ols1', 'ols2', 'us1', 'us2', 'lmer1', 'lmer2')){
    out.FWE <- vector()
    mn <- unique((unlist(lapply(colnames(all.p), function(x) strsplit(x, '_')[[1]][1]))))
  for (nn in mn) {
      pcomp <- all.p[ , grep(nn, colnames(all.p))] 
        print( table(apply(pcomp < 0.05, 1, function(x) sum(x) > 0)) )
        out.FWE <- rbind(out.FWE, 
                         mean(apply(pcomp < 0.05, 1, function(x) sum(x) > 0)) )
        #cat(nn ,table(apply(pcomp < 0.05, 1, function(x) sum(x) > 0)), '\n' )
    }

data.frame('Model' = mn, out.FWE)
#
  # Common multiplicity correction is Benjamini & Hochberg’s FDR
  # Evaluate how FDR affects the power
    
### 
  ## for(nn in c('ols1', 'ols2', 'us1', 'us2', 'gls1', 'gls2', 'lmer1', 'lmer2')){
  #labels <- vector()
    all.out <- vector()
  
  for (nn in mn) {
      out <- all.p[ , grep(nn, colnames(all.p))] 
      out.FDR <- matrix(0, nrow = nrow(out), ncol = ncol(out))
        #labels <- rbind(labels, colnames(out))
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
        all.out <- rbind(all.out, colMeans(out.FDR) )
      #print(colMeans(out.FDR))
      print(table( apply(out.FDR, 1, function(x) sum(x) > 0) ))

  }# end loop over model types

# Power after FDR Correction:
  data.frame('Model' = mn, all.out)
# FDR works well - next step is to compare to bootstrap approach
