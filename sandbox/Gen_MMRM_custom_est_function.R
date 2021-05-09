
########################################################
#
#
#           GENERATE LONGITUDINAL DATA
#
#
############################################################

rm(list = ls())
gc()


N = 100 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 4


    dat <- data.frame(
                      'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                      'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                      'Y_comp' = rep(NA, N*number.timepoints), 
                      #'Bio' = rep(rnorm(N, mean = 0, sd = 1), number.timepoints),
                      'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
                      stringsAsFactors=F)
    
    # Create Beta parameters for these design matrix:
    X <- model.matrix( ~ Group + Time + Time*Group , data = dat) 
    
    Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    # Small separation between txa:
    #Beta[grepl('Group_2', rownames(Beta)) & grepl('Time', rownames(Beta)), ] <- seq(0.25, 1, length.out = number.timepoints - 1)
    Beta[grepl('Time', rownames(Beta)), ] <- seq(0.25, 1, length.out = sum(grepl('Time', rownames(Beta))))
    # Matrix multiply:
    XB <- X %*% Beta
    ########################################################################
###
    # Distribution of residuals:
    mu <- rep(0, number.timepoints)
    cor.mat <- diag(1, nrow = number.timepoints, ncol = number.timepoints) 
    for(i in 1:number.timepoints){
      for(j in 1:i){
        
        cor.mat[i , j] <- 0.8^(i -j) # exponential decay - very high correlations
        cor.mat[j, i] <- cor.mat[i, j]
        
      }
    }
    
    # Sigma 
    var.mat <- seq(1, 2, length.out = number.timepoints)
    var.mat <- diag(sqrt(var.mat), nrow = number.timepoints, ncol = number.timepoints) 
    sigma <- var.mat %*% cor.mat %*% var.mat
    # Simulate the errors:
    error <- MASS:::mvrnorm(n = N, mu = mu, Sigma = sigma)
    # Re-arrange in long format:
    error.long <- vector()
    for(time in 1:number.timepoints){
      error.long <- rbind(error.long, 
                              error[, time, drop = F]
      )
    }
    
    Y <- XB + error.long
    dat$Y_comp <- as.vector(Y)

#-------------------------------------------------------------------
    # Fit Model:
    library(mvtnorm)
    source('C:/Users/ciaconangelo/Documents/RESEARCH_NEW_LAPTOP/R_CODE_Long_Mixed_Models/Fit_MMRM_custom_est_function.R')

    mod <- custom_MMRM_estimation(Y_comp ~ Group + Time + Time*Group, subject_id = 'USUBJID',  time_var = 'Time',  dat = dat)

    mod.ols <- lm(Y_comp ~ Group + Time + Time*Group, dat)
    # summary(mod.ols)
    # sqrt(diag(vcov(mod.ols)))
    # coef(mod.ols)
    
    cbind('True_Beta' = Beta, 
          'lm()_Beta' = coef(mod.ols), 
          'Custom_Beta' = mod$Beta.hat, 
          'lm()_SE' = sqrt(diag(vcov(mod.ols))), 
          'Custom_SE' = mod$SE)
    sigma
    mod$sigma.hat
    
    