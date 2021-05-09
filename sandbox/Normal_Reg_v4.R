

# Normal-distribution Regression
# 3.16.21 - Need to extend to repeated measures
# Py <- dnorm(ee, mean = 0, sd = 10) 
# Py <- dnorm(ee, mean = 0, sd = sqrt(s2))
# Why does it work either way? 
# Now extend to repeated measures - except we want the marginal model, 
# not conditional model 

  
rm(list = ls())
gc()

library(MASS)

N = 1000 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
#set.seed(2182021)



    dat <- data.frame(
                      'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                      'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                      'Y_comp' = rep(NA, N*number.timepoints), 
                      #'Bio' = rep(rnorm(N, mean = 0, sd = 1), number.timepoints),
                      'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
                      stringsAsFactors=F)
    
        
        # Design Matrix
    X <- model.matrix( ~ Group , data = dat) 
    Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    Beta[] <- c(0.2, 1)
    sigma2 <- 3 
    
    
    # Parameters:
    XB <- X %*% Beta
    dat$XB <- as.vector(XB)
    error <- rnorm(n = N, mean = 0, sd = sqrt(sigma2))
    dat$Y <- dat$XB + error
    
    # check
    hist(dat$Y)
    aggregate(Y ~ 1, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = var, data = dat)
    
    # Model?
    mod <-lm(Y ~ Group, data = dat)
    summary(mod)
    summary(mod)$sigma
    summary(mod)$sigma^2
    # g2g
  
    
#-------------------------------------------------------------  
Normal_loglike <- function(vP, dat){
  
# Model-based prob: 
  # s2 <- exp(vP['sigma2'])
  s2 <- vP['sigma2']
  Beta.hat <- vP[c('b0', 'b1')]
    Beta.hat <- matrix(Beta.hat, ncol = 1)
  X <- model.matrix( ~ Group , data = dat) 
  Y.hat <- X %*% Beta.hat 
  Y <- matrix(dat$Y, ncol = 1)
  ee <- Y - Y.hat
  #Py <- dnorm(ee, mean = 0, sd = 1) 
  Py <- dnorm(ee, mean = 0, sd = sqrt(s2))
  LL <- -1*sum(log(Py), na.rm = T)
  return(LL)

}
# -----------------------------------------------------



    
vP <- c('b0' = 0, 'b1' = 0, 'sigma2' = 1)
delta <- 1

while(delta > 0.0001){

  est0 <- vP

  # M -step:
    out <- optim(par = vP, fn = Normal_loglike, method = 'BFGS', dat = dat)
    Beta.hat <- out$par[c('b0', 'b1')]
    Beta.hat <- matrix(Beta.hat, ncol = 1)
    X <- model.matrix( ~ Group , data = dat) 
    Y.hat <- X %*% Beta.hat 
    Y <- matrix(dat$Y, ncol = 1)
    ee <- Y - Y.hat
    SSE <- sum(ee^2)
  # E-step update:
    vP[c('b0', 'b1')] <- out$par[c('b0', 'b1')]
    vP['sigma2'] <- SSE/(N-1)
    delta <- max(abs(est0 - vP))
    est0 <- vP
  
} #

vP

# Compare parameter estimates:
cbind(
  'Gen param' = c(Beta, sigma2), 
  'R package' = c(coef(mod), summary(mod)$sigma^2), 
  'Custom Function' = vP 
)

