

# Normal-distribution Regression
# Works, g2g
# Minimize SSE
# 2.19.21
# TODO: pull the hessian and estimate SE

rm(list = ls())
gc()

library(MASS)

N = 100 #this should be divisible by however many groups you use!
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
    sigma2 <- 9 
    
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
    # g2g
  
    

vP <- c('b0' = 0.2, 'b1' = 1)

Normal_loglike <- function(vP, dat){
  
  Beta.hat <- vP[c('b0', 'b1')]
  Beta.hat <- matrix(Beta.hat, ncol = 1)
  X <- model.matrix( ~ Group , data = dat) 
  Y.hat <- X %*% Beta.hat 
  SSE <- sum((Y.hat - dat$Y)^2)
  return(SSE)

}



# optimize
vP <- c('b0' = 0.2, 'b1' = 1)
out <- optim(par = vP, fn = Normal_loglike, method = 'BFGS', dat = dat)

# Use Beta.hat to estimate the sigma:
Beta.hat <- out$par
Beta.hat <- matrix(Beta.hat, ncol = 1)
X <- model.matrix( ~ Group , data = dat) 
Y.hat <- X %*% Beta.hat 
sigma2.hat <- sum((Y.hat - dat$Y)^2)/(N-2)  # N-p-1, where p is number of predictors

# Compare observed vs predicted:
mean(Y.hat)
mean(dat$Y)
# Compare parameter estimates:
cbind(
  'Gen param' = c(Beta, sigma2), 
  'Custom Function' = c(out$par, sigma2.hat), 
  'R package' = c(coef(mod), summary(mod)$sigma^2)
)

