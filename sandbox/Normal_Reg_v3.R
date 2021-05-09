


# Normal Distribution - Regression 
# Works, g2g
# Not interesting, v1 is cleaner and clearer
# 2.19.21


rm(list = ls())
gc()

library(MASS)

N = 1e4 #this should be divisible by however many groups you use!
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
    sigma <- 1
    
    
    # Parameters:
    XB <- X %*% Beta
    dat$XB <- as.vector(XB)
    error <- rnorm(n = N, mean = 0, sd = sigma)
    dat$Y <- dat$XB + error
    
    # check
    hist(dat$Y)
    aggregate(Y ~ 1, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = var, data = dat)
    
    # Model?
    mod <-lm(Y ~ Group, data = dat)
    #summary(mod)
    cbind(Beta, coef(mod))
    # #####
  

vP <- c('b0' = 0.2, 'b1' = 1) #, 'sigma' = 1)

Normal_loglike <- function(vP, dat){
  
  
  # the 'par' (parameters to estimate) are only passed as a vector
  Beta.hat <- vP[c('b0', 'b1')]
  Beta.hat <- matrix(Beta.hat, ncol = 1)
  #sigma <- vP['sigma']
  
  out <- aggregate(Y ~ Group, FUN = function(x) 'var' = var(x), data = dat)
  s2 <- out$Y # observed
  
  out <- aggregate(Y ~ Group, FUN = function(x) 'mean' = mean(x), data = dat)
  xbar <- out$Y # observed 
  
  X <- model.matrix( ~ Group , data = dat) 
  dat$mu <- as.vector(X %*% Beta.hat)
  
  out <- aggregate(mu ~ Group, FUN = function(x) 'mean' = mean(x), data = dat)
  mu <- out$mu

  out <- aggregate(Y ~ Group, FUN = function(x) x, data = dat)
  out <- out[ , -1]
  Yg <- t(out)
  tmp <- Yg - matrix(xbar, nrow = nrow(Yg), ncol = ncol(Yg), byrow = T)
  sigma2.hat <- apply(tmp, 2, function(x) sum(x^2)/(N-2))
  
  o <- rbind(xbar, s2)
  e <- rbind(mu, sigma2.hat)

  #LL <- sum((o-e)^2)
  LL <- sum((xbar - mu)^2) # Don't even need the sigma...
  # LL <- (o * log(e))
  # LL <- -1*LL
  # LL <- sum(LL)
  return(LL)
  
}

# optimize
out <- optim(par = vP, fn = Normal_loglike, method = 'BFGS', dat = dat)
cbind(
  'Gen param' = c(Beta, sigma),
  'Custom' = c(out$par, NA), # N-p
  'R package' = c(coef(mod), summary(mod)$sigma^2)
)