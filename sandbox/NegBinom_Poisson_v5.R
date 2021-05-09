
# Custom estimation code
# Estimates Poisson model 
# This works

# 2.9.21 - this Poisson estimation routine reduces the dimensions
# Its' just two groups by two means - should be faster
# It's not faster, but it mimics a typical EM algorithm approach for Latent variable model

# Custom written estimation routine:
Poisson_loglike_v2 <- function(vP, dat, Y){
  
  # the 'par' (parameters to estimate) are only passed as a vector
  # re-create what you need from that vector
  Beta <- vP[c('b0', 'b1')]
  Beta <- matrix(Beta, ncol = 1)

  X <- model.matrix( ~ Group , data = dat) 
  out <- aggregate(as.formula(paste0(Y, ' ~ Group')), FUN = mean, data = dat)
  Y <- out[ , Y]
  # Parameters:
  XB <- exp(X %*% Beta)
  XB <- aggregate(XB ~ Group, FUN = mean, data = dat)
  XB <- XB$V1
  #print(mu)

  loglike <-  Y * log(XB) - XB - lgamma(Y + 1)
  loglike <- -1*loglike
  loglike <- sum(loglike)
  #print(loglike)
  return(loglike)
  
}


# optimize
#Y <- 'Y_pois'
vP <- c('b0' = 0, 'b1' = .5)
start <- Sys.time()
out <- optim(par = vP, fn = Poisson_loglike_v2, method = 'BFGS', 
             dat = dat, Y = 'Y_pois')

end <- Sys.time()

end - start
# Compare Generating Parameter to R package estimate and the custom code estimate
cbind('**Gen Param**' = c(as.vector(Beta)), '**Custom**' = out$par, '**R package**' = c(coef(mod.pois)))
