

# Custom written code 
# Estimate Neg Binom 
# Estimate Poisson
# Both align with R packages
# g2g
# NB: They use the full response pattern and full Group variable
# No dimension reduction

# Custom written estimation routine:
Neg_Binom_loglike <- function(vP, dat, Y){
  
  # the 'par' (parameters to estimate) are only passed as a vector
  # re-create what you need from that vector
  Beta <- vP[c('b0', 'b1')]
  Beta <- matrix(Beta, ncol = 1)
  size <- vP['theta']

  X <- model.matrix( ~ Group , data = dat) 
  Y <- dat[ , Y]
  # Parameters:
  mu <- exp(X %*% Beta)
  prob <- size/(size + mu)
  
  loglike <- lgamma(Y + size) - lgamma(size) - lgamma(Y + 1) + size*log(prob) + (Y)*log(1-prob)
  loglike <- -1*loglike
  loglike <- sum(loglike)
  return(loglike)
  
}


# optimize
Beta; theta # Gen param
vP <- c('b0' = 0.2, 'b1' = 1, 'theta' = 1)
out <- optim(par = vP, fn = Neg_Binom_loglike, method = 'BFGS', 
             dat = dat, Y = 'Y_nb')
# Compare Generating Parameter to R package estimate and the custom code estimate

cbind('Gen Param' = c(as.vector(Beta), theta), 'Custom' = out$par, 'R package' = c(coef(mod.nb), mod.nb$theta))


# Custom written estimation routine:
Poisson_loglike <- function(vP, dat, Y){
  
  # the 'par' (parameters to estimate) are only passed as a vector
  # re-create what you need from that vector
  Beta <- vP[c('b0', 'b1')]
  Beta <- matrix(Beta, ncol = 1)

  X <- model.matrix( ~ Group , data = dat) 
  Y <- dat[ , Y]
  # Parameters:
  mu <- exp(X %*% Beta)

  loglike <-  Y * log(mu) - mu - lgamma(Y + 1)
  loglike <- -1*loglike
  loglike <- sum(loglike)
  return(loglike)
  
}


# optimize
vP <- c('b0' = 0, 'b1' = 0.5)
  start <- Sys.time()
out <- optim(par = vP, fn = Poisson_loglike, method = 'BFGS', 
             dat = dat, Y = 'Y_pois')
  end <- Sys.time()

  end - start
# Compare Generating Parameter to R package estimate and the custom code estimate
cbind('**Gen Param**' = c(as.vector(Beta)), '**Custom**' = out$par, '**R package**' = c(coef(mod.pois)))
