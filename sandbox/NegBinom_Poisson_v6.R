


# Custom estimation code
# Estimate Negative Binomial
# Attempted dimension reduction to mimic EM algorithm
# Appears to work...
# Try implementing this in LCM/DCM 


# Custom written estimation routine:
Neg_Binom_loglike_v2 <- function(vP, dat, Y){
  
  # the 'par' (parameters to estimate) are only passed as a vector
  Beta <- vP[c('b0', 'b1')]
  Beta <- matrix(Beta, ncol = 1)
  size <- vP['theta']

  X <- model.matrix( ~ Group , data = dat) 
  # Parameters:
  muu <- exp(X %*% Beta)
  prob <- size/(size + muu)
  dat <- cbind.data.frame(dat, prob)
  # Stratify the probabilities on 2 groups (would be latent groups in LCM)
  out <- aggregate(prob ~ Group, FUN = mean, data = dat)
  # okay so you have a probability for each of the two groups - there that would be latent groups
  x <- unique(dat$Y_nb)
  x <- sort(x)
  
  prob <- matrix(out$prob, nrow = length(x), ncol = length(out$prob), byrow = T)
  x <- matrix(x, nrow = length(x), ncol = ncol(prob), byrow = F)
  
  # Negative Binomial distribution (PMF) 
  ## https://stat.ethz.ch/R-manual/R-devel/library/stats/html/NegBinomial.html
  #Px <- gamma(x + size)/(gamma(size) * factorial(x)) * prob^size * (1 - prob)^x
  logPx <- lgamma(x + size) - lgamma(size) - lgamma(x + 1) + size*log(prob) + (x)*log(1-prob)
 ## Px is the likelihood, and logPx is the loglikelihood

  # Now create the "Expected Values" or pseudo-counts:
  dat <- cbind.data.frame(dat, factor(dat[ , Y]))
  ev <- xtabs(~ Y_nb + Group, data = dat)
  ev <- as.matrix(ev)

  # Use same loglikelihood as you're accustomed to:
  #loglike <- ev * log(Px) # if computed Px above
    loglike <- ev * logPx   # if computed logPx above
  loglike <- -1*loglike
  loglike <- sum(loglike)

  return(loglike)
  
  
}


# optimize
Beta; theta # Gen param
Y = 'Y_nb'
vP <- c('b0' = 0.5, 'b1' = .5, 'theta' = .5)

out_v2 <- optim(par = vP, fn = Neg_Binom_loglike_v2, method = 'BFGS', 
             dat = dat, Y = 'Y_nb')

round(out_v2$par, 2)
# Compare Generating Parameter to R package estimate and the custom code estimate
cbind('Gen Param' = c(as.vector(Beta), theta), 'Custom' = out_v2$par, 'R package' = c(coef(mod.nb), mod.nb$theta))

