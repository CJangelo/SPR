
# Use with "Ordinal_probit_v1"

# Ordinal_logit_v1 and Ordinal_logit_v2
# Some differences in the parameterization
# v1 aligns with polr estimates
# v2 flips sign of the intercepts 




# Custom written estimation routine:
ord_loglike <- function(vP, X, Y, k){
  
  N <- nrow(X)
  Y <- Y + 1 # Must start at 1!
  XX <- X[ , -1, drop = F] # drop intercept
  XB <- XX %*% vP[1] 
  XB <- XB %*% matrix(1, nrow = ncol(XB), ncol = k)
  XB <- XB + matrix(vP[-1], nrow = nrow(XB), ncol = k, byrow = T)
  p <- exp(XB)/(1 + exp(XB))
  p <- cbind(1, p, 0)

  like <- p[cbind(1:N, Y)] - p[cbind(1:N, Y + 1)]  
  loglike <- log(like)
  loglike <- -1*sum(loglike, na.rm = F)
  return(loglike)
  
}


# optimize
# Set starting values:
Y <- dat$Y_ord + 1
intercepts.hat = matrix(NA, nrow = 1, ncol = length(unique(Y)) - 1)
  for(score in 2:length(unique(Y))){
intercepts.hat[ , score -1] = log(mean(Y >= score, na.rm =T))
  }

vP <- c(1, intercepts.hat)
out <- optim(par = vP, fn = ord_loglike, method = 'BFGS', 
             X = X, Y = dat$Y_ord,  k = length(thresholds))
# Compare Generating Parameter to R package estimate and the custom code estimate
summary(mod.logit)
mod.logit$coefficients
out$par
