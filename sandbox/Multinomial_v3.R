



# Custom written estimation routine:
multinom_loglike <- function(vP, X, Y, data, k){
  
  vP <- matrix(vP, nrow = ncol(X), ncol = k - 1)
  Y <- model.matrix( ~ - 1 + Y, data = data)
  XB <- X %*% vP
  sum.expXB <- apply(exp(XB), 1, sum)
  p <- exp(XB)/(1 + sum.expXB)
  param0 <-  1 - rowSums(p)
  p <- cbind(param0, p)

  # Loglikelihood:
  loglike <- Y*log(p)
  loglike <- -1*loglike # optimization will minimize function
  loglike <- sum(loglike)
  return(loglike)
  
}

# optimize
vP <- Beta
out <- optim(par = vP, fn = multinom_loglike, method = 'BFGS', 
             X = X, Y = dat$Y_nom_factor, data = dat, k = k)
# Compare Generating Parameter to nnet package estimate and the custom code estimate
t(Beta)
coef(mod)
t(out$par)
