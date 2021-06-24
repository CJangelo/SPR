



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
out <- optim(par = vP,
             fn = multinom_loglike,
             method = 'BFGS',
             hessian = T,
             X = X,
             Y = dat$Y_nom_factor,
             data = dat,
             k = k)
# Compare Generating Parameter to nnet package estimate and the custom code estimate
cbind.data.frame('Gen_param' = t(Beta),
                 'nnet R package' = coef(mod),
                 'Custom code' = t(out$par))

# Compare SE:
SE <- sqrt(diag(solve(out$hessian)))
cbind.data.frame('nnet R package SE' = summary(mod)$standard.errors,
                 'Custom code SE' = matrix(SE, ncol = 2, byrow = T))
