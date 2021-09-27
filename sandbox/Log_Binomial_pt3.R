
# Custom Estimation
# Note: X %*% Beta must be less than or equal to zero
# see pg 6 of the https://www.jstatsoft.org/article/view/v086i09
# unclear how to implement using constrained BFGS, because the constraints
# have to apply to parameter being estimated (Beta) not X %*% Beta
#
#

loglike <- function(vP, X, Y){

  Y <- cbind(1-Y, Y)
  XB <- X %*% vP
  p <- exp(XB) # log link
  #p <- 1/(1 + exp(-1*(XB))) # logistic link
  P <- cbind(1-p, p)
  # Loglikelihood:
  loglike <- Y*log(P)
  loglike <- -1*loglike # optimization will minimize function
  loglike <- sum(loglike)
  return(loglike)

}


vP <- rep(-1, length(Beta))
out <- optim(par = vP, fn = loglike, hessian = T,  method = 'BFGS', X = X, Y = dat$Y_binom)
# Compare Generating Parameter to my estimate and the glm() estimate
cbind(Beta, round(out$par,4), round(coef(mod),4))
# Compare the standard errors:
sqrt(diag(solve(out$hessian)))
