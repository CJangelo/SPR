


# https://cran.r-project.org/web/packages/betareg/vignettes/betareg.pdf
# Our interest in what follows will be more closely related to the recent literature, i.e., modeling continous
# random variables that assume values in (0, 1), such as rates, proportions, and concentration
# or inequality indices (e.g., Gini).

# Page 5:
# It is noteworthy that the beta regression model described above was developed to allow
# practitioners to model continuous variates that assume values in the unit interval such as
# rates, proportions, and concentration or inequality indices (e.g., Gini). However, the data
# types that can be modeled using beta regressions also encompass proportions of “successes”
# from a number of trials, if the number of trials is large enough to justify a continuous model. In
# this case, beta regression is similar to a binomial generalized linear model (GLM) but provides
# some more flexibility – in particular when the trials are not independent and the standard
# binomial model might be too strict. In such a situation, the fixed dispersion beta regression
# is similar to the quasi-binomial model (McCullagh and Nelder 1989) but fully parametric.
# Furthermore, it can be naturally extended to variable dispersions.


rm(list = ls())
gc()


N = 300 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(6282021)



    dat <- data.frame(
                      'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                     'Group' = c(rep('Group_1', 0.75*N), rep('Group_2', 0.25*N)),
                      'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
                      stringsAsFactors=F)


    # Design Matrix
    X <- model.matrix( ~ Group , data = dat)
    Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    Beta[] <- c(-.5, 1)


    # Parameters:
    mu <- X %*% Beta
    mu <- as.vector(mu)
    mu <- exp(mu)/(1 + exp(mu))
    phi <- 10

    # Re-parameterize, see pg 3 of PDF
    shape1.i <- mu*phi
    shape2.i <- phi - mu*phi

    # Confirm that parmeterization is equivalent:
    unique(shape1.i/(shape1.i + shape2.i))
    unique(mu)

    # Generate Data:
    dat$Y <- stats::rbeta(n = N, shape1 = shape1.i, shape2 = shape2.i)


    # check
    hist(dat$Y, xlim = c(0,1))
    hist(dat$Y[dat$Group == 'Group_1'], xlim = c(0,1))
    hist(dat$Y[dat$Group == 'Group_2'], xlim = c(0,1))
    aggregate(Y ~ 1, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = var, data = dat)


#---------------
# Fit Model
library(betareg)
mod <- betareg::betareg(Y ~ Group, data = dat)
summary(mod)



# Parameter recovery:
unique(shape1.i/(shape1.i + shape2.i)) # generating parameters
aggregate(Y ~ Group, FUN = mean, data = dat) # observed


#----
Beta_loglike <- function(vP, dat){

  # the 'par' (parameters to estimate) are only passed as a vector
  Beta.hat <- vP[c('b0', 'b1')]
  Beta.hat <- matrix(Beta.hat, ncol = 1)
  phi <- vP['phi']

  X <- model.matrix( ~ Group , data = dat)
  mu <- X %*% Beta.hat
  mu <- as.vector(mu)
  mu <- exp(mu)/(1 + exp(mu))
  Y <- dat$Y

  loglike <- lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi) +
    (mu*phi - 1)*log(Y) +
    ( (1-mu)*phi - 1)*log(1-Y)

  loglike <- -1*sum(loglike, na.rm = T)
  return(loglike)

}

#---------
# optimize
vP <- c('b0' = 0, 'b1' = 1, 'phi' = 2)
out <- optim(par = vP,
             fn = Beta_loglike,
             method = 'BFGS',
             hessian = T,
             dat = dat)
SE <- sqrt(diag(solve(out$hessian)))
cbind('Custom code Beta' = out$par,
      'Custom code SE' = SE,
      'Model Beta' = coef(mod),
      'Model SE' = sqrt(diag(vcov(mod))))
