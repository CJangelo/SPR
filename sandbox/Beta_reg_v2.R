


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
    param <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    param[] <- c(3, 7)


    # Parameters:
    shape2 <- 5
    shape1.i <- X %*% param
    dat$XB <- as.vector(shape1.i)
    dat$Y <- stats::rbeta(n = N, shape1 = shape1.i, shape2 = shape2)
    unique(shape1.i/(shape1.i + shape2))


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

# alternative parameterization:
Beta.hat <- coef(mod)[c('(Intercept)', 'GroupGroup_2')]
eta.hat <- c(Beta.hat['(Intercept)'], sum(Beta.hat))
  names(eta.hat) <- c('Group_1', 'Group_2')
mu.hat <- exp(eta.hat)/(1 + exp(eta.hat)) # model-based mu
phi.hat <- coef(mod)['(phi)']

# Parameter recovery:
mu.hat # model-based estimate of the mean
unique(shape1.i/(shape1.i + shape2)) # generating parameters
aggregate(Y ~ Group, FUN = mean, data = dat) # observed

mu.hat*(1 - mu.hat)/(1 + phi.hat) # model-based estimate of the variance
unique( (shape1.i*shape2)/ ((shape1.i + shape2 + 1)*(shape1.i + shape2)^2) ) # gen param
aggregate(Y ~ Group, FUN = var, data = dat) # observed



#----
Beta_loglike <- function(vP, dat){

  # the 'par' (parameters to estimate) are only passed as a vector
  BB <- vP[c('b0', 'b1')]
  BB <- matrix(BB, ncol = 1)
  phi <- vP['phi']

  X <- model.matrix( ~ Group , data = dat)
  mu <- X %*% BB
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
cbind(out$par, SE)
summary(mod)
