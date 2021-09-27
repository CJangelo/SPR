

# Count Data
# Poisson & Zero-Inflated Poisson
#

rm(list = ls())
gc()

library(MASS)

N = 1e4 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(7062021)

dat <- data.frame(
  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
  'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
  stringsAsFactors=F)


# Design Matrix
X <- model.matrix( ~ Group , data = dat)
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(0.2, 1)

# Parameters:
mu <- exp(X %*% Beta)
dat$mu <- as.vector(mu)

size <- theta <- 1 # dispersion parameter
prob <- theta/(theta + mu)


# Zero-inflation ADDED:
zi <- 1
P_zi <- 1/(1 + exp(-zi)) #scalar
Y <- 0:100

#--------------------------
#
  Y_zip <- rep(NA, N)
  Y_pois <- rep(NA, N)
  i <- 1
for (i in 1:nrow(prob)){
  #prob.i <- prob[i,]
  mu.i <- mu[i,]

  # -----  Description:
  # If response = 0
  # Can be 0 two different ways, either inflation OR from the count process
  # 'OR' means you add the probabilities (then take the log)
    prob0 <- log(P_zi + (1 - P_zi)*exp(-mu.i))
  # prob0 <- log(P_zi + (1-P_zi)*prob.i^size)


  # If response > 0
  # multiply probability that it's NOT inflated zero times count process
  # This is an "AND" statement, requires multiplying probabilities

  # Poisson, Not Negative Binomial
  prob1 <- log(1 - P_zi) +  Y * log(mu.i) - mu.i - lgamma(Y + 1)

  # prob1 <- log(1 - P_zi) +
  #   lgamma(Y + size) -
  #   lgamma(size) -
  #   lgamma(Y + 1) +
  #   size*log(prob.i) +
  #   (Y)*log(1-prob.i)

  logP <- (Y == 0) * prob0  +   (Y != 0) *  prob1

  # The probabilities need to sum to 1:
  P <- exp(logP)
  P <- P/sum(P)

  Y_zip[i] <- sample(x = Y, size = 1, prob = P)
  Y_pois[i] <- rpois(n = 1, lambda = mu.i)

}


# Check Data Generation:
table(Y_zip == 0)
table(Y_pois ==0)
hist(Y_zip)
hist(Y_pois)
aggregate(Y_zip ~ Group, FUN = mean, data = dat)
aggregate(mu ~ Group, FUN = mean, data = dat)
aggregate(Y_pois ~ Group, FUN = mean, data = dat)

#--------
# Fit Models
library(glmmTMB)
# Poisson
mod.pois <- glmmTMB(Y_pois ~ Group,
                    ziformula = ~ 0,
                    data = dat,
                    family = poisson)
summary(mod.pois)


# Zero-Inflated Poisson
mod.zip <- glmmTMB(Y_zip ~ Group,
                    ziformula = ~ 1,
                    data = dat,
                    family = poisson)
summary(mod.zip)
