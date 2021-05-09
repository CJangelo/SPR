

# Count Data
# Poisson & Negative binomial 
##################################################

# Compare the Poisson and Negative Binomial distributions
# Develop illustration of that


rm(list = ls())
gc()

library(MASS)

N = 1e5 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(2012021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  #'Y_nb' = rep(NA, N*number.timepoints), 
                  #'Y_pois' = rep(NA, N*number.timepoints), 
                  #'Bio' = rep(rnorm(N, mean = 0, sd = 1), number.timepoints),
                  'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
                  stringsAsFactors=F)


# Design Matrix
X <- model.matrix( ~ Group , data = dat) 
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(0.2, 1)

# Parameters:
mu <- exp(X %*% Beta)
dat$mu <- as.vector(mu)

theta <- 1 # dispersion parameter
dat$theta <- as.vector(theta)

upsilon <- mu + (mu^2)/theta
dat$upsilon <- as.vector(upsilon)

prob <- theta/(theta + mu)

dat$Y_nb <- rnbinom(n = N, size = theta, mu = mu)
  hist.nb <- hist(dat$Y_nb, plot = F)
  hist.nb$counts <- 100*hist.nb$counts/N
  plot(hist.nb, ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Negative Binomial')


dat$Y_pois <- rpois(n = N, lambda = mu)
  hist.pois <- hist(dat$Y_pois, plot = F)
  hist.pois$counts <- 100*hist.pois$counts/N
  plot(hist.pois, ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Poisson')

# Poisson:
# Mean
aggregate(Y_pois ~ Group, FUN = mean, data = dat)
aggregate(mu ~ Group, FUN = mean, data = dat)
# Variance
aggregate(Y_pois ~ Group, FUN = var, data = dat)

# Negative Binomial:
# Mean:
aggregate(Y_nb ~ Group, FUN = mean, data = dat)
aggregate(mu ~ Group, FUN = mean, data = dat)

# Variance:
aggregate(Y_nb ~ Group, FUN = var, data = dat)
aggregate(upsilon ~ Group, FUN = mean, data = dat)


# Fit Models:

mod.pois <- glm(Y_pois ~ Group, data = dat, family = 'poisson')
summary(mod.pois)

mod.nb <- MASS::glm.nb(Y_nb ~ Group, data = dat)
summary(mod.nb)

# Relationship between Negative Binomial and Poisson
# TODO: Illustrate how the Negative Binomial distribution
# 1. Changes shape as mu changes
# 2. Changes shape as theta changes
# 3. When does NB approach Poisson distribution?
# Please show multiple plots together to visualize the change

