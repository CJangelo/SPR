

rm(list = ls())
gc()

N = 300 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(7222021)

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
Beta[] <- c(-0.2, 0.0)

# Parameters:
mu <- exp(X %*% Beta)
dat$mu <- as.vector(mu)

theta <- 0.5 # dispersion parameter
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
#aggregate(mu ~ Group, FUN = mean, data = dat)
# Variance
aggregate(Y_pois ~ Group, FUN = var, data = dat)

# Negative Binomial:
# Mean:
aggregate(Y_nb ~ Group, FUN = mean, data = dat)
#aggregate(mu ~ Group, FUN = mean, data = dat)

# Variance:
aggregate(Y_nb ~ Group, FUN = var, data = dat)
#aggregate(upsilon ~ Group, FUN = mean, data = dat)


# Fit Models:

# Poisson
mod1 <- glm(Y_pois ~ Group, data = dat, family = 'poisson')
mod2 <- MASS::glm.nb(Y_pois ~ Group, data = dat)
summary(mod1)
summary(mod2)


# NB:
mod3 <- glm(Y_nb ~ Group, data = dat, family = 'poisson')
mod4 <- MASS::glm.nb(Y_nb ~ Group, data = dat)
summary(mod3)
summary(mod4)

dat <- dat[, c('USUBJID', 'Group', 'Y_nb', 'Y_pois')]
colnames(dat) <- c('USUBJID', 'Group', 'Y1', 'Y2')

  write.table(dat, na = '.', quote = F, sep = ', ', col.names = F, row.names = F, file = 'dataset3.txt')
