

# Count Data
# Poisson & Negative binomial 
##################################################

# Compare the Poisson and Negative Binomial distributions
# Relationship between Negative Binomial and Poisson
# Illustrate how the Negative Binomial distribution
# 1. Changes shape as mu changes
# 2. Changes shape as theta changes
# 3. When does NB approach Poisson distribution?
# Please show multiple plots together to visualize the change


rm(list = ls())
gc()

library(MASS)
library(ggplot2)
library(gridExtra)

N = 1e5 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(4202021)

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



for (theta in c(1000, 100, 10, 1, 0.1)) {
  
  #theta <- 1000# dispersion parameter
  dat$theta <- as.vector(theta)

  upsilon <- mu + (mu^2)/theta
  dat$upsilon <- as.vector(upsilon)

  prob <- theta/(theta + mu)

 # Generate Data:
  dat$Y_nb <- rnbinom(n = N, size = theta, mu = mu)
  dat$Y_pois <- rpois(n = N, lambda = mu)

  p1 <- ggplot2::ggplot(dat, aes(x = Y_nb)) + 
            geom_bar(aes(y = 100*(..count..)/sum(..count..))) +
            labs(subtitle = paste0('Theta = ', theta), title='Negative Binomial', x ="Value", y = "Percent") +
            ylim(c(0, 25)) + 
            theme_minimal() + 
            theme(plot.title = element_text(hjust = 0.5), 
                  plot.subtitle = element_text(hjust = 0.5)) 
            # Why the fuck would I want a left aligned title

  
  p2 <- ggplot2::ggplot(dat, aes(x = Y_pois)) + 
            geom_bar(aes(y = 100*(..count..)/sum(..count..))) +
            labs(subtitle = paste0('Theta = ', theta), title='Poisson', x ="Value", y = "Percent") +
            ylim(c(0, 50)) + 
            theme_minimal() + 
            theme(plot.title = element_text(hjust = 0.5), 
                  plot.subtitle = element_text(hjust = 0.5)) 
            # Why the fuck would I want a left aligned title

  
  gridExtra::grid.arrange(p1, p2, ncol = 2)   


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

}

