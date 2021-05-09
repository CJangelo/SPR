

# Count Data
# Poisson & Negative binomial 
#########################################################

# Evaluate Type I error and Power
# To discuss: Power is meaningless in the context of inflated Type I error rates
# 


rm(list = ls())
gc()

library(MASS)

N = 50 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(2012021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
                  stringsAsFactors=F)


# Design Matrix
X <- model.matrix( ~ Group , data = dat) 
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
#Beta[] <- c(0.2, 0)  # Type I error
Beta[] <- c(0.2, 1)  # Power

# Parameters:
mu <- exp(X %*% Beta)
dat$mu <- as.vector(mu)

theta <- 0.5 # dispersion parameter
dat$theta <- as.vector(theta)

upsilon <- mu + (mu^2)/theta
dat$upsilon <- as.vector(upsilon)


#######
# Simulation
out <- vector()

for(repl in 1:1000){
  
  # Generate Data:
    dat$Y_nb <- rnbinom(n = N, size = theta, mu = mu)
    #dat$Y_pois <- rpois(n = N, lambda = mu)
 

  # Fit Models -both Poisson and NB to data that is NB
    mod.pois <- glm(Y_nb ~ Group, data = dat, family = 'poisson')
    mod.nb <- MASS::glm.nb(Y_nb ~ Group, data = dat)

    out <- rbind(out, c( 
                 summary(mod.pois)$coef['GroupGroup_2', 'Pr(>|z|)'], 
                 summary(mod.nb)$coef['GroupGroup_2', 'Pr(>|z|)']
    ))
  
  cat(paste0('Replication: ', repl, '\n'))
  
  }

colMeans(out < 0.05)
