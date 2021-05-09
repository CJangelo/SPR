########################################################
#
#
#           GENERATE 
#           Ordered Data
#           Cross sectional
##            Evaluate Type I error and Power
#
############################################################
#

rm(list = ls())
gc()

library(MASS)
library(polycor)
#source("C:/Users/ciaconangelo/Documents/RESEARCH/R_CODE_Long_Mixed_Models/Function_findRhoBin.R")

N = 30 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(2012021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group  , data = dat) 


# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(0.0, 0) # Type I error
Beta[] <- c(0.0, 1) # Power

# Matrix multiply:
XB <- X %*% Beta


# Define thresholds: (3 thresholds, 4 categories)
thresholds <- c(0.2, 0.4, 0.6, 0.8) # probabilities of the normal distribution
p <- c(0, thresholds, 1)
diff(p) # here's the proportion in each category

mu <- mean(XB) # mean of the latent variable
var.theta <- 1 # variance of the latent variable
sigma2 <- var(XB) + var.theta # total error of the distribution
zeta <- qnorm(thresholds, mean = mu, sd = sqrt(sigma2))
  zeta <- matrix(zeta, nrow = N, ncol = length(thresholds), byrow = T)

  
##################################
  # Replications:
out <- vector()
# Theta is the latent variable
  for(repl in 1:1000){
    
     # Generate Data:
      theta <- rnorm(n = N, mean = XB, sd = sqrt(var.theta))
      theta <- matrix(theta, nrow = N, ncol = ncol(zeta), byrow = F)
      dat$Y_ord <- apply(theta > zeta, 1, sum)

      # Fit Models:
      mod0 <- MASS:::polr(as.factor(Y_ord) ~ 1, data= dat, method = 'probit', Hess = T)
      mod <- MASS:::polr(as.factor(Y_ord) ~ Group , data= dat, method = 'probit', Hess = T)
      tmp <- anova(mod0, mod)
      out <- c(out, tmp$`Pr(Chi)`[2])
      cat(paste0('Replication: ', repl, '\n'))

  }

mean(out < 0.05)

