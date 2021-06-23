########################################################
#
#
#           GENERATE
#           Ordered Data
#           Cross sectional
#
#
############################################################
#

rm(list = ls())
gc()

library(MASS)
library(polycor)
#source("C:/Users/ciaconangelo/Documents/RESEARCH/R_CODE_Long_Mixed_Models/Function_findRhoBin.R")

N = 5000 #this should be divisible by however many groups you use!
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
Beta[] <- c(0.0, 1)

# Matrix multiply:
XB <- X %*% Beta


# Define thresholds:
thresholds <- c(0.1, 0.4, 0.7, 0.9) # probabilities of the normal distribution
p <- c(0, thresholds, 1)
diff(p) # here's the proportion in each category

mu <- mean(XB) # mean of the latent variable
var.theta <- 1 # variance of the latent variable
sigma2 <- var(XB) + var.theta # total error of the distribution
zeta <- qnorm(thresholds, mean = mu, sd = sqrt(sigma2))
  zeta <- matrix(zeta, nrow = N, ncol = length(thresholds), byrow = T)

# Theta is the latent variable
theta <- rnorm(n = N, mean = XB, sd = sqrt(var.theta))
  theta <- matrix(theta, nrow = N, ncol = ncol(zeta), byrow = F)

dat$Y_ord <- apply(theta > zeta, 1, sum)
polyserial(theta[ , 1], dat$Y_ord)

aggregate(Y_ord ~ Group, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
barplot(100*table(dat$Y_ord)/sum(table(dat$Y_ord)), ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Ordinal')


# Fit Models:
mod <- MASS:::polr(as.factor(Y_ord) ~ Group , data= dat, method = 'probit', Hess = T)
summary(mod)
mod$coefficients # can't estimate an intercept here
Beta
mod$zeta
zeta[1,]

# Logistic, not probit
# Note that the variance of the logistic disribution is pi/3 rather than 1
mod.logit <- MASS:::polr(as.factor(Y_ord) ~ Group , data= dat, method = 'logistic', Hess = T)
summary(mod.logit)
mod.logit$coefficients
mod.logit$zeta
