

rm(list = ls())
gc()

library(MASS)
library(polycor)
#source("C:/Users/ciaconangelo/Documents/RESEARCH/R_CODE_Long_Mixed_Models/Function_findRhoBin.R")

N = 5000 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(6152021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group  , data = dat)


# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(0.0, 2)

# Matrix multiply:
XB <- X %*% Beta

# Thresholds:
thr <- c(-2, -1, 0, 1)
eta <- matrix(XB, nrow = nrow(XB), ncol = length(thr), byrow = F) +
  matrix(thr, nrow = nrow(XB), ncol = length(thr), byrow = T)

p <- exp(eta)/(1 + exp(eta))
dat$p <- p
dat$Y <- apply(runif(n = N)  > p, 1, sum)


# Quick barplot, base R plot functions
barplot(100*table(dat$Y)/sum(table(dat$Y)),
        ylim = c(0, 100), ylab = 'Percentage',
        col = 'grey', main = 'Ordinal')


# Fit Models:
mod <- MASS:::polr(as.factor(Y) ~ Group , data= dat, method = 'logistic', Hess = T)
summary(mod)
mod$coefficients # can't estimate an intercept here
Beta
mod$zeta
thr

# Compute odds ratios by hand
# Examine the proportional odds assumption
# Notice that it's not exact
addmargins(xtabs(~ Y + Group, data = dat), 1)
# P <= 0
G1 <- (289)/(2500-289)
G2 <- (1226)/(2500-1226)
exp(2)
G2/G1
# P <= 1
G1 <- (289+378)/(2500-289-378)
G2 <- (1226+574)/(2500-1226-574)
G2/G1
# P <= 2
G1 <- (289+378+574)/(2500-289-378-574)
G2 <- (1226+574+394)/(2500-1226-574-394)
G2/G1
#... you can do the rest...
