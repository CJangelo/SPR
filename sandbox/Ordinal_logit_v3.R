

rm(list = ls())
gc()

library(MASS)
library(polycor)
#source("C:/Users/ciaconangelo/Documents/RESEARCH/R_CODE_Long_Mixed_Models/Function_findRhoBin.R")

N = 100 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(6242021)

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

# Thresholds:
thr <- c(0, 1, 2, 3)
eta <- matrix(thr, nrow = nrow(XB), ncol = length(thr), byrow = T) -
  matrix(XB, nrow = nrow(XB), ncol = length(thr), byrow = F)

p <- exp(eta)/(1 + exp(eta))
dat$p <- p
dat$Y <- apply(runif(n = N)  > p, 1, sum)


write.csv(dat[, c('USUBJID', 'Group', 'Y')],
            quote = F,
            file = './sandbox/dataset2.csv')

# Quick barplot, base R plot functions
barplot(100*table(dat$Y)/sum(table(dat$Y)),
        ylim = c(0, 100), ylab = 'Percentage',
        col = 'grey', main = 'Ordinal')

# Key here is that this data is ORDERED


# Fit Models:
mod <- MASS::polr(as.factor(Y) ~ Group , data= dat, method = 'logistic', Hess = T)
summary(mod)
mod$coefficients # can't estimate an intercept here
Beta
mod$zeta
thr


library(ordinal)
mod.clm <- clm(as.factor(Y) ~ Group, data = dat)
summary(mod.clm)
mod.clm$beta
# Compute odds ratios by hand
# Examine the proportional odds assumption
# Notice that it's not exact
tab <- addmargins(xtabs(~ Y + Group, data = dat), 1)
#
# >0
G1 <- (tab['Sum', 'Group_1'] - tab['0','Group_1'])/tab['0', 'Group_1']
G2 <- (tab['Sum', 'Group_2'] - tab['0','Group_2'])/tab['0', 'Group_2']
G2/G1
exp(mod.clm$beta)
#
# > 1 (includes 0 and 1 categories)
G1 <- (tab['Sum', 'Group_1'] - sum(tab[c('0', '1'), 'Group_1']))/
  sum(tab[c('0', '1'), 'Group_1'])
G2 <- (tab['Sum', 'Group_2'] - sum(tab[c('0', '1'), 'Group_2']))/
  sum(tab[c('0', '1'), 'Group_2'])
G2/G1
exp(mod.clm$beta)
#
# > 2 (includes 0, 1, 2)
G1 <- (tab['Sum', 'Group_1'] - sum(tab[c('0', '1', '2'), 'Group_1']))/
  sum(tab[c('0', '1', '2'), 'Group_1'])
G2 <- (tab['Sum', 'Group_2'] - sum(tab[c('0', '1', '2'), 'Group_2']))/
  sum(tab[c('0', '1', '2'), 'Group_2'])
G2/G1
exp(mod.clm$beta)
#
# > 3 (includes 0, 1, 2, 3)
G1 <- (tab['Sum', 'Group_1'] - sum(tab[c('0', '1', '2', '3'), 'Group_1']))/
  sum(tab[c('0', '1', '2', '3'), 'Group_1'])
G2 <- (tab['Sum', 'Group_2'] - sum(tab[c('0', '1', '2', '3'), 'Group_2']))/
  sum(tab[c('0', '1', '2', '3'), 'Group_2'])
G2/G1
exp(mod.clm$beta)
#

# Custom written estimation routine:
ord_loglike <- function(vP, X, Y, k){

  N <- nrow(X)
  Y <- Y + 1 # Must start at 1!
  XX <- X[ , -1, drop = F] # drop intercept
  XB <- XX %*% vP[1] # Not sure why polr parameterizes it this way
  XB <- XB %*% matrix(1, nrow = ncol(XB), ncol = k)
  eta <- matrix(vP[-1], nrow = nrow(XB), ncol = k, byrow = T) - XB
  p <- exp(eta)/(1 + exp(eta))
  p <- cbind(0, p, 1)

  like <- p[cbind(1:N, Y + 1)] - p[cbind(1:N, Y)]
  loglike <- log(like)
  loglike <- -1*sum(loglike, na.rm = F)
  return(loglike)

}


# optimize
# Set starting values:
Y <- dat$Y + 1
intercepts.hat = matrix(NA, nrow = 1, ncol = length(unique(Y)) - 1)
  for (score in 2:length(unique(Y))){
    intercepts.hat[ , score -1] = log(mean(Y >= score, na.rm =T))
  }

# vP is the vector of parameters passed to the optim() function
vP <- c(1, -1*intercepts.hat)
# optim() function optimizes the objective function
out <- optim(par = vP,
             hessian = T,
             fn = ord_loglike,
             method = 'BFGS',
             X = X,
             Y = dat$Y,
             k = length(intercepts.hat))

# Compare Generating Parameter to R package estimate and the custom code estimate
cbind('polr() estimate' = c(mod$coefficients, mod$zeta),
      'clm() estimate' = c(mod.clm$beta, mod.clm$alpha),
      'custom code' = out$par)

# Standard errors:
SE <- sqrt(diag(solve(out$hessian)))
zval <- out$par/SE
2*(1-pnorm(tval))
summary(mod.clm)

