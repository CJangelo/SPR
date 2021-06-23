########################################################
#
#
#           GENERATE
#           Binomial Data
#           Cross sectional
#
#
############################################################
#

rm(list = ls())
gc()

library(ggplot2)
#library(MASS)
#library(polycor)
#source("C:/Users/ciaconangelo/Documents/RESEARCH/R_CODE_Long_Mixed_Models/Function_findRhoBin.R")

N = 5000 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(2092021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group, data = dat)

# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(-0.5, 2)

# Matrix multiply:
XB <- X %*% Beta
p <- exp(XB)/(1 + exp(XB))
hist(p, xlim = c(0, 1))

# Dichotomize
dat$Y_binom <- 1*(p > runif(n = N))

# Plot
aggregate(Y_binom  ~ Group, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
barplot(100*table(dat$Y_binom )/sum(table(dat$Y_binom )), ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Binomial')

# Added 5.27.21
library(ggplot2)
ggplot(data = dat, aes(x= as.factor(Y_binom), fill = Group)) +
      geom_bar(aes(y = (..count..)/sum(..count..)), color="#e9ecef", alpha=0.6, position="dodge", stat="count") +
      scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
      scale_fill_manual(name = 'Groups', values=c("blue2", "red2")) +
      theme_minimal() +
      theme(legend.position = 'bottom') +
      labs(x = 'Score', y = 'Percentage',
           title = 'Scores with a Binomial Distribution',
           caption = 'Note: Simulated data')

# Estimate
mod <- glm(Y_binom ~ Group, data = dat, family = 'binomial')
summary(mod)
exp(1.99821)
xtabs(~ Y_binom + Group, data = dat)
G2 <- 2062/438
G1 <- 974/1526
G2/G1

# Probability that Group 1 is sick:
exp(-0.448)/(1 + exp(-0.448))
# Probability that Group 2 is sick:
exp(-0.448 + 1.99821)/(1 + exp(-0.448 + 1.99821))
# Compare to observed proportions:
prop.table(xtabs(~ Y_binom + Group, data = dat), 2)




# Custom rolled glm code:
loglike <- function(vP, X, Y){

  Y <- cbind(1-Y, Y)
  XB <- X %*% vP
  p <- 1/(1 + exp(-1*(XB)))
  P <- cbind(1-p, p)
  # Loglikelihood:
  loglike <- Y*log(P)
  loglike <- -1*loglike # optimization will minimize function
  loglike <- sum(loglike)
  return(loglike)

}

# optimize
vP <- rep(0, length(Beta))
out <- optim(par = vP, fn = loglike, hessian = T,  method = 'BFGS', X = X, Y = dat$Y_binom)
# Compare Generating Parameter to my estimate and the glm() estimate
cbind(Beta, round(out$par,4), round(coef(mod),4))
# Compare the standard errors:
sqrt(diag(solve(out$hessian)))


# Custom rolled glm code:
loglike_v2 <- function(vP, X, Y){

  XB <- X %*% vP
  p <- 1/(1 + exp(-1*(XB)))
  # Stratify the probabilities on 2 groups (would be latent groups in LCM)
  p <- aggregate(p ~ dat$Group, FUN = mean)
  p <- p$V1
  # Create probability of each observed response (0/1):
  P <- cbind(1-p, p)

  # Stratify the average observed response on 2 groups (would be latent groups in LCM)
  Y <- aggregate(Y ~ dat$Group, FUN = mean)
  Y <- Y$param
  # Create average for each observed response (0/1):
  Y <- cbind(1-Y, Y)

  # Loglikelihood:
  loglike <- Y*log(P)
  loglike <- -1*loglike # optimization will minimize function
  loglike <- sum(loglike)
  return(loglike)

}

# optimize
vP <- rep(0, length(Beta))
#Y = dat$Y_binom
#vP <- c(-0.5, 2)
out_v2 <- optim(par = vP, fn = loglike_v2, method = 'BFGS', X = X, Y = dat$Y_binom)
round(out_v2$par, 4)
# Compare Generating Parameter to my estimate and the glm() estimate
cbind('Generating param' = Beta, 'R package' =  round(coef(mod),4),
      'Custom v1' = round(out$par,4), 'Custom v2' = round(out_v2$par,4))

