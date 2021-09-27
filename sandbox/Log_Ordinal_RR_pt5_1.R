
# Extension of Log Binomial
# Another illustration of relative risk
# This time computing it in ordinal data

# include a covariate
# when reduced down to two categories, it aligns perfectly



#-------------
rm(list = ls())
gc()
N = 1e5 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(9152021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)

# Unbalanced:
dat$Z <- rbinom(n = N, size = 1, prob = 0.4 + 0.4*(dat$Group == 'Group_1'))
table(dat$Z)
xtabs(~ Group + Z, data = dat)

# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group + Z , data = dat)


# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(0, -0.25, -0.25)

# Matrix multiply:
XB <- X %*% Beta

# Thresholds:
thr <- c(-2, -1, 0)
eta <- matrix(thr, nrow = nrow(XB), ncol = length(thr), byrow = T) -
  matrix(XB, nrow = nrow(XB), ncol = length(thr), byrow = F)

p <- exp(eta)/(1 + exp(eta))
#p <- exp(eta)
dat$p <- p
dat$Y <- apply(runif(n = N)  > p, 1, sum)


# Quick barplot, base R plot functions
aggregate(Y  ~ Group, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
xtabs( ~ Y + Group, dat = dat, na.action = na.pass)
barplot(100*table(dat$Y)/sum(table(dat$Y)),
        ylim = c(0, 100), ylab = 'Percentage',
        col = 'grey', main = 'Ordinal')

library(ggplot2)
ggplot(data = dat, aes(x= as.factor(Y), fill = Group)) +
      geom_bar(aes(y = (..count..)/sum(..count..)), color="#e9ecef", alpha=0.6, position="dodge", stat="count") +
      scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
      scale_fill_manual(name = 'Groups', values=c("blue2", "red2")) +
      theme_minimal() +
      theme(legend.position = 'bottom') +
      labs(x = 'Score', y = 'Percentage',
           title = 'Scores with a Binomial Distribution',
           caption = 'Note: Simulated data')

#----
# Key here is that this data is ORDERED

# Observed Data:
xtabs(~ Y + Group, data = dat)
ptab <- prop.table(xtabs(~ Y + Group, data = dat), 2)
# Percentage achieving each category across Groups:
ptab

# Fit ordinal regression
# use the predicted probabilities to compute relative risk
library(ordinal)
mod.clm <- clm(as.factor(Y) ~ Group + Z, data = dat)
summary(mod.clm)
beta <- summary(mod.clm)$beta
alpha <- summary(mod.clm)$alpha


# Group 2:
X1 <- cbind('Group' = rep(1, N), 'Z' = dat$Z)
BB <- matrix(summary(mod.clm)$beta)
XB <- X1 %*% BB
#XB <- matrix(summary(mod.clm)$beta)
thr <- summary(mod.clm)$alpha
eta <- matrix(thr, nrow = nrow(XB), ncol = length(thr), byrow = T) -
  matrix(XB, nrow = nrow(XB), ncol = length(thr), byrow = F)
p1 <- exp(eta)/(1 + exp(eta))

# Group 1
X0 <- cbind('Group' = rep(0, N), 'Z' = dat$Z)
BB <- matrix(summary(mod.clm)$beta)
XB <- X0 %*% BB
#XB <- matrix(0) # obviously this is more complicated with more predictors!
thr <- summary(mod.clm)$alpha
eta <- matrix(thr, nrow = nrow(XB), ncol = length(thr), byrow = T) -
  matrix(XB, nrow = nrow(XB), ncol = length(thr), byrow = F)
p0 <- exp(eta)/(1 + exp(eta))

p0 <- cbind(0, p0, 1)
p1 <- cbind(0, p1, 1)
p0 <- t(apply(p0, 1, diff))
p1 <- t(apply(p1, 1, diff))
p0 <- colMeans(p0)
p1 <- colMeans(p1)
ptab
cbind(p0, p1)
ptab[, 'Group_1']/ptab[, 'Group_2']
p0/p1

#------------------------------------
# > ptab[, 'Group_1']/ptab[, 'Group_2']
#         0         1         2
# 0.8732260 0.8552308 1.1510231

# > p0/p1
# [1] 0.6966564 0.8638131 1.2358708
