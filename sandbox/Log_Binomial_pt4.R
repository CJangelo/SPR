
# Simulation study to see what we're dealing with here:


rm(list = ls())
gc()

library(SPR)
library(MASS)

N = 30 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(9012021)

out1 <- vector()
out2 <- vector()

for (repl in 1:1000){

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ 1  , data = dat)

# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(-0.1)

# Matrix multiply:
XB <- X %*% Beta
p <- exp(XB) # prevalence!
#p <- exp(XB)/(1 + exp(XB))


# Dichotomize
dat$Y_binom1 <- 1*(p > runif(n = N))
dat$Y_binom2 <- 1 - dat$Y_binom1


# Plot
# aggregate(Y_binom  ~ Group, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
#
# barplot(100*table(dat$Y_binom1 )/sum(table(dat$Y_binom1 )), ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Binomial')
# barplot(100*table(dat$Y_binom2 )/sum(table(dat$Y_binom2 )), ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Binomial')


# -----------------------
mod1 <- glm(Y_binom1 ~ 1, data = dat, start = -1, family = binomial(link ='log'))
mod2 <- glm(Y_binom2 ~ 1, data = dat, start = -1, family = binomial(link ='log'))

out1 <- rbind(out1, summary(mod1)$coef)
out2 <- rbind(out2, summary(mod2)$coef)

cat(paste0('iteration ', repl, '\n'))

}


#-----------------------------------------------------------
#
#
#

# hist(out1[, 'Estimate'])
# hist(exp(out1[, 'Estimate']))
# hist(out2[, 'Estimate'])
# hist(exp(out2[, 'Estimate']))
# plot(x = out1[, "Pr(>|z|)"], y = out2[, "Pr(>|z|)"])
#
#
# hist(out1[, "Pr(>|z|)"])
# hist(out2[, "Pr(>|z|)"])
# quantile(out1[, "Pr(>|z|)"])
# quantile(out2[, "Pr(>|z|)"])
#
#
# inv.out1 <- c(quantile(out1[, 'Estimate'], 0.025),
#              median(out1[, 'Estimate']),
#              quantile(out1[, 'Estimate'], 0.975))
# exp(inv.out1)
# inv.out2 <- c(quantile(out2[, 'Estimate'], 0.025),
#              median(out2[, 'Estimate']),
#              quantile(out2[, 'Estimate'], 0.975))
# exp(inv.out2)

mean(out1[, "Pr(>|z|)"] < 0.05)
mean(out2[, "Pr(>|z|)"] < 0.05)
# so they're not just the complement

tmp <- cbind(out1[, "Pr(>|z|)"] < 0.05, out2[, "Pr(>|z|)"] < 0.05)
table(apply(tmp, 1, sum))
# with n=30, only 623/1000 were one or the other

tmp <- cbind(tmp, apply(tmp, 1, sum))
round(head(cbind(out1, out2)), 2)

summary(mod1)$coef
summary(mod2)$coef
exp(Beta)
est <- exp(summary(mod1)$coef)[,'Estimate']
margin <- as.numeric(exp(1.96*sqrt(vcov(mod1))))
est - margin
confint(mod1)
ci <- (confint(mod1))
plnorm(q = ci, mean = 0, sd = 1)
2*pnorm(q = summary(mod1)$coef[,'z value'], mean = 0, sd = 1)
2*pnorm(q = summary(mod2)$coef[,'z value'], mean = 0, sd = 1)

# can we match?
# exp(coef(mod1))
# exp(sqrt((vcov((mod1)))))
# exp(confint(mod1))

#2*pnorm(q = exp(-0.9999522), mean = exp(0), sd = exp(1))

#2*plnorm(q = exp(-0.9999522), meanlog = 0, sdlog = 1)
2*plnorm(q = summary(mod1)$coef[,'z value'], meanlog = 0, sdlog = 1)
2*plnorm(q = summary(mod2)$coef[,'z value'], meanlog = 0, sdlog = 1)

#2*pnorm(q = exp(-0.9999522), mean = 0, sd = 1)
#2*pnorm(q = exp(-0.9999522), mean = 0, sd = exp(1))

#mean(out)
#quantile(out, c(0.025, 0.50, 0.975))
# inv.out <- c(quantile(out, 0.025), mean(out), quantile(out, 0.975))
# inv.out
# exp(inv.out)
# 1-exp(Beta)
# > inv.out
#      2.5%               97.5%
# -2.207275 -1.720357 -1.347074
# > exp(inv.out)
#      2.5%               97.5%
# 0.1100000 0.1790022 0.2600000
# > 1-exp(Beta)
#                 param
# (Intercept) 0.1812692
# >
