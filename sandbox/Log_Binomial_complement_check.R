
# Simulation study to see what we're dealing with here:
# TODO: include covariates to simulate the adjusted models
# Switch to Poisson model rather than log binomial


rm(list = ls())
gc()

library(SPR)
library(MASS)

N = 1e3 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1

out <- vector()

for (repl in 1:1000){

  set.seed(9072021 + repl)


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
p <- exp(XB)
range(p)


# Dichotomize
dat$Y_binom1 <- 1*(p > runif(n = N))
dat$Y_binom2 <- 1 - dat$Y_binom1

# -----------------------
mod1 <- glm(Y_binom1 ~ 1, data = dat, start = -1, family = binomial(link ='log'))
mod2 <- glm(Y_binom2 ~ 1, data = dat, start = -1, family = binomial(link ='log'))

ci1 <- suppressMessages( confint(mod1) )
ci2 <- suppressMessages( confint(mod2) )

# Are the estimates complements?
check1 <- round(exp(coef(mod1)), 3) == round(1-exp(coef(mod2)), 3)

# Are the conf inter complements?
check2 <- all(round(exp(ci1), 3) == round(rev(1-exp(ci2)), 3))

# p-value null is P = 0, not P = 1
pval <- summary(mod2)$coef[,'Pr(>|z|)']
# if zero is NOT in the confidence interval:
pval2 <- all(!(0 >= 1-exp(ci2)))
check3 <- (pval < 0.05 & pval2 == TRUE)

out <- rbind(out, c(check1, check2, check3))


cat(paste0('replication ', repl, '\n'))

}

# check the checks:
apply(out, 2, table)

