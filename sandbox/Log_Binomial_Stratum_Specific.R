

rm(list = ls())
gc()
library(SPR)
library(MASS)

N = 1e5 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(9122021)

dat <- data.frame('USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints), 'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints), stringsAsFactors=F)

# Generate a Covariate
# Option 3: Categorical, Unbalanced:
dat$Z <- rbinom(n = N, size = 1, prob = 0.4 + 0.4*(dat$Group == 'Group_1'))

# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group + Z , data = dat)

# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(-2, -1, 1.75)

# Matrix multiply:
XB <- X %*% Beta
p <- exp(XB) # prevalence - log link
round(range(p), 2)
#> [1] 0.17 0.78
#p <- exp(XB)/(1 + exp(XB)) - this is the logistic link

# Dichotomize
dat$Y_binom <- as.vector(1*(p > runif(n = N)))

#-----------------------------------------------------------------

# Observed Data:
ptab <- prop.table(xtabs(~ Y_binom + Group, data = dat), 2)

# Compare to model output:
mod <- glm(Y_binom ~ Group, data = dat, family = binomial(link ='log'))

# Adjusted Relative Risk - log binomial model with covariate (confound)
mod2 <- glm(Y_binom ~ Group + Z,
            start = as.vector(Beta),
            data = dat, family = binomial(link ='log'))

#--------------------------------
# Logistic Regression:
mod.lr <- glm(Y_binom ~ Group + as.factor(Z), data = dat,
              family = binomial(link ='logit'))

# Stratum specific RR
# Intercept from the logistic regression model:
# Pulled this from Stata function: http://fmwww.bc.edu/repec/bocode/l/logittorisk.ado
odds0 <- exp( coef(mod.lr)['(Intercept)'] )
  odds0/(1 + odds0)
odds0 <- exp( sum(coef(mod.lr)[c('(Intercept)', 'as.factor(Z)1')]) )
  odds0/(1 + odds0)
# These align with Predicted values below
# Group 1, Z = 0, and Group 1, Z = 1

predict(mod.lr, type = 'response', newdata = data.frame('Group' = 'Group_1','Z' = 0))
predict(mod.lr, type = 'response', newdata = data.frame('Group' = 'Group_2','Z' = 0))

predict(mod.lr, type = 'response', newdata = data.frame('Group' = 'Group_1','Z' = 1))
predict(mod.lr, type = 'response', newdata = data.frame('Group' = 'Group_2','Z' = 1))

# Log Reg Stratum Specific RR:
predict(mod.lr, type = 'response', newdata = data.frame('Group' = 'Group_2','Z' = 1))/
  predict(mod.lr, type = 'response', newdata = data.frame('Group' = 'Group_1','Z' = 1))

predict(mod.lr, type = 'response', newdata = data.frame('Group' = 'Group_2','Z' = 0))/
  predict(mod.lr, type = 'response', newdata = data.frame('Group' = 'Group_1','Z' = 0))


predict(mod2, type = 'response', newdata = data.frame('Group' = 'Group_2','Z' = 0))/
  predict(mod2, type = 'response', newdata = data.frame('Group' = 'Group_1','Z' = 0))

# Non-collapsibility
# Observed Stratum-Specific RR:
ptabZ0 <- prop.table(xtabs(~ Y_binom + Group, data = dat[dat$Z == 0, ]), 2)
ptabZ0
ptabZ0['1', 'Group_2']/ptabZ0['1', 'Group_1']

ptabZ1 <- prop.table(xtabs(~ Y_binom + Group, data = dat[dat$Z == 1, ]), 2)
ptabZ1
ptabZ1['1', 'Group_2']/ptabZ1['1', 'Group_1']

# So stratum-specific RR (i.e. RR within each level of Z)  is equal
# to the adjusted marginal RR?
# This does not line up with the OR estimates at all

p1i <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_2',
                                  'Z' = dat$Z))

p0i <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_1',
                                   'Z' = dat$Z))

mean(p1i)/mean(p0i)
exp(Beta['GroupGroup_2', ])
exp(coef(mod2)['GroupGroup_2'])
