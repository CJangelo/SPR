
# Log_Binomial_pt2 - Only with binary coviarate
# Mimic the papers




#-------------
rm(list = ls())
gc()

library(SPR)
library(MASS)

N = 1e5 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(9122021)

dat <- data.frame('USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)

# Balanced:
# dat$Z <- sample(x = c(0,1), size = N, replace = T, prob = c(0.5, 0.5))
# Unbalanced:
dat$Z <- rbinom(n = N, size = 1, prob = 0.4 + 0.4*(dat$Group == 'Group_1'))
table(dat$Z)
xtabs(~ Group + Z, data = dat)

# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group + Z , data = dat)

# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(-0.25, -1, -0.5)

# Matrix multiply:
XB <- X %*% Beta
p <- exp(XB) # prevalence - log link
round(range(p), 2)
#p <- exp(XB)/(1 + exp(XB)) - this is the logistic link

# Dichotomize
dat$Y_binom <- as.vector(1*(p > runif(n = N)))

# Plot
aggregate(Y_binom  ~ Group, FUN = function(x) round(mean(x, na.rm = T),2), dat = dat, na.action = na.pass)
barplot(100*table(dat$Y_binom )/sum(table(dat$Y_binom )), ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Binomial')


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

# Observed Data:
xtabs(~ Y_binom + Group, data = dat)
ptab <- prop.table(xtabs(~ Y_binom + Group, data = dat), 2)

# Compare to model output:
mod <- glm(Y_binom ~ Group, data = dat, family = binomial(link ='log'))

# Adjusted Relative Risk - log binomial model with covariate (confound)
mod2 <- glm(Y_binom ~ Group + Z,
            start = as.vector(Beta),
            data = dat, family = binomial(link ='log'))

#--------------------------------
# Logistic Regression:
mod.lr <- glm(Y_binom ~ Group + Z, data = dat,
              family = binomial(link ='logit'))

# Now we have the odds ratio from the model, we can compute
# the relative risk using the baseline prevalence.

#
# Can think of several ways of getting the baseline prevalance:
#
# 1. observed rate in control group:
p0_v1 <- ptab['1', 'Group_1']
# 2. Intercept from the logistic regression model:
# Pulled this from Stata function:
# http://fmwww.bc.edu/repec/bocode/l/logittorisk.ado
odds0 <- exp( coef(mod.lr)['(Intercept)'] )
p0_v2 <- odds0/(1 + odds0)
# but this assumes the covariates are at reference value:
# 3. Prediction aligns with intercept:
p0_v3 <- predict(mod.lr, type = 'response',
                 newdata = data.frame('Group' = 'Group_1',
                                      'Z' = 0))

# 4. Average predictions over covariates:
p0i <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_1',
                                   'Z' = dat$Z))

p0_v4 <- mean(p0i)

#----
# Compare

# Zhu and Yang, 1998 Method: https://jamanetwork.com/journals/jama/fullarticle/188182
# Harrell et all method as well:
OR <- exp(coef(mod.lr)['GroupGroup_2'])
P0 <- p0_v2
RR <- OR/( (1-P0) + (P0 * OR) )
RR
exp(Beta)
# LOOK at what you're doing here
# You are using OR that ignore Z
# So you are implicitly computing OR of Group assuming
# that Z = 0, reference level
# That's why these values align with
# the p1 and p0 below, where Z=0 in both
# Suddenly it makes more sense to use the standardized
# approach!
# Note that this is the same as the marginal vs conditional
# not aligning with longitudinal binomial data
# it aligns with continuous -
# pg 17 of http://eprints.lse.ac.uk/84163/1/On%20group%20comparisons_Final.pdf

# Compute relative risk using logit link model
p1 <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_2',
                                    'Z' = 0))

p0 <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_1',
                                   'Z' = 0))

# adjusted RR from logit link model:
p1/p0
# Note that this is the same as the v2 and v3 above,
# where you assume the covariate Z is at reference level
# for baseline risk
# OR <- exp(coef(mod.lr)[c('GroupGroup_2')])
# P0 <- p0
# OR/( (1-P0) + (P0 * OR) )
#
# OR <- exp(sum(coef(mod.lr)[c('GroupGroup_2')]))
# P0 <- p1
# OR/( (1-P0) + (P0 * OR) )
# 0.4* 0.431069 + 0.6 * 0.2383973


# 9.9.21
# probability of Y given Group membership, averaged over Covariate
p1i <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_2',
                                  'Z' = dat$Z))

p0i <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_1',
                                   'Z' = dat$Z))

# "Standardized measures are constructed by taking averages over C
# before comparisons (e.g., ratios or differences) across X
#...Recalling Jensen’s inequality (an average of a nonlinear
# function does not equal the function applied to the averages),
# it should not be surprising to find
# divergences between collapsibility conditions depending on
# the step at which averaging is done
# (Samuels, 1981, sec. 3)."
# https://escholarship.org/content/qt3ng2r0sm/qt3ng2r0sm.pdf
mean(p1i)/mean(p0i)
# This seems do-able because you've only got 1 covariate to sum over
# otherwise you'd be taking a lot of averages
# I guess the other option is just to add all the observed
# combinations of the covariates
# but your p1 and p0 will depend on those values

# epitools R package
# Function is 'probratio()'
# See documentation - marginalizes over observed covariates
library(epitools)
pr <- probratio(mod.lr, method='delta', scale='linear')
pr

# Observed:
ptab # Group 1 & Group 2
ptab['1', 'Group_2']/ptab['1', 'Group_1'] # RR

# Unadjusted log model:
exp(coef(mod)['(Intercept)' ]) # Group 1
exp(sum(coef(mod))) # Group 2
exp(coef(mod)['GroupGroup_2']) # RR

# Adjusted log model
exp(coef(mod2)['GroupGroup_2']) # RR

# Logistic Regression (also adjusted)
RR
p1/p0
mean(p1i)/mean(p0i)

# Generating parameter:
exp(Beta)

# bias
exp(coef(mod2)['GroupGroup_2']) - exp(-1)
(mean(p1i)/mean(p0i)) - exp(-1)
(ptab['1', 'Group_2']/ptab['1', 'Group_1']) - exp(-1)
RR - exp(-1)
