
# This is extension to part 1
# include unbalanced covariate
# adjusted model gets you back to your correct value
# adjusted model aligns with pt 1, where there's no covariate influence
# unadjusted model yields misleading result



#-------------
rm(list = ls())
gc()

library(SPR)
library(MASS)

N = 1e5 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(9242021)

dat <- data.frame('USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)

dat$Covariate <- rlnorm(n = N, mean = dat$Group == 'Group_2', sd = 0.25)
range(dat$Covariate)
aggregate(Covariate  ~ Group,
          FUN = function(x) round(mean(x, na.rm = T),2),
          dat = dat, na.action = na.pass)
# Covariate is unbalanced across groups, could serve as a confound

# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group + Covariate , data = dat)

# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(-0.25, -1, -0.2)

# Matrix multiply:
XB <- X %*% Beta
p <- exp(XB) # prevalence - log link
round(range(p), 2)
#p <- exp(XB)/(1 + exp(XB)) - this is the logistic link

# Dichotomize
dat$Y_binom <- 1*(p > runif(n = N))

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
# Percentage achieving "1" across Groups:
ptab

# Compare to model output:
mod <- glm(Y_binom ~ Group, data = dat, family = binomial(link ='log'))
summary(mod)
exp(coef(mod)['(Intercept)' ]) # Group 1
exp(sum(coef(mod))) # Group 2
# equal to observed - won't be equal if adjusted (i.e., other variables in model)


# What is the relative risk?
ptab['1', 'Group_2']/ptab['1', 'Group_1']
# Group 1 is almost 4 times as likely to have event as Group 2

# model estimate of relative risk:
exp(coef(mod)['GroupGroup_2'])
# aligned with observed - remember this only lines up exactly because
# you have no other variables in the model. If you include other variables,
# then it's the "adjusted relative risk" and it won't line up exactly with
# the observed values. This is the same as the other examples.

# Adjusted Relative Risk - log binomial model with covariate (confound)
mod2 <- glm(Y_binom ~ Group + Covariate,
            start = c('(Intercept)' = 0, 'GroupGroup_2' = -1, 'Covariate' = -1),
            data = dat, family = binomial(link ='log'))
summary(mod2)
Beta

# Observed:
ptab # Group 1 & Group 2
ptab['1', 'Group_2']/ptab['1', 'Group_1'] # RR

# Unadjusted model:
exp(coef(mod)['(Intercept)' ]) # Group 1
exp(sum(coef(mod))) # Group 2
exp(coef(mod)['GroupGroup_2']) # RR

# Adjusted model
exp(coef(mod2)['GroupGroup_2']) # RR

#--------------------------------
# Logistic Regression:
mod.lr <- glm(Y_binom ~ Group + Covariate, data = dat, family = binomial(link ='logit'))
summary(mod.lr)


# Zhu and Yang, 1998 Method: https://jamanetwork.com/journals/jama/fullarticle/188182
OR <- exp(coef(mod.lr)['GroupGroup_2'])
P0 <- ptab['1', 'Group_1'] # requires the proportion of control subjects who experience the outcome.
denom <- (1-P0) + (P0 * OR)
RR <- OR/denom
RR

# Compute relative risk using logit link model
p1 <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_2',
                                  'Covariate' = mean(dat$Covariate)))

p0 <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_1',
                                   'Covariate' = mean(dat$Covariate)))

# adjusted RR from logit link model:
p1/p0

# this is pretty close but not exactly the value of the RR
# from the adjusted log-binomial model:
exp(coef(mod2)['GroupGroup_2'])
# TO DO: WHY DOES THIS NOT LINE UP EXACTLY!?!?!
RR


# 9.9.21
# probability of Y given Group membership, averaged over Covariate
p1i <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_2',
                                  'Covariate' = dat$Covariate))

p0i <- predict(mod.lr, type = 'response',
              newdata = data.frame('Group' = 'Group_1',
                                   'Covariate' = dat$Covariate))

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
#mean(p1i/p0i)

# Compare all of this to the code in "Log_Binomial_pt1.R"
# Look at the results of the model in pt1
# there is no covariate, unadjusted model output:
# > exp(coef(mod)['(Intercept)' ]) # Group 1
# (Intercept)
#     0.60488
# > exp(sum(coef(mod))) # Group 2
# [1] 0.22462
# > exp(coef(mod)['GroupGroup_2'])
# GroupGroup_2
#    0.3713464
# >
# The point is that your adjusted log-binomial model estimates what the
# prevalence would be if there was no covariate effect
