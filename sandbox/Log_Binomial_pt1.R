
# 9.1.21: put this together due to COURAGE study
# review relative risk with team
# here is how to model prevalence data - relative risk - find paper
# This is just modified version of the Binomial data vignette
# instead of OR, you are using the relative risk and percentages
# This is more intuitive for audience
# include in the SPR training sessions

# Section 2 of the Jstat soft paper is a nice walk thru:
# https://www.jstatsoft.org/article/view/v086i09


# Zhu and Yang, 1998 Method: https://jamanetwork.com/journals/jama/fullarticle/188182


#-------------
rm(list = ls())
gc()

library(SPR)
library(MASS)

N = 1e5 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(9062021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group  , data = dat)

# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(-0.5, -1)

# Matrix multiply:
XB <- X %*% Beta
p <- exp(XB) # prevalence - log link
#p <- exp(XB)/(1 + exp(XB)) - this is the logistic link
round(range(p), 2)

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
# Group 1 is almost 3 times as likely to have event as Group 2

# model estimate of relative risk:
exp(coef(mod)['GroupGroup_2'])
# aligned with observed - remember this only lines up exactly because
# you have no other variables in the model. If you include other variables,
# then it's the "adjusted relative risk" and it won't line up exactly with
# the observed values. This is the same as the other examples.


#--------------------------------
# You can get the relative risk from the logistic regression:

# Logistic Regression:
mod.lr <- glm(Y_binom ~ Group, data = dat, family = binomial(link ='logit'))
summary(mod.lr)

# Zhu and Yang, 1998 Method: https://jamanetwork.com/journals/jama/fullarticle/188182
# Should align with Frank Harrell papers
OR <- exp(coef(mod.lr)['GroupGroup_2'])
P0 <- ptab['1', 'Group_1'] # requires the proportion of control subjects who experience the outcome.
denom <- (1-P0) + (P0 * OR)
RR <- OR/denom
RR

# One way is just to output the probabilities:
p1 <- predict(mod.lr, type = 'response', newdata = data.frame('Group' = c('Group_2')))
p0 <- predict(mod.lr, type = 'response', newdata = data.frame('Group' = c('Group_1')))
p1/p0

# Another way is based on the links
# https://cran.r-project.org/web/packages/logisticRR/vignettes/logisticRR.html
# also slide 17 here: http://www.omicron.dk/relative-risk.pdf
# Hold all covariates equal, then compute. Here it is easy, no covariates

# break it out into the different terms so it's easier to read:
t1 <- exp(coef(mod.lr)['GroupGroup_2'])
t2 <- 1 + exp(coef(mod.lr)['(Intercept)'])
t3 <- 1 + exp(sum(coef(mod.lr)))
(t1*t2)/t3
# aligns with observed and model estimate
