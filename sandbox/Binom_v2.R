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
# 5.27.21:  added example showing that with additional predictors this
# can no longer aligned because it's a conditional estimate

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
                  'Covariate' = rnorm(n = N, mean = 0, sd = 1),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group + Covariate , data = dat)

# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(-0.5, 2, 1)

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
      scale_fill_manual(name = 'Latent Classes', values=c("blue2", "red2")) +
      theme_minimal() +
      theme(legend.position = 'bottom') +
      labs(x = 'Score', y = 'Percentage',
           title = 'Scores with a Binomial Distribution',
           caption = 'Note: Simulated data')

# Estimate
mod <- glm(Y_binom ~ Group, data = dat, family = 'binomial')
summary(mod)
exp(1.70170)
xtabs(~ Y_binom + Group, data = dat)
G2 <- 1963/537
G1 <- 1000/1500
G2/G1

# Estimate the Group parameter adjusting for covariate
mod2 <- glm(Y_binom ~ Group + Covariate, data = dat, family = 'binomial')
summary(mod2)
exp(2.01047)
# Cannot recreate this conditional/adjusted odds ratio for Groups using
# hand computations.
