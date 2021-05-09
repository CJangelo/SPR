########################################################
#
#
#           GENERATE 
#           Multinomial Model
#           Cross sectional
#
#
############################################################
#

rm(list = ls())
gc()

library(MASS)
library(polycor)
#source("C:/Users/ciaconangelo/Documents/RESEARCH/R_CODE_Long_Mixed_Models/Function_findRhoBin.R")

N = 1e4 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(02032021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group  , data = dat) 

k <- 5 # Number of categories in the nominal item

# Create Beta
Beta <- matrix(0, nrow = ncol(X), ncol = k - 1, dimnames=list(colnames(X), paste0('param', 1:(k-1))))
Beta[1, ] <- c(0.2, 0.8, 0.4, 0.6)
Beta[2, ] <- c(0.2, -1.0, -0.6, 0.8)

# Matrix multiply:
XB <- X %*% Beta

sum.expXB <- apply(exp(XB), 1, sum)
p <- exp(XB)/(1 + sum.expXB)
param0 <-  1 - rowSums(p)
p <- cbind(param0, p)

  out <- vector()
for(i in 1:nrow(p)){
  out <- c(out, 
           sample(x = c(0, 1, 2, 3, 4), size = 1, prob = p[i, ])  
  )
  
} #end loop
 prop.table(table(out))
 colMeans(p)
  
 dat$Y_nom <- out

aggregate(Y_nom ~ Group, FUN = function(x) table(x), dat = dat, na.action = na.pass)
xtabs( ~ Y_nom + Group, data = dat)
prop.table(xtabs( ~ Y_nom + Group, data = dat))
addmargins(prop.table(xtabs( ~ Group + Y_nom, data = dat), 1), 2)
barplot(100*table(dat$Y_nom)/sum(table(dat$Y_nom)), ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Nominal')


# Fit Models:
library(nnet)
dat$Y_nom_factor <- as.factor(dat$Y_nom)
mod <- nnet:::multinom(Y_nom_factor ~ Group, data = dat)
summary(mod)
coef(mod)
t(Beta)

# TODO: Try other R packages for fitting multinomial models
mod0 <- nnet:::multinom(Y_nom ~ 1, data = dat)
mod1 <- nnet:::multinom(Y_nom ~ Group, data = dat)
anova(mod0, mod1)
