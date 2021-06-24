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
           sample(x = c('A', 'B', 'C', 'D', 'E'), size = 1, prob = p[i, ])
  )

} #end loop

 dat$Y_nom <- out

#----
# Compare observed proportions to the probabilities:
 prop.table(table(out))
 colMeans(p)


# Plot data:
 barplot(100*table(dat$Y_nom)/sum(table(dat$Y_nom)), ylim = c(0, 100), ylab = 'Percentage', col = 'grey', main = 'Nominal')


# Fit Models:
library(nnet)
dat$Y_nom_factor <- as.factor(dat$Y_nom)
mod <- nnet:::multinom(Y_nom_factor ~ Group, data = dat)
summary(mod)
coef(mod)
t(Beta)

# What are the estimates?
# Compute the odds ratios

# Cross tabs
tab <- addmargins(xtabs(~ Y_nom + Group, data = dat), 1)
coef(mod)
# Odds of selecting B over reference category A
G1 <- tab['B', 'Group_1'] /tab['A', 'Group_1']
G2 <- tab['B', 'Group_2'] /tab['A', 'Group_2']
G2/G1
exp(coef(mod)['B','GroupGroup_2'])

# Odds of selecting C over reference category A
G1 <- tab['C', 'Group_1'] /tab['A', 'Group_1']
G2 <- tab['C', 'Group_2'] /tab['A', 'Group_2']
G2/G1
exp(coef(mod)['C','GroupGroup_2'])

# Odds of selecting D over reference category A
G1 <- tab['D', 'Group_1'] /tab['A', 'Group_1']
G2 <- tab['D', 'Group_2'] /tab['A', 'Group_2']
G2/G1
exp(coef(mod)['D','GroupGroup_2'])

# Odds of selecting C over reference category A
G1 <- tab['E', 'Group_1'] /tab['A', 'Group_1']
G2 <- tab['E', 'Group_2'] /tab['A', 'Group_2']
G2/G1
exp(coef(mod)['E','GroupGroup_2'])

