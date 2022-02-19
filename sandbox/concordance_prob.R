

# Concordance Proability
# Ordinal Data
# 1 if greater, 0 if less, 0.5 if a tie
# next steps
# 1. make it work for longitudinal data, CLMM and the GEE

# 2. logistic regression and
# 3. OLS regression

rm(list = ls())
gc()

library(SPR)
library(ordinal)


N = 2000 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(02032022)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group  , data = dat)


# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
Beta[] <- c(0.0, 1)

# Matrix multiply:
XB <- X %*% Beta

# Thresholds:
thr <- c(0, 1, 2, 3)
eta <- matrix(thr, nrow = nrow(XB), ncol = length(thr), byrow = T) -
  matrix(XB, nrow = nrow(XB), ncol = length(thr), byrow = F)

p <- exp(eta)/(1 + exp(eta))
dat$p <- p


#----------------------
#
  out <- vector()

for( repl in 1:1000){

  dat$Y <- apply(runif(n = N)  > p, 1, sum)
  mod.clm <- ordinal::clm(as.factor(Y) ~ Group, data = dat)

  # Concordance Probability:
  true.beta <- Beta[2,]
  C_true <- 1/(1 + exp(-1*true.beta/1.52))

  beta.hat <- mod.clm$beta
  C_hat <- 1/(1 + exp(-1*beta.hat/1.52))

  g1 <- dat[dat$Group == 'Group_1', 'Y' ]
  g2 <- dat[dat$Group == 'Group_2', 'Y' ]

  g1.re <- sample(x = g1, size = 100, replace = T)
  g2.re <- sample(x = g2, size = 100, replace = T)

  C_half <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
  C_half <- mean(C_half)


  #
  num <- sum(g2.re > g1.re)
  denom <- sum(g2.re > g1.re) + sum(g1.re > g2.re)
  C_no_ties <- num/denom

  #
  C_emp1 <- mean(sample(x = g2, size = 100, replace = T) >
                   sample(x = g1, size = 100, replace = T))

  #
   C_emp2 <- mean(sample(x = g2, size = 100, replace = T) >=
                   sample(x = g1, size = 100, replace = T))


   out <- rbind(out, c(C_true, C_hat, C_emp1, C_emp2, C_no_ties, C_half))
   cat('Replication: ', repl, '\n')

}

apply(out, 2, mean)
