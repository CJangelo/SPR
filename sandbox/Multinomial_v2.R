########################################################
#
#
#           GENERATE 
#           Multinomial Model
#           Cross sectional
#
#
############################################################
# Simulation Study

rm(list = ls())
gc()

#library(MASS)
#library(polycor)
library(nnet)

N = 100
number.groups <- 2
number.timepoints <- 1
set.seed(2032021)

dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)

# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group  , data = dat) 

k <- 5 # Number of categories in the nominal item

# Create Beta
Beta <- matrix(0, nrow = ncol(X), ncol = k - 1, dimnames=list(colnames(X), paste0('param', 1:(k-1))))
Beta[1, ] <- c(0.2, 0.8, 0.4, 0.6) # Intercepts
#Beta[2, ] <- 0 # Type I error
Beta[2, ] <- c(0.2, -1.0, -0.6, 0.8)  # Power

# Matrix multiply:
XB <- X %*% Beta

sum.expXB <- apply(exp(XB), 1, sum)
p <- exp(XB)/(1 + sum.expXB)
param0 <-  1 - rowSums(p)
p <- cbind(param0, p)

##

out <- vector()
for(repl in 1:1000){
      
      Y <- vector()
    for(i in 1:nrow(p)){
      Y <- c(Y, sample(x = c(0, 1, 2, 3, 4), size = 1, prob = p[i, ]))
      } #end loop
    
     dat$Y_nom <- Y
      # Fit Models:
      mod0 <- nnet:::multinom(Y_nom ~ 1, data = dat, trace = F)
      mod1 <- nnet:::multinom(Y_nom ~ Group, data = dat, trace = F)
    
      tmp <- anova(mod0, mod1)
      out <- c(out, tmp$`Pr(Chi)`[2])
      
      cat('Replication: ', repl, '\n')
      
} #end loop

mean(out < 0.05)
