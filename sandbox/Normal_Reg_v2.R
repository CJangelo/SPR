

# Normal-distribution Regression
# Works, g2g
# This COULD be interesting - try with different probability density functions
# 2.19.21


rm(list = ls())
gc()

library(MASS)

N = 1e4 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
#set.seed(2182021)



    dat <- data.frame(
                      'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                      'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                      'Y_comp' = rep(NA, N*number.timepoints), 
                      #'Bio' = rep(rnorm(N, mean = 0, sd = 1), number.timepoints),
                      'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
                      stringsAsFactors=F)
    
        
        # Design Matrix
    X <- model.matrix( ~ Group , data = dat) 
    Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    Beta[] <- c(0.2, 1)
    sigma <- 1
    
    
    # Parameters:
    XB <- X %*% Beta
    dat$XB <- as.vector(XB)
    error <- rnorm(n = N, mean = 0, sd = sigma)
    dat$Y <- dat$XB + error
    
    # Log normal:
    #error <- rlnorm(n = N, meanlog = 0, sd = 1)
    #dat$Y <- dat$XB + error

    
    # check
    hist(dat$Y)
    aggregate(Y ~ 1, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = var, data = dat)
    
    # Model?
    mod <-lm(Y ~ Group, data = dat)
    #summary(mod)
    cbind(Beta, coef(mod))
    # #####
  

vP <- c('b0' = 0.2, 'b1' = 1, 'sigma' = 1)

Normal_loglike <- function(vP, dat){
  
  
  # the 'par' (parameters to estimate) are only passed as a vector
  Beta.hat <- vP[c('b0', 'b1')]
  Beta.hat <- matrix(Beta.hat, ncol = 1)
  sigma <- vP['sigma']

  X <- model.matrix( ~ Group , data = dat) 
  dat$mu <- X %*% Beta.hat
  out <- aggregate(mu ~ Group, FUN = function(x) 'mean' = mean(x), data = dat)
  mu <- out$V1
  Y1 <- dat$Y[dat$Group == 'Group_1']
  Y2 <- dat$Y[dat$Group == 'Group_2']
  p1 <- sort(Y1)
  p2 <- sort(Y2)
  Fn1 <- ecdf(p1)
  Fn2 <- ecdf(p2)
  x <- seq(from = min(dat$Y), to = max(dat$Y), length.out = 1000)
  p1 <- Fn1(x)
  p2 <- Fn2(x)
  p <- cbind(p1, p2)
  
  xx <- matrix(x, nrow = length(x), ncol = 2)
  mu <- matrix(mu, nrow = length(x), ncol = 2, byrow = T)

  q <- 1/(sigma * sqrt(2*pi)) * exp( -0.5*((x - mu)/sigma)^2 )
    #q <- as.vector(q)
  qq <- apply(q, 2, function(x) cumsum(x)/sum(x))
  #qq <- cumsum(q)/sum(q)
  
  # plot(qq[, 1])
  #   plot(p[, 1])
  # plot(qq[, 2])
  #   plot(p[, 2])

  DD <- sum((p - qq)^2, na.rm = F)
  ##DD <- sum(p * log(p/qq), na.rm = T)
  #KL <- sum(p * (log(p) - log(qq)), na.rm = T)
  #KL <- -1*KL
  #print(vP)
  #KL <- sum(q * (log(qq) - log(p)), na.rm = T)
  #KL <- -1*sum(p*log(qq))
  return(DD)
  

}

# optimize
#vP <- c('b0' = 0.2, 'b1' = 0, 'sigma' = 1)
out <- optim(par = vP, fn = Normal_loglike, method = 'BFGS', dat = dat)
out$par
cbind(Beta, coef(mod) )
cbind(
  'Gen param' = c(Beta, sigma),
  'Custom' = c(out$par), # N-p
  'R package' = c(coef(mod), summary(mod)$sigma^2)
)