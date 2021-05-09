

# Beta-distribution Regression
# Works, g2g
# Trying different probability distributions
# Seems promising!
# 2.19.21
# TODO: clean this shit up so it's a matrix
# Then go on to generate proportions this way:
# https://cran.r-project.org/web/packages/ciTools/vignettes/ciTools-binomial-vignette.html
# try glm() approach WITH DIFFERENT SIZE GROUPS
# Compare to binom.test, prop.test, and Fisher's exact test:
# https://stats.stackexchange.com/questions/123609/exact-two-sample-proportions-binomial-test-in-r-and-some-strange-p-values



rm(list = ls())
gc()

library(MASS)

N = 1e5 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
#set.seed(2182021)



    dat <- data.frame(
                      'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                     # 'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                     'Group' = c(rep('Group_1', 0.75*N), rep('Group_2', 0.25*N)),
                      'Y_comp' = rep(NA, N*number.timepoints), 
                      #'Bio' = rep(rnorm(N, mean = 0, sd = 1), number.timepoints),
                      'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
                      stringsAsFactors=F)
    
        
        # Design Matrix
    X <- model.matrix( ~ Group , data = dat) 
    param <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    param[] <- c(3, 7)

    
    # Parameters:
    mu <- X %*% param
    dat$XB <- as.vector(mu)
    dat$Y <- rbeta(n = N, shape1 = mu, shape2 = 5)
    unique(mu/(mu + 5))

    
    # check
    hist(dat$Y, xlim = c(0,1))
    hist(dat$Y[dat$Group == 'Group_1'], xlim = c(0,1))
    hist(dat$Y[dat$Group == 'Group_2'], xlim = c(0,1))
    aggregate(Y ~ 1, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = mean, data = dat)
    aggregate(Y ~ Group, FUN = var, data = dat)
    


vP <- c('shape1_0' = 3, 'shape1_1' = 7, 'shape2' = 5)

Beta_loglike <- function(vP, dat){
  
  
  Y1 <- dat$Y[dat$Group == 'Group_1']
  Y2 <- dat$Y[dat$Group == 'Group_2']
  p1 <- sort(Y1)
  p2 <- sort(Y2)
  Fn1 <- ecdf(p1)
  Fn2 <- ecdf(p2)
  x <- seq(from = 0, to = 1, length.out = 1000)
  p1 <- Fn1(x)
  p2 <- Fn2(x)
  p <- cbind(p1, p2)
  
  
  # the 'par' (parameters to estimate) are only passed as a vector
  param.hat <- vP[c('shape1_0', 'shape1_1')]
  param.hat <- matrix(param.hat, ncol = 1)
  #shape2 <- vP['shape2']

  X <- model.matrix( ~ Group , data = dat) 
  dat$mu <- X %*% param.hat
  out <- aggregate(mu ~ Group, FUN = function(x) 'mean' = mean(x), data = dat)
  mu <- out$V1
  
  
  xx <- matrix(x, nrow = length(x), ncol = 2)
  mu <- matrix(mu, nrow = length(x), ncol = 2, byrow = T)

  # PDF:
  tmp1 <- xx^(mu-1) 
  tmp2 <- (1-xx)^(vP['shape2'] - 1)
  tmp3 <- beta(mu, vP['shape2'])
  q <- (tmp1*tmp2)/tmp3
  #q <- 1/(sigma * sqrt(2*pi)) * exp( -0.5*((x - mu)/sigma)^2 )
  qq <- apply(q, 2, function(x) cumsum(x)/sum(x))

  # plot(qq[, 1])
  #   plot(p[, 1])
  # plot(qq[, 2])
  #   plot(p[, 2])

  DD <- sum((p - qq)^2, na.rm = F)

  return(DD)
  

}

# optimize
#vP <- c('b0' = 0.2, 'b1' = 0, 'sigma' = 1)
out <- optim(par = vP, fn = Beta_loglike, method = 'BFGS', dat = dat)
out$par
param
vP
