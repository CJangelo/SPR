
#
#
#
# Custom written estimation routine:
Neg_Binom_loglike_Beta <- function(vP, dat, Y, theta){
  

  Beta <- vP[c('b0', 'b1')]
  Beta <- matrix(Beta, ncol = 1)
  size <- theta

  X <- model.matrix( ~ Group , data = dat) 
  # Parameters:
  mu <- exp(X %*% Beta)
  prob <- size/(size + mu)
  # Stratify the probabilities on 2 groups (would be latent groups in LCM)
  prob <- aggregate(prob ~ Group, FUN = mean, data = dat)
    prob <- prob$V1
  out <- aggregate(as.formula(paste0(Y, ' ~ Group')), FUN = mean, data = dat)
    Y <- out[ , Y]

  loglike <- lgamma(Y + size) - lgamma(size) - lgamma(Y + 1) + size*log(prob) + (Y)*log(1-prob)
  loglike <- -1*loglike
  loglike <- sum(loglike)

  return(loglike)
  
}


# Custom written estimation routine:
Neg_Binom_loglike_size <- function(vP, dat, Y, Beta){
  

  Beta <- matrix(Beta, ncol = 1)
  size <- vP['theta']

  X <- model.matrix( ~ Group , data = dat) 
  # Parameters:
  mu <- exp(X %*% Beta)
  prob <- size/(size + mu)
  # Stratify the probabilities on 2 groups (would be latent groups in LCM)
  prob <- aggregate(prob ~ Group, FUN = mean, data = dat)
    prob <- prob$V1
  out <- aggregate(as.formula(paste0(Y, ' ~ Group')), FUN = mean, data = dat)
    Y <- out[ , Y]

  check <- vector()
for(size in 1:100){
  loglike <- lgamma(Y + size) - lgamma(size) - lgamma(Y + 1) + size*log(prob) + (Y)*log(1-prob)
  loglike <- -1*loglike
  loglike <- sum(loglike)
  check <- rbind(check, c(size, loglike))
}
  check
  plot(check)
  return(loglike)
  
}


# optimize
Y = 'Y_nb'
vP1 <- c('b0' = 0.2, 'b1' = 1)
vP2 <- c('theta' = 1)
delta <- 0

while(diff < 0.01){
  
  est0 <- c(vP1, vP2)

  out_beta <- optim(par = vP1, fn = Neg_Binom_loglike_Beta, method = 'BFGS', 
             dat = dat, Y = 'Y_nb', theta = vP2)

  vP1 <- out_beta$par

  out_theta <- optim(par = vP2, fn = Neg_Binom_loglike_size, method = 'BFGS', 
             dat = dat, Y = 'Y_nb', Beta = vP1)
  
  vP2 <- out_theta$par
  
  delta <- max(abs(c(vP1, vP2) - est0)))
       