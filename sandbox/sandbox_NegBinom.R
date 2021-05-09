
#
#
#
# Custom written estimation routine:
Neg_Binom_loglike_v2 <- function(vP, dat, Y){
  
  #print(vP)
  # the 'par' (parameters to estimate) are only passed as a vector
  # re-create what you need from that vector
  Beta <- vP[c('b0', 'b1')]
  Beta <- matrix(Beta, ncol = 1)
  size <- vP['theta']
  #size <- 1 #works if you set it to 1!
  #print(size)
  
  X <- model.matrix( ~ Group , data = dat) 
  # Parameters:
  mu <- exp(X %*% Beta)
  prob <- size/(size + mu)
  # Stratify the probabilities on 2 groups (would be latent groups in LCM)
  prob <- aggregate(prob ~ Group, FUN = mean, data = dat)
  prob <- prob$V1
  print(prob)
  
  # New:
  # Stratify the average observed response on 2 groups (would be latent groups in LCM)
  out <- aggregate(as.formula(paste0(Y, ' ~ Group')), FUN = mean, data = dat)
  Y <- out[ , Y]
  print(Y)
  #print(prob)
  loglike <- Y * log(prob)
  #loglike <- lgamma(Y + size) - lgamma(size) - lgamma(Y + 1) + size*log(prob) + (Y)*log(1-prob)
  loglike <- -1*loglike
  loglike <- sum(loglike)

  #cat( paste0( vP), '\n')
  return(loglike)
  
}

#size <- 1

# optimize
 Beta; theta # Gen param
Y = 'Y_nb'
vP <- c('b0' = 0.2, 'b1' = 1, 'theta' = 1)
#vP <- c('b0' = 0.2, 'b1' = 1) # Test - this works!! Why is 'theta' breaking?!?!
out_v2 <- optim(par = vP, fn = Neg_Binom_loglike_v2, method = 'BFGS', 
             dat = dat, Y = 'Y_nb')
out_v2$par
# Compare Generating Parameter to nnet package estimate and the custom code estimate

cbind('Gen Param' = c(as.vector(Beta), theta), 
      'R package' = c(coef(mod.nb), mod.nb$theta), 
      'Custom v1' = out$par,
      'Custom v2' = out_v2$par)



