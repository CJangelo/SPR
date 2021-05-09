    

# 3.17.21
# g2g
# Resolve the dmvnorm `Sigma` details - why does it work without the correct variance?
# estimation of Beta hat works regardless of the sigma value passed to dmvnorm, same as cross-sectional (see Normal_Reg_v4.R)
# Maybe test with drop-out?



# TODO 
# 1. use the "Function_Generate_Long_Data.R" code, no need for "Gen_MMRM_custom_est_function.R"
# 3. Fuckery with dat vs dat.est, some uses one some uses the other
# also how often do we need to compute some of these things?
# 2. make more flexible - pass the Y you want estimated, not just "dat$Y"
  # Requires "Y" to be the name of the variable you are fitting!
  # pass the following: 
  # USUBJID
  # Y
  # Time
  # Group

# -------------------------------------------------------------------------------------


custom_MMRM_estimation <- function(mod.formula, 
                                   subject_id = 'USUBJID', 
                                   time_var = 'Time',  
                                   var_matrix = 'unstructured', 
                                   dat, 
                                   converg.crit = 1e-4){
  
  #library(mvtnorm)
  # Data prep:
  # mod.formula must be a formula!
  if (class(mod.formula) != 'formula') { stop('`mod.formula` must be a formula') }
    form.split <- sub(" ~.* ", " ", mod.formula)
    reg <- as.formula(paste0('~ ', form.split[3]))
    colnames(dat)[ colnames(dat) == form.split[2] ] <- 'Y'
  # requires "USUBJID" and "Time" variables in the dataframe!
  colnames(dat)[ colnames(dat) == subject_id ] <- 'USUBJID'
  colnames(dat)[ colnames(dat) == time_var ] <- 'Time'

  if ( !all(c('USUBJID', 'Time') %in% colnames(dat)) ) { stop('Subject ID variable needs to be labeled `USUBJID` and your discrete timepoins variable needs to be labeled `Time`') }
    
  N = length(unique(dat$USUBJID))
  number.timepoints <- length(unique(dat$Time))
  # pass the following: 
  # USUBJID
  # Y - or some kinda dependent variable
  # Time
  # Group - or some kinda independent variable
  


#-------------------------------------------------------------  
Normal_loglike <- function(vP, dat, sigma.hat){
  
  Beta.hat <- as.matrix(vP)
  X <- model.matrix( reg , data = dat)  
  Y.hat <- X %*% Beta.hat 
  Y <- matrix(dat$Y, ncol = 1)
  ee <- Y - Y.hat
  dat.est <- cbind.data.frame(dat.est, ee)
    ee.wide <- vector()
  for (tt in unique(dat$Time)) {
    ee.wide <- cbind(ee.wide, dat.est$ee[ dat.est$Time == tt ])
  }
  
  ## TEST:
  if (var_matrix == 'unstructured') {
    
    # LL <- 0
    # for (tt in unique(dat$Time)) {
    #     tmp1 <- ee.wide[, 1:tt, drop = F]
    #     tmp2 <- sigma.hat[1:tt, 1:tt, drop = F]
    #     Py <- mvtnorm::dmvnorm(x = tmp1, mean = rep(0, ncol(tmp1)) , sigma = tmp2)
    #     LL <- LL + -1*sum(log(Py), na.rm = T) 
    # }


    # dat.est$ee[ dat.est$Time == tt ]
    # sigma.hat <- diag(c(1, 1.5, 2), 3, 3)
    # sigma.hat[1, 2:3] <- c(0.8, 0.64)
    # sigma.hat[2, c(1, 3)] <- c(0.8, 0.8)
    # sigma.hat[3, c(1, 2)] <- c(0.64, 0.8)
    # tmp1 <- ee.wide
    # tmp1[is.na(tmp1)] <- 0
    # Py <- mvtnorm::dmvnorm(x = tmp1, mean = rep(0, ncol(ee.wide)), sigma = sigma.hat)
    # LL <- -1*sum(log(Py), na.rm = T)
    
    # Slowest and closest to replicating MMRM:
      LL <- 0
    for (ii in 1:N) {
       tmp1 <- ee.wide[ii, , drop = F]
       nit <- sum(!is.na(tmp1))
       tmp1 <- tmp1[, 1:nit, drop = F]
       tmp2 <- sigma.hat[1:nit, 1:nit, drop = F]
       Py <- mvtnorm::dmvnorm(x = tmp1, mean = rep(0, nit) , sigma = tmp2)
       # LL <- LL + -1*sum(log(Py), na.rm = T)  # T or F, doesn't matter
       LL <- LL + -1*mean(log(Py), na.rm = T)  # T or F, doesn't matter
      }
 
    # # Try this:
    #   LL <- 0
    # for (tt in 1:number.timepoints) {
    #   tmp1 <- ee.wide[ , 1:tt, drop = F]
    #   tmp2 <- sigma.hat[ 1:tt, 1:tt, drop = F]
    #   Py <- mvtnorm::dmvnorm(x = tmp1, mean = rep(0, tt) , sigma = tmp2)
    #   # LL <- LL + -1*sum(log(Py), na.rm = T)
    #   LL <- LL + -1*mean(log(Py), na.rm = T)
    # }#end tt loop
    #
    
    # Py <- mvtnorm::dmvnorm(x = ee.wide, mean = rep(0, ncol(ee.wide)), sigma = sigma.hat)
    #     LL <- -1*sum(log(Py), na.rm = T)

    
      
  }# end unstructured 
    
  if (var_matrix == 'independent') {
    Py <- mvtnorm::dmvnorm(x = ee.wide, mean = rep(0, ncol(ee.wide)), sigma = diag(1, nrow = number.timepoints, ncol = number.timepoints))
    LL <- -1*sum(log(Py), na.rm = T)
   }
    
  #plot(x = Py1, y = Py2)
  
  #LL <- -1*sum(log(Py), na.rm = T)
  return(LL)

}
# -----------------------------------------------------




# -----------------------------------------------------------------
# Initialize
dat.est <- dat
XX <- model.matrix( reg , data = dat)  
Beta.hat <- matrix(0, nrow = ncol(XX), dimnames=list(colnames(XX), 'param'))
SE <- Beta.hat
  SE[] <- 0
vP <- as.vector(Beta.hat)
  names(vP) <- row.names(Beta.hat)
# sigma.hat <- diag(1, nrow = number.timepoints, ncol = number.timepoints)
sigma.hat <- matrix(0.5, nrow = number.timepoints, ncol = number.timepoints)
diag(sigma.hat) <- 1
delta <- 1
it <- 0

# -------------------------------------------------------------------------
# Estimation 
while(delta > converg.crit){

      est0 <- vP
  # M -step:
    out <- optim(par = vP, 
                 fn = Normal_loglike, 
                 method = 'BFGS',
                 hessian = T,
                 dat = dat, 
                 sigma.hat = sigma.hat)
    
    ## Standard Errors:
    if (det(out$hessian) > 0) {
      SE <- sqrt(diag(solve(out$hessian)))
        SE <- matrix(SE, ncol = 1, dimnames = list(names(out$par), 'SE'))
    }
    
    ## Parameters: 
    Beta.hat <- matrix(out$par, ncol = 1, dimnames = list(names(out$par), 'param.hat'))
    X <- model.matrix( reg , data = dat)  
    Y.hat <- X %*% Beta.hat 
    Y <- matrix(dat$Y, ncol = 1)
    
    ## Residuals - Used for Variance-Covariance Matrix Estimation:
    ee <- Y - Y.hat
    dat.est <- cbind.data.frame(dat.est, 'ee' = as.vector(ee)) 
      ee.wide <- vector()
    for(tt in unique(dat$Time)){
      ee.wide <- cbind(ee.wide, dat.est$ee[ dat.est$Time == tt ])
      }

  # E-step update:
      ## Parameters:
    vP <- as.vector(Beta.hat)
      names(vP) <- row.names(Beta.hat)
      ## Variance-Covariance 
    # sigma.hat <- var(ee.wide, use = 'pairwise.complete') # unclear if this is appropriate!
    sigma.hat <- var(ee.wide, use = 'complete.obs') # unclear if this is appropriate!
    # Try others! 'complete.obs' --> then missing values are handled by casewise deletion (and if there are no complete cases, that gives an error).
    
  ## Convergence check:
    delta <- max(abs(est0 - vP)) # only beta parameters drive convergence - TODO
    it <- it + 1
    
} #end while loop 

#----------------------------------------------------------------------------------------

    out <- list('Beta.hat' = Beta.hat, 
                'SE' = SE, 
               'sigma.hat' = sigma.hat, 
               'iterations' = it)

    return(out)

} #end function custom_MMRM_estimation() 



