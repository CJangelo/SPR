#' Simulate Longitudinal Ordinal Data
#'
#' TODO: Build in check if the correlations are too high! Where is that code?
#' Note that the variances aren't used if passed, 
#' Requires MASS R package
#' @param reg.formula this is a regression formula to pass to the data generation; null is
#' formula(~ Group + Time + Time*Group),
#' @param Beta this is the Beta for the regression equation; numeric matrix of Beta values OR a scalar value != 0
#' that will be the final value of the interaction parameters, default is all(Beta == 0) for Type I error simulations
#' @param thresholds the ordinal generation uses a probit approach, so these are the thresholds
#' @param corr options for correlation structure  c('ind', 'ar1', 'cs'), default is 'ar1'
#' @param cor.value numeric, the first corr in ar1, the corr in cs option
#' @param cond.mcar logical; do you want a conditional MCAR data generation
#' @param Covariate logical; do you want a random simulated covariate?
#'
#' @return returns a dataframe containing the simulated data
#' @export

sim_dat_ord <- function(N = 100 , 
                    number.groups = 2, 
                    number.timepoints = 4, 
                    reg.formula = NULL, 
                    Beta = 0, 
                    thresholds = c(0.25, 0.50, 0.75), 
                    corr = 'ar1', 
                    cor.value = NULL,
                    cond.mcar = F,
                    Covariate = F
){
  
 # checks: 
  # if (!is.null(var.values) & length(var.values != number.timepoints))  stop('variance values not equal to number of timepoints')
  if (is.null(reg.formula) & cond.mcar == F) { reg.formula <-  formula(~ Group + Time + Time*Group) }
  if (!is.null(reg.formula) & cond.mcar == T) { 
    if (all(!grepl('Covariate', paste0(reg.formula)))) stop('cannot pass a regression formula without the covariate in it if you want to generate conditional mcar')
  }

 
#------------------------------------------------------------------------------------------------------------------ 
  
  dat <- data.frame(
    'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
    'id.geepack' = rep(1:N, length.out= N*number.timepoints),
    'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
    'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
    'Y_comp' = rep(NA, N*number.timepoints), 
    stringsAsFactors=F)
  

  # Biomarker:
if (Covariate == T) {
   dat$Covariate <- rnorm(N*number.timepoints, mean = 0, sd= 1) # No differences in Biomarker across Groups
  # Note: This is similar to having randomized Biomarker levels across arms in a RCT
  }
  


# Beta parameter default is zero - permits Type I error simulations
if (cond.mcar == F) {

    
  # ----------------------------------------------------------------------------------------
  # Design Matrix
  X <- model.matrix( reg.formula, data = dat) 
  #-------------------------------------------------------------------------------------
  
  if (length(Beta) == 1) {
    if (Beta == 0) {
    
    Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    
    } else {
     
      # If pass scalar, then that is the value of the final interaction parameter 
      beta.values <- seq(0.25, Beta, length.out = number.timepoints - 1)

      Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
      
      Beta[grepl('Group_2', rownames(Beta)) & grepl('Time', rownames(Beta)), ] <- beta.values

    }
  }
}
  

  

  if (cond.mcar == T) {

  # Generate Biomarker:
   dat$Covariate <- rnorm(N*number.timepoints, mean = 1*(dat$Group == 'Group_2'), sd= 1) # Biomarker differs across Groups
  # Design Matrix - Condition MCAR 
   reg.formula <- formula(~ Group + Time + Covariate + Time*Group + Covariate*Time)
   X <- model.matrix( reg.formula, data = dat) # include Biomarker drop-out!
  # Conditional MCAR - biomarker affects Y 
   beta.values <- seq(0.25, Beta, length.out = number.timepoints - 1)
   Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
   Beta[grepl('Covariate', rownames(Beta)) & grepl('Time', rownames(Beta)), ] <- -0.5*beta.values
   Beta[grepl('Group_2', rownames(Beta)) & grepl('Time', rownames(Beta)), ] <- beta.values
  
  }     
    
  
  # --------------------------------------------------------------------------------------------------
  # Matrix multiply:
  rownames(Beta) <- colnames(X)
  XB <- X %*% Beta
  
  
  # -------------------------------------------------------------------
  # DISTRIBUTION OF RESIDUALS
  #
  if (corr == 'ind') {
    
    cor.mat <- diag(1, nrow = number.timepoints, ncol = number.timepoints) 
    
  }# end independent structure
  
  
  if (corr == 'cs') {
    
    if (is.null(cor.value)) { cor.value <- 0.10 }
    cor.mat <- matrix(cor.value, nrow = number.timepoints, ncol = number.timepoints) 
    diag(cor.mat) <- 1
    
  } # end Compound Symmetry correlation
  
  
  if (corr == 'ar1') {
    
    if (is.null(cor.value)) { cor.value <- 0.4 }
    
    cor.mat <- diag(1, nrow = number.timepoints, ncol = number.timepoints) 
    for (i in 1:number.timepoints) {
      for (j in 1:i) {
        
        cor.mat[i , j] <- cor.value^(i -j) # AR1
        cor.mat[j, i] <- cor.mat[i, j]
        
      }
    }
    
  }# end exponential decay correlations
  
# ---------------------------------------------------------------------
# GENERATE ORDINAL DATA
#  This chunk here is the only part that is distinct from "sim_dat.R"

#thresholds <- c(0.25, 0.5, 0.75) # probabilities of the normal distribution
zeta <- qnorm(thresholds, mean = 0, sd =1)
  zeta <- matrix(zeta, nrow = nrow(XB), ncol = length(thresholds), byrow = T)

# Make XB a matrix of N by T
  XB.t <- vector()
  for (tt in unique(dat$Time)) {
    XB.t  <- cbind(XB.t , XB[dat$Time == tt, , drop = F])
    
  }
  
 
# Simulate theta, a matrix of N by T 
  theta <- vector()
  for (ii in 1:nrow(XB.t)) {
    mu <- XB.t[ii,]
    theta.i <- MASS:::mvrnorm(n = 1, mu = mu, Sigma = cor.mat)
    theta <- rbind(theta, theta.i)
    
  }
  
  # Re-arrange theta in long format (length N*T)
  theta.long <- vector()
  for (time in 1:number.timepoints) {
    theta.long <- rbind(theta.long, theta[, time, drop = F])
    
  }
  

# Compare theta to zeta
 theta.long.cat <- matrix(theta.long, nrow = nrow(theta.long), ncol = length(thresholds), byrow = F)
 Y <- rowSums(theta.long.cat > zeta)
 dat$Y_comp <- Y

  # Note: re-arrange the dataset so the GEE can deal with it:
  dat <- dat[order(dat$USUBJID, dat$Time),] 
  
#-----------------------------------------------------------------------------------------  
  
  
  out <- list('dat' = dat, 
              'reg.formula' = reg.formula,
              'Beta' = Beta, 
              'thresholds' = thresholds,
              'cor.mat' = cor.mat 
  )
  
  return(out)
  
  
}## End "sim_data_ord.R" Code
