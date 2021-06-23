#' Simulate Conditional Longitudinal Ordinal Data with Logistic Specification
#' This is a conditional model with a random intercept
#'
#' TODO: Build in check if the correlations are too high! Where is that code?
#' Note that the variances aren't used if passed,
#' Requires MASS R package
#' @param reg.formula this is a regression formula to pass to the data generation; null is
#' formula(~ Group + Time + Time*Group),
#' @param Beta this is the Beta for the regression equation; numeric matrix of Beta values OR a scalar value != 0
#' that will be the final value of the interaction parameters, default is all(Beta == 0) for Type I error simulations
#' @param thresholds the ordinal generation uses a logistic approach,
#' and these are the tresholds
#' @param subject.var the subject-level variance, also known as the variance
#' associated with the random intercept
#' @param cond.mcar logical; do you want a conditional MCAR data generation
#' @param Covariate logical; do you want a random simulated covariate?
#'
#' @return returns a dataframe containing the simulated data
#' @export



sim_dat_ord_logistic_conditional <- function(N = 500 ,
                                             number.groups = 2,
                                             number.timepoints = 4,
                                             reg.formula = NULL,
                                             Beta = 0,
                                             thresholds = c(-2, -1, 0, 1),
                                             subject.var = 1,
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



# ---------------------------------------------------------------------
# GENERATE ORDINAL DATA - LOGISTIC
#  This chunk here is the only part that is distinct from "sim_dat.R"

mat.XB <- matrix(XB, nrow = nrow(XB), ncol = length(thresholds), byrow = F)
mat.thr <- matrix(thresholds, nrow = nrow(XB), ncol = length(thresholds), byrow = T)
#eta <- mat.thr - mat.XB

# - Conditional Model:
# How many of each? Is the subject level variance only 1 per subject?
# And the residual is 1 per subject per timepoint?
  Zi <- rnorm(n = N, mean = 0, sd = sqrt(subject.var))
    mat.Zi <- matrix(Zi, nrow = nrow(XB), ncol = length(thresholds), byrow = F)
 # ei <- rnorm(n = N*number.timepoints, mean = 0, sd = sqrt(residual.var))
  #   mat.ei <- matrix(ei, nrow = nrow(XB), ncol = length(thresholds), byrow = F)

eta <- mat.thr - mat.XB + mat.Zi #+ mat.ei
p <- exp(eta)/(1 + exp(eta))
dat$Y_comp <- as.vector(apply(runif(n = N*number.timepoints) > p, 1, sum))

#------------------------------------------------------------------------


  out <- list('dat' = dat,
              'reg.formula' = reg.formula,
              'Beta' = Beta,
              'thresholds' = thresholds
  )

  return(out)


}##

