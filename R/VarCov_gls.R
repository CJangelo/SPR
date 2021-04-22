#' Extract Variance-Covariance Matrix from nlme object
#'
#' Requires nlme object
#' @param fit nlme model, specifically a MMRM model estimated using `gls()` function
#' @return returns a matrix of the estimated variance covariance values
#' @export



VarCov_gls <- function(fit){
  
  vars <- coef(fit$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE)^2 * fit$sigma^2
  r = coef(fit$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)
  cors = matrix(NA, ncol = length(vars), nrow = length(vars))
	cors[lower.tri(cors)] = r
	cors[upper.tri(cors)] = t(cors)[upper.tri(t(cors))]
	diag(cors) = rep(1, length(vars))
  
  covs <- diag(sqrt(vars) ) %*% cors %*% diag(sqrt(vars))
  rownames(covs) <- names(vars)
  colnames(covs) <- names(vars)
  return(covs)
  
}