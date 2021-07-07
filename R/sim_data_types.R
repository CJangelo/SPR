#' Simulate Conditional Longitudinal Data
#' This is a conditional model with a random intercept
#'
#' Select the data type you want to generate - this is best suited for
#' fitting data in `glmmTMB()` among other R packages that are based on
#' the conditional model with a random intercept
#'
#' @param data.type Specify the type of data you want to generate. Select
#' one from the following: 'Ordinal', 'Beta', 'Gaussian', 'NegBinom',
#' 'ZIP', 'ZINB', 'Multinom'. This must be a character vector of length 1.
#' @param reg.formula this is a regression formula to pass to the data generation; null is
#' formula(~ Group + Time + Time*Group),
#' @param Beta this is the Beta for the regression equation; numeric matrix of Beta values OR a scalar value != 0
#' that will be the final value of the interaction parameters, default is all(Beta == 0) for Type I error simulations
#' @param thresholds the ordinal generation uses a logistic approach,
#' and these are the tresholds
#' @param subject.var the subject-level variance, also known as the variance
#' associated with the random intercept
#' @param sigma2 the residual error in Gaussian data
#' @param shape2 parameter in Beta distribution
#' @param phi parameter in Beta distribution
#' @param theta size parameter in NB
#' @param zip zero inflation parameter
#' @param thresholds thresholds in ordinal data
#' @param k number of categories in ordinal/nominal
#' @return returns a dataframe containing the simulated data
#'
#' @export


sim_dat_types <- function(N = 1000,
                          data.type = NULL,
                          number.groups = 2,
                          number.timepoints = 4,
                          reg.formula = NULL,
                          Beta = NULL,
                          sigma2 = NULL,
                          shape2 = NULL,
                          phi = NULL,
                          theta = NULL,
                          zip = NULL,
                          thresholds = NULL,
                          k = NULL,
                          subject.var = 1
){

  # checks:
  if (is.null(data.type)) stop("Specify the data type you want and pass parameters")
  if (is.null(reg.formula)) { reg.formula <-  formula(~ Group + Time + Time*Group) }



#--------------------------------------

  dat <- data.frame(
    'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
    'id.geepack' = rep(1:N, length.out= N*number.timepoints),
    'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
    'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
    stringsAsFactors=F)


  # Biomarker:
  # if (Covariate == T) {
  #   dat$Covariate <- rnorm(N*number.timepoints, mean = 0, sd= 1) # No differences in Biomarker across Groups
  #   # Note: This is similar to having randomized Biomarker levels across arms in a RCT
  # }



  # Beta parameter default is zero - permits Type I error simulations
  # if (cond.mcar == F) {


  # ----------------------------------------------------------------------------------------
  # Design Matrix
  X <- model.matrix( reg.formula, data = dat)
  #-------------------------------------------------------------------------------------

  #---------------------------------------------------------------
  # Set some defaults for these parameter values if they're not passed:

  #
  if (data.type == 'Ordinal') {
    if (is.null(thresholds)) thresholds <- c(-2, 0, 2)
  }

  #
  if (data.type == 'Beta') {
    if (is.null(shape2)) shape2 <- 5
    if (is.null(phi)) phi <- 10
  }

  #
  if (data.type == 'Gaussian') {
    if (is.null(sigma2)) sigma2 <- 3
  }

  #
  if (data.type == 'NegBinom') {
    if (is.null(theta)) theta <- 1
  }

  if (data.type %in% c('ZIP', 'ZINB')) {
    if (is.null(theta)) theta <- 1
    if (is.null(zip)) zip <- 1
  }

  if (data.type == 'Multinom') {
    if (is.null(k)) k <- 5
  }

  #-------------
  # Beta
  if (is.null(Beta)) {
    if (data.type != 'Multinom') {
      bv <- seq(0, 1, length.out = ncol(X))
      Beta <- matrix(bv, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    }

    if (data.type == 'Multinom') {
      bv <- seq(0, 1, length.out = ncol(X))
      Beta <- matrix(bv, nrow = ncol(X), ncol = k-1,
                     byrow = F,
                     dimnames=list(colnames(X),
                                   paste0('param', 1:(k-1))))

    }
  }# end null Beta
  #---------------------------------------------------------------------

  # if (cond.mcar == T) {
  #
  #   # Generate Biomarker:
  #   dat$Covariate <- rnorm(N*number.timepoints, mean = 1*(dat$Group == 'Group_2'), sd= 1) # Biomarker differs across Groups
  #   # Design Matrix - Condition MCAR
  #   reg.formula <- formula(~ Group + Time + Covariate + Time*Group + Covariate*Time)
  #   X <- model.matrix( reg.formula, data = dat) # include Biomarker drop-out!
  #   # Conditional MCAR - biomarker affects Y
  #   beta.values <- seq(0.25, Beta, length.out = number.timepoints - 1)
  #   Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
  #   Beta[grepl('Covariate', rownames(Beta)) & grepl('Time', rownames(Beta)), ] <- -0.5*beta.values
  #   Beta[grepl('Group_2', rownames(Beta)) & grepl('Time', rownames(Beta)), ] <- beta.values
  #
  # }


  #----------------------------------------------
  # Continuous Data - Normal error distribution
  if (data.type == 'Gaussian') {

    XB <- X %*% Beta
    mat.XB <- XB # already a matrix??
    Zi <- rnorm(n = N, mean = 0, sd = sqrt(subject.var))
    mat.Zi <- matrix(Zi, nrow = nrow(mat.XB), ncol = ncol(mat.XB), byrow = F)

    eta <- mat.XB + mat.Zi
    eta <- as.vector(eta)
    error <- rnorm(n = N*number.timepoints, mean = 0, sd = sqrt(sigma2))
    dat$Y_gaussian <- eta + error

  }

  #-----------------------------------------------
  # Ordinal Data
  if (data.type == 'Ordinal') {
    # must pass thresholds!
    XB <- X %*% Beta
    mat.XB <- matrix(XB, nrow = nrow(XB), ncol = length(thresholds), byrow = F)
    mat.thr <- matrix(thresholds, nrow = nrow(XB), ncol = length(thresholds), byrow = T)
    Zi <- rnorm(n = N, mean = 0, sd = sqrt(subject.var))
    mat.Zi <- matrix(Zi, nrow = nrow(mat.XB), ncol = ncol(mat.XB), byrow = F)
    eta <- mat.thr - mat.XB + mat.Zi
    p <- exp(eta)/(1 + exp(eta))
    dat$Y_ord <- as.vector(apply(runif(n = N*number.timepoints) > p, 1, sum))

  }



  #-----------------------------------------------
  # Multinomial Data
  if (data.type == 'Multinom') {
    # must pass thresholds!
    XB <- X %*% Beta
    mat.XB <- XB
    Zi <- rnorm(n = N, mean = 0, sd = sqrt(subject.var))
    mat.Zi <- matrix(Zi, nrow = nrow(mat.XB), ncol = ncol(mat.XB), byrow = F)
    eta <- mat.XB + mat.Zi
    sum.expXB <- apply(exp(XB), 1, sum)
    p <- exp(XB)/(1 + sum.expXB)
    param0 <-  1 - rowSums(p)
    p <- cbind(param0, p)
    out <- vector()
    for(i in 1:nrow(p)){
      out <- c(out,
               sample(x = paste0('Cat_', 1:k), size = 1, prob = p[i, ])
      )
    }#end loop

    dat$Y_nom <- out

  }# end multinomial




  #----------------------------------------------------
  # Poisson and Negative Binomial:
  if (data.type %in% c('Poisson', 'NegBinom', 'ZINB', 'ZIP')) {

    XB <- X %*% Beta
    mat.XB <- XB # already a matrix??
    Zi <- rnorm(n = N, mean = 0, sd = sqrt(subject.var))
    mat.Zi <- matrix(Zi, nrow = nrow(mat.XB), ncol = ncol(mat.XB), byrow = F)

    eta <- mat.XB + mat.Zi
    mu <- exp(eta)
    mu <- as.vector(mu)
    dat$mu <- mu

    if (data.type == 'NegBinom') {
      dat$Y_nb <- rnbinom(n = length(dat$mu), size = theta, mu = mu)
    }

    if (data.type == 'Poisson') {
      dat$Y_pois <- rpois(n = length(dat$mu), lambda = mu)
    }

    if (data.type == 'ZIP') {
      # Zero-inflation ADDED:
      P_zi <- 1/(1 + exp(-zip)) #scalar
      Y <- 0:100
      Y_zip <- rep(NA, length(mu))
      for (i in 1:length(mu)){
        mu.i <- mu[i]
        prob0 <- log(P_zi + (1 - P_zi)*exp(-mu.i))
        # Poisson, Not Negative Binomial
        prob1 <- log(1 - P_zi) +  Y * log(mu.i) - mu.i - lgamma(Y + 1)
        logP <- (Y == 0) * prob0  +   (Y != 0) *  prob1
        P <- exp(logP)
        P <- P/sum(P)
        Y_zip[i] <- sample(x = Y, size = 1, prob = P)
      }
      dat$Y_zip <- as.vector(Y_zip)

    } #end if ZIP


    if (data.type == 'ZINB') {

      size <- theta  # dispersion parameter
      prob <- theta/(theta + mu)

      # Zero-inflation ADDED:
      P_zi <- 1/(1 + exp(-zip)) #scalar
      Y <- 0:100 # should be 0 to infinity
      #--------------------------
      Y_zinb <- rep(NA, length(prob))
      for (i in 1:length(prob)){
        prob.i <- prob[i]
        # If response = 0
        # Can be 0 two different ways, either inflation OR from the count process
        # 'OR' means you add the probabilities (then take the log)
        prob0 <- log(P_zi + (1-P_zi)*prob.i^size)
        # If response > 0
        # multiply probability that it's NOT inflated zero times count process
        # This is an "AND" statement, requires multiplying probabilities
        prob1 <- log(1 - P_zi) + lgamma(Y + size) - lgamma(size) -
          lgamma(Y + 1) + size*log(prob.i) + (Y)*log(1-prob.i)

        logP <- (Y == 0) * prob0  +   (Y != 0) *  prob1
        P <- exp(logP)
        P <- P/sum(P) # probabilities sum to 1
        Y_zinb[i] <- sample(x = Y, size = 1, prob = P)

      }# end loop over subjects by timepoints

      dat$Y_zinb <- as.vector(Y_zinb)
    } # end if ZINB
  }# end  ('Poisson', 'NegBinom', 'ZINB', 'ZIP')



  #-----------------------------------------------
  # Beta Distribution
  if (data.type == 'Beta') {

    XB <- X %*% Beta
    mat.XB <- XB # already a matrix??
    Zi <- rnorm(n = N, mean = 0, sd = sqrt(subject.var))
    mat.Zi <- matrix(Zi, nrow = nrow(mat.XB), ncol = ncol(mat.XB), byrow = F)

    mu.i <- as.vector(mat.XB + mat.Zi)
    mu.i <- exp(mu.i)/(1 + exp(mu.i))

    # Re-parameterize, see pg 3 of PDF
    shape1.i <- mu.i*phi
    shape2.i <- phi - mu.i*phi

    # Generate data:
    dat$Y_beta <- stats::rbeta(n = length(shape1.i),
                               shape1 = shape1.i,
                               shape2 = shape2.i)

  }



  #------------------------------------------------------------------------


  out <- list('dat' = dat,
              'N' = N,
              'data.type' = data.type,
              'Beta' = Beta,
              'subject.var' = subject.var,
              'sigma2' = sigma2,
              'phi' = phi,
              'shape2' = shape2,
              'theta' = theta,
              'thresholds' = thresholds,
              'k' = k
  )

  return(out)


}##

