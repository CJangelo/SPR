

# 11.30.21
# Simulate PRO data
# random slopes and random intercepts
# N = 1e3; number.groups = 2; number.timepoints = 4;
# reg.formula = formula( ~ Group*Time); sigma2 = 1

sim_dat_PRO_random_slopes <- function(N = NULL,
                                      number.groups = NULL,
                                      number.timepoints = NULL,
                                      reg.formula = NULL, #formula(~ Group + Time + Time*Group),
                                      Beta = NULL,
                                      sigma2 = NULL,
                                      random.sigma = NULL
){


  #--------------------------------------
  dat <- data.frame(
    'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
    'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
    'Time' = rep(seq(1, 4, length.out = number.timepoints), each = N),
    'Time_factor' = rep(paste0('Time_', 1:number.timepoints), each = N),
    stringsAsFactors=F)


  # ----------------------------------------------------------------------------------------
  # Design Matrix
  X <- model.matrix( reg.formula, data = dat)

  # Create Beta
  Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    Beta[] <- c(0, 0, -0.5, 0.5)
  #-------------------------------------------------------------------------------------

  # Fixed Effects
  XB <- X %*% Beta
  mat.XB <- XB

  # Random Effects - Time within Subject, Random slopes & random intercepts
  Z <- model.matrix( ~ Time, data = dat)
  random.sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = T)
  bb <- MASS::mvrnorm(n = N, mu = c(0, 0), Sigma = random.sigma)
  df.bb <- data.frame('USUBJID' = unique(dat$USUBJID), 'bb' = bb)
  dat$rn <- 1:nrow(dat)
  df.bb <- merge(x = dat, y = df.bb, by = 'USUBJID', all.x = T)
  df.bb <- df.bb[order(df.bb$rn), ]
  mat.bb <- as.matrix(df.bb[ , grep('bb', colnames(df.bb), value = T)])
  mat.Zi <- Z * mat.bb

  # Fixed and Random
  mt <- cbind(mat.XB, mat.Zi)
  mt <- apply(mt, 1, sum)
  mt <- as.vector(mt)

  # Error term
  error <- rnorm(n = N*number.timepoints, mean = 0, sd = sqrt(sigma2))

  # output the mt for the survival model:
  dat$mt <- mt
  dat$error <- error
  dat$Y <- mt + error


  out <- list('dat' = dat,
              'Beta' = Beta,
              'random.sigma' = random.sigma,
              'reg.formula' = reg.formula)

  return(out)


}## # end function

