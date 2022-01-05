

# 11.24.2021
# Simulate data for the Joint Model
# This is the PRO data part
# output the mt


sim_dat_PRO <- function(N = NULL,
                          number.groups = NULL,
                          number.timepoints = NULL,
                          reg.formula = NULL, #formula(~ Group + Time + Time*Group),
                          Beta = NULL,
                          sigma2 = NULL,
                          subject.var = NULL
){


#--------------------------------------
  dat <- data.frame(
    'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
   # 'id.geepack' = rep(1:N, length.out= N*number.timepoints),
    'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
    # 'Time' = rep(seq(0, 180, by = 60), each = N),
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

# Generate data:
    XB <- X %*% Beta
    mat.XB <- XB # already a matrix??
    Zi <- rnorm(n = N, mean = 0, sd = sqrt(subject.var))
    mat.Zi <- matrix(Zi, nrow = nrow(mat.XB), ncol = ncol(mat.XB), byrow = F)

    mt <- mat.XB + mat.Zi
    mt <- as.vector(mt)
    error <- rnorm(n = N*number.timepoints, mean = 0, sd = sqrt(sigma2))

   # output the mt for the survival model:
    dat$mt <- mt
    dat$error <- error
    dat$Y <- mt + error


  out <- list('dat' = dat,
              'Beta' = Beta,
              'reg.formula' = reg.formula)

  return(out)


}## # end function

