

sim_dat_surv <- function(dat.PRO = NULL, reg.formula = NULL){

  st <- Sys.time()
  # ----------------------------------------------------------------------------------------
  # Design Matrix
  X <- model.matrix( reg.formula, data = dat.PRO)

  # Create Beta
  Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    Beta[] <- c(-5, 1)

  # Matrix multiply:
  XB <- X %*% Beta
  dat.PRO$wi <- as.vector(XB)
  # use wi and mt in the survival process
  # will need to incorporate the PRO dataset, so put it there

  # Hazard Ratio
  h0 <- 1  # Baseline HR: I guess just set this value? TODO
  alpha <- 0 # Associated between PRO process and survival process

  # Create survival time dataset
  dat.surv <- dat.PRO[!duplicated(dat.PRO$USUBJID), ]
  dat.surv$time.surv <- NA
  dat.surv$status.surv <- NA


  # Hazard Function, with density function as follows:
  compute_ht <- function(dat.PRO, i, tt, h0, alpha){
    dat.it <- dat.PRO[dat.PRO$USUBJID == i & dat.PRO$Time <= tt, ]
    mt <- tail(dat.it$mt, 1)
    wi <- tail(dat.it$wi, 1)
    ht <- h0*exp(wi + alpha*mt)
    return(ht)
  }


# Probability of event over given time span:
Time <- seq(min(dat.PRO$Time), max(dat.PRO$Time), length.out = 100)
id <- unique(dat.PRO$USUBJID)


# Loop over subjects
for (i in id ) {

    pp <- vector() # prob for subject i
  for (tt in Time) {
    pp <- c(pp, compute_ht(dat.PRO = dat.PRO, i = i, tt = tt, h0 = h0, alpha = alpha) )
  }

  St <- exp(-1*cumsum(pp)) # pg 38 of JM textbook
  surv <- runif(n = 1) < (1 - St)

  # Censored:
  if (sum(surv) == 0) {
    dat.surv[dat.surv$USUBJID == i, 'status.surv'] <- 1 # censored
    dat.surv[dat.surv$USUBJID == i, 'time.surv'] <- tail(Time, 1) # time of censoring
  }

  # Dead:
  if (sum(surv) > 0) {
    dat.surv[dat.surv$USUBJID == i, 'status.surv'] <- 2 # dead
    dat.surv[dat.surv$USUBJID == i, 'time.surv'] <- Time[which(cumsum(surv) == 1)] #survival time
  }


}# end loop over i subjects


  out <- list('dat.PRO' = dat.PRO, 'dat.surv' = dat.surv)
  et <- Sys.time()
  cat(paste0('Run time: ', round(et - st), '\n'))


  return(out)

}


