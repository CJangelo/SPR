
#---
#rm(list = ls())
#gc()

#library(SPR)

sim_surv_data <- function(){

    st <- Sys.time()

  # Update this later, also need to pass Beta
  N = 260
  number.groups = 2
  number.timepoints = 4
  reg.formula =  formula(~ Group)
  subject.var = 1

#------------------------------------------------------------------------------------------------------------------

  dat <- data.frame(
    'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
    'id.geepack' = rep(1:N, length.out= N*number.timepoints),
    'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
    # 'Time' = rep(seq(1, 4, length.out = number.timepoints), each = N),
    'Time' = rep(seq(0, 3, length.out = number.timepoints), each = N),
    # 'Time_factor' = rep(paste0('Time_', 1:number.timepoints), each = N),
    'Visit' = rep(paste0('Visit_', 1:number.timepoints), each = N),
   # 'Y_comp' = rep(NA, N*number.timepoints),
    stringsAsFactors=F)


  # ----------------------------------------------------------------------------------------
  # Design Matrix
  X <- model.matrix( reg.formula, data = dat)

  # Create Beta
  Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
    Beta[] <- c(-5, -0.7)

  # Matrix multiply:
  XB <- X %*% Beta
  dat$XB_surv <- as.vector(XB)
  # - Conditional Model:
  dat.Zi <- data.frame('USUBJID' = unique(dat$USUBJID),
                       'Zi' = rnorm(n = N, mean = 0, sd = sqrt(subject.var)))
  dat <- merge(x = dat, y = dat.Zi, by = 'USUBJID')
  #dat$Zi <- as.vector(Zi)
  dat$mt <- dat$XB_surv + dat$Zi
  # Eh this might be a bit much, just re-use the random effect:
  #dat$mt <- dat$Zi
  #mat.Zi <- matrix(Zi, nrow = nrow(XB), ncol = length(thresholds), byrow = F)

  # Hazard Ratio
  h0 <- 1  # Baseline HR: I guess just set this value? TODO
  #alpha <- 1 # Associated between PRO process and survival process

  # Create survival time dataset
  dat.surv <- dat[!duplicated(dat$USUBJID), ]
  dat.surv$time.surv <- NA
  dat.surv$status.surv <- NA


  # Hazard Function, with density function as follows:
  compute_ht <- function(dat, i, tt, h0){
    dat.it <- dat[dat$USUBJID == i & dat$Time <= tt, ]
    mt <- tail(dat.it$mt, 1)
    ht <- h0*exp(mt) # Now it's reversed, you have to pass mt to the PRO generation
    #wi <- tail(dat.it$wi, 1)
    #ht <- h0*exp(wi + alpha*mt)
    return(ht)
  }


# Probability of event over given time span:
Time <- seq(min(dat$Time), max(dat$Time), length.out = 100)
id <- unique(dat$USUBJID)


# Loop over subjects
for (i in id ) {

    pp <- vector() # prob for subject i
  for (tt in Time) {
    pp <- c(pp, compute_ht(dat = dat, i = i, tt = tt, h0 = h0) )
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


  #View(dat[order(dat$USUBJID), ])
  #View(dat.surv[order(dat.surv$USUBJID), ])
  out <- list('dat.PRO' = dat, 'dat.surv' = dat.surv)
  et <- Sys.time()
  tmp <- et - st
  print(round(et - st))
  #cat(paste0('Run time: ', print(round(et - st)), '\n'))


  return(out)

}


