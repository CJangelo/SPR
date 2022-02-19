

sim_PRO_data <- function(){


  # simulate the data
sim.out <- sim_surv_data()

dat.surv <- sim.out$dat.surv
dat.PRO <- sim.out$dat.PRO

N = length(unique(dat.PRO$USUBJID))
#number.groups = 2
#number.timepoints = length(unique(dat.PRO$Time_factor))
number.timepoints = length(unique(dat.PRO$Visit))
#reg.formula =  formula(~ Group + Time_factor + Time_factor*Group)
reg.formula =  formula(~ Group + Visit + Visit*Group)
#Beta = 0
thresholds = c(-2,0, 2)
#subject.var = 1

 # ----------------------------------------------------------------------------------------
  # Design Matrix
  X <- model.matrix( reg.formula, data = dat.PRO)
  #-------------------------------------------------------------------------------------
beta.values <- seq(-0.2, -0.7, length.out = number.timepoints - 1)
Beta <- matrix(0,nrow = ncol(X),dimnames = list(colnames(X), 'param'))
#Beta[grepl('Group_2', rownames(Beta)) & grepl('Time', rownames(Beta)),] <- beta.values
Beta[grepl('Group_2', rownames(Beta)) & grepl('Visit', rownames(Beta)),] <- beta.values


  # --------------------------------------------------------------------------------------------------
  # Matrix multiply:
  rownames(Beta) <- colnames(X)
  XB <- X %*% Beta

  mat.XB <- matrix(XB, nrow = nrow(XB), ncol = length(thresholds), byrow = F)
  mat.thr <- matrix(thresholds, nrow = nrow(XB), ncol = length(thresholds), byrow = T)

# - Conditional Model:
#Zi <- rnorm(n = N, mean = 0, sd = sqrt(subject.var))
# Use the random effect shared with the survival model instead here:
mat.Zi <- matrix(dat.PRO$Zi, nrow = nrow(XB), ncol = length(thresholds), byrow = F)
eta <- mat.thr - mat.XB - mat.Zi
#eta <- mat.thr - mat.XB - mat.Zi
p <- exp(eta)/(1 + exp(eta))
dat.PRO$Y_comp <- as.vector(apply(runif(n = N*number.timepoints) > p, 1, sum))

out <- list('dat.PRO' = dat.PRO, 'dat.surv' = dat.surv,
            'Beta_PRO' = Beta, 'thresholds' = thresholds )

return(out)

}
