

# TODO: figure out how Weibull is parameterized
# I'm using this as:
# https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Cox-Regression.pdf
# Set rho = 0 to get the Cox Prop Haz regression mode


#---
rm(list = ls())
gc()
N = 1e4 #this should be divisible by however many groups you use!
number.groups <- 2
number.timepoints <- 1
set.seed(10182021)

dat <- data.frame('USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Group' = rep(paste0('Group_', 1:number.groups), length.out = N*number.timepoints),
                  stringsAsFactors=F)


# Create Beta parameters for these design matrix:
X <- model.matrix( ~ Group  , data = dat)

# Create Beta
Beta <- matrix(0, nrow = ncol(X), dimnames=list(colnames(X), 'param'))
  Beta[] <- c(-5, -1)

# Matrix multiply:
XB <- X %*% Beta

# Hazard Ratio
#ht <- exp(XB)
#range(ht)
dat$XB <- as.vector(XB)

# Data mgmt
#dat$ht <- as.vector(ht)
dat$time <- NA
dat$status <- NA


# Hazard Function, with density function as follows:
pdf_surv <- function(tt, vv, rho) {

  ht <- exp(vv + rho*log(tt))
  pdf <- ht*exp(-1*ht*tt)
  return(pdf)

}

# tt is the time
# ht is the subject-specific hazard ratio
# https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Cox-Regression.pdf

# Probability of event over 6 months, 180 days:
tt <- 1:180
i <- 1

# Loop over subjects
for (i in 1:nrow(dat)) {

  pp <- pdf_surv(tt = tt, rho = 0, vv = dat$XB[i]) # prob for subject i
  pp <- cumsum(pp) # it's a cumulative process
  out <- runif(n = 1) < pp # random draw

  # Censored:
  if (sum(out) == 0) {
    dat$status[i] <- 1 # censored
    dat$time[i] <- tail(tt, 1) # time of censoring
  }

  # Dead:
  if (sum(out) > 0) {
    dat$status[i] <- 2 # dead
    dat$time[i] <- tt[which(cumsum(out) == 1)] #survival time in months
  }


}# end loop over i subjects


#-------------------
# Analyze - confirm can recover gen param
library(survival)
fit <- survfit(Surv(time, status) ~ Group, data = dat)
print(fit)
library(ggplot2)
library(survminer)
ggsurvplot(fit,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw()) # Change ggplot2 theme


mod1 <- coxph(Surv(time, status) ~ Group, data = dat)
summary(mod1)

test.ph <- cox.zph(mod1)
test.ph

mod2 <- survreg(Surv(time, status) ~ Group,
                dist= 'weibull',
                data = dat)
summary(mod2)
