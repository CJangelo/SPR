
# TODO:
# 1. sort out Weibull parameterization in other code file
# 2. consider adding simulation of discrete time survival analysis
# 3. Missing data!

# show that the beta values make it harder to recover
# Try the following:
# Beta[] <- c(-3, -1), c(-4, -1), and c(-5, -1)
# Look at the N and the number of events in summary(mod1)
# Estimation is wrong when the number of events is almost equal to N

# Resources:
# I think John Fox is great: https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Cox-Regression.pdf
# Very understandable, ties it into the R output

# Full textbook: http://sistemas.fciencias.unam.mx/~ediaz/Cursos/Estadistica3/Libros/0a9X.pdf
# This is a decent reference

# As far as R resources, I think this is a great series:
# 1. http://www.sthda.com/english/wiki/survival-analysis-basics
# 2. http://www.sthda.com/english/wiki/cox-proportional-hazards-model
# 3. http://www.sthda.com/english/wiki/cox-model-assumptions
# This makes it super easy to just mash a few buttons and get output



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
  Beta[] <- c(-4, -1)

# Matrix multiply:
XB <- X %*% Beta

# Hazard Ratio
ht <- exp(XB)
range(ht)

# Data mgmt
dat$ht <- as.vector(ht)
dat$time <- NA
dat$status <- NA


# Hazard Function, with density function as follows:
pdf_surv <- function(tt, ht) ht*exp(-1*ht*tt)
# tt is the time
# ht is the subject-specific hazard ratio
# https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Cox-Regression.pdf

# Probability of event over 6 months, 180 days:
tt <- 1:180
i <- 1

# Loop over subjects
for (i in 1:nrow(dat)) {

  pp <- pdf_surv(tt = tt, ht = dat$ht[i]) # prob for subject i
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


#dat1 <- dat[dat$Group == 'Group_1', ]
#xtabs(~ status + Group + time, data = dat)

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
          #palette = c("#E7B800", "#2E9FDF"))

#mod0 <- coxph(Surv(time, status) ~ 1, data = dat)
mod1 <- coxph(Surv(time, status) ~ Group, data = dat)
#anova(mod0, mod1, test = 'LRT')
summary(mod1)

tmp <- data.frame(fit$time, fit$n.risk, fit$n.event, fit$n.censor, fit$surv)
tmp$p <- tmp$fit.n.event/tmp$fit.n.risk
tmp1 <- tmp[tt,]
tmp2 <- tmp[(tt+tail(tt,1)),]
mean(tmp2$p/tmp1$p)
median(tmp2$p/tmp1$p)
exp(-1)


#---- Proportional hazards ratio holds? Test it
test.ph <- cox.zph(mod1)
test.ph
ggcoxzph(test.ph)
