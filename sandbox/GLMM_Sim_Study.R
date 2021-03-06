

# 5.28.21
# First step is to recover estimates with full data
# use GEE
# use GLMM and then convert back to population average estimate
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650580/#R11
# https://projecteuclid.org/download/pdf_1/euclid.ss/1009212671

# Then figure out a way to compute the standard errors of the population average
# unclear how to do that - that's only thing preventing you from abandoning
# the GEE!

rm(list = ls())
gc()

library(MASS)
library(SPR)
library(geepack)
library(lme4)

# Initialize outputs: out1 is for the means/medians
est.names <- c('GLM', 'GEE', 'GLMM', 'GLMM_pop_avg_logit','GLMM_pop_avg_probit')
number.param <- 8 # see Beta
out <- replicate(n = length(est.names),
          expr = as.data.frame(matrix(ncol=number.param,nrow=0)),
          simplify = F)
  names(out) <- est.names

  number.repl <- 100

for (repl in 1:number.repl) {


  set.seed(5282021 + repl)

 sim.out <- SPR::sim_dat_binom(N = 1000,
                               number.groups = 2 ,
                               number.timepoints = 4,
                               reg.formula = formula( ~ Time + Group + Time*Group),
                               Beta = 1,
                               corr = 'ar1',
                               cor.value = 0.4)


  dat <- sim.out$dat

# GLM
  mod.glm1 <- stats::glm(Y_comp ~ Time + Group + Time*Group,
                         family = 'binomial',
                         data = dat)

  # Generalized Estimating Equations:
  mod.gee1 <- geepack::geeglm(Y_comp ~ Time + Group + Time*Group,
                              id = id.geepack,
                              data = dat,
                              family = binomial, corstr = "ar1")


  # GLMM
  mod.glmm <- lme4::glmer(Y_comp ~ Time + Group + Time*Group + (1|USUBJID),
                          family = 'binomial',
                          data = dat)

  sigma2 <- lme4::VarCorr(mod.glmm)
  sigma2 <- as.data.frame(sigma2)$vcov
  denom.logit <- sqrt( (sigma2 + (pi^2/3))/(pi^2/3) )
  denom.probit <- sqrt( 1 + sigma2)

      #rbind all this shit
      out[['GLM']][repl, ] <- coef(mod.glm1)
              #rbind.data.frame(out[['GLM']], coef(mod.glm1))
      out[['GEE']][repl, ] <- coef(mod.gee1)
        #rbind(out[['GEE']], coef(mod.gee1))
      out[['GLMM']][repl, ] <- fixef(mod.glmm)
        #rbind(out[['GLMM']], fixef(mod.glmm))
      out[['GLMM_pop_avg_logit']][repl, ] <- fixef(mod.glmm)/denom.logit
        #rbind(out[['GLMM_pop_avg_logit']], as.vector(fixef(mod.glmm)/denom.logit))
      out[['GLMM_pop_avg_probit']][repl, ] <- fixef(mod.glmm)/denom.probit
        #rbind(out[['GLMM_pop_avg_probit']], fixef(mod.glmm)/denom.probit)

  cat(paste0('Replication: ', repl, '\n'))

}# end replications



  out2 <- lapply(out, function(x) setNames(object = x, nm = rownames(sim.out$Beta)))
  tmp <- lapply(out2, colMeans)
  out3 <- cbind(sim.out$Beta, rbind.data.frame(tmp))
