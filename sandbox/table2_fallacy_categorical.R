
# R code corresponding to Table 2 Fallacy paper
# categorical data

# Notes on generation of binary data
# 1. With large coefficients (2, 3, 4) all estimates are bad. Just doesn't work.
# 2. See note below about how age and U have to have a strong influence on
# smoking in order to influence the estimates in model 3
# Westreich D, Greenland S. The table 2 fallacy: presenting and interpreting
# confounder and modifier coefficients. Am J Epidemiol.
# 2013 Feb 15;177(4):292-8. doi: 10.1093/aje/kws412.
# Epub 2013 Jan 30. PMID: 23371353; PMCID: PMC3626058.

# https://pubmed.ncbi.nlm.nih.gov/23371353/

#--------------
rm(list = ls())
gc()

set.seed(5232021)
n = 1e4
age <- rnorm(n = n, mean = 0, sd = 1)
  slope <- .2*age #+ rnorm(n = n, mean = 0, sd = .1)
  p <- 1/(1+exp((slope - mean(slope))))
  smoking <- 1*(runif(n = n) > p)
    slope <- .3*age + .2*smoking# + rnorm(n = n, mean = 0, sd = .1)
    p <- 1/(1+exp((slope - mean(slope))))
    hiv <- 1*(runif(n = n) > p)
      slope <- 1*age + 1*smoking + 1*hiv #+ rnorm(n = n, mean = 0, sd = .1)
      p <- 1/(1+exp((slope - mean(slope))))
      stroke <- 1*runif(n = n) > p

# Model 1:
mod1 <- glm(stroke ~  age + smoking + hiv, family = 'binomial')
coef(mod1)
# g2g, recovers gen param

# Include U in the generation of smoking, rest are same:
U <- rnorm(n = n, mean = 0, sd = 1)
age <- rnorm(n = n, mean = 0, sd = 1) #standardize age so mean of 40y/o is zero
  #slope <- .2*age + .2*U #+ rnorm(n = n, mean = 0, sd = .1)
# This line is the other interesting point - if the coefficients here
# are small, then the U and age don't influence stroke enough to matter!
  slope <- 4*age + 4*U #+ rnorm(n = n, mean = 0, sd = .1)
  p <- 1/(1+exp((slope - mean(slope))))
  smoking <- 1*(runif(n = n) > p)
    slope <- 1*age + 1*smoking #+ rnorm(n = n, mean = 0, sd = .1)
    p <- 1/(1+exp((slope - mean(slope))))
    hiv <- 1*(runif(n = n) > p)
      slope <- 1*age + 1*smoking + 1*hiv #+ rnorm(n = n, mean = 0, sd = 0.1)
      p <- 1/(1+exp((slope - mean(slope))))
      stroke <- 1*(runif(n = n) > p)

# Model 2:
mod.check <- glm(stroke ~ age + smoking + hiv, family = 'binomial')
coef(mod.check)
# Also g2g, note that model not in the manuscript! This is to check.

# save these variables:
df1 <- data.frame(stroke, age, smoking, hiv, stringsAsFactors = F)
colnames(df1) <- paste0(colnames(df1), '_Fig1')

# Figure 2: U also contributes to stroke
slope <- 1*U + 1*age + 1*smoking + 1*hiv #+ rnorm(n = n, mean = 0, sd = 0.1)
p <- 1/(1+exp((slope - mean(slope))))
stroke2 <- 1*runif(n = n) > p

# Model 3
mod3 <- glm(stroke2 ~ age + smoking + hiv, family = 'binomial')
coef(mod3)
# smoking p = 1/(1 + exp(3*age + 3*U)); these have to be 3, can't be like 0.2
# if they're small, then there's only some influence on smoking, and then
# there's less influence on stroke via smoking

# Include the variable U
mod4 <- glm(stroke2 ~ U + age + smoking + hiv, family = 'binomial')
coef(mod4)
# estimates correct

#---------------------------------
# write it out:

df2 <- data.frame(stroke, U, age, smoking, hiv, stringsAsFactors = F)
colnames(df2) <- paste0(colnames(df2), '_Fig2')
df <- cbind.data.frame(df1, df2)

write.table(df, na = '.', quote = F,
            sep = ', ', col.names = T,
            row.names = F, file = 'data_table2_fallacy.txt')
# covers Figures 1 and 2, up through Page 3 to the beginning
# of the section "INTERPRETATION GETS HARDER WITH HETEROGENEITY"


#----------------------------------------
# heterogeneity
# Model 3:
slope <-
  .2*age + .2*smoking + .3*hiv + # main effects
  .1*hiv*smoking + .1*hiv*age + .1*smoking*age  #interactions
p <- 1/(1+exp((slope - mean(slope))))
stroke3 <- 1*runif(n = n) > p
mod5 <- glm(stroke3 ~ age + smoking + hiv, family = 'binomial')
coef(mod5)


mod6 <- glm(stroke3 ~ age + smoking + hiv +
              hiv*smoking + hiv*age + smoking*age, family = 'binomial')
coef(mod6)
