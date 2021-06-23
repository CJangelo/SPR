
# R code corresponding to Table 2 Fallacy paper
# Continuous Variables

# Westreich D, Greenland S. The table 2 fallacy: presenting and interpreting
# confounder and modifier coefficients. Am J Epidemiol.
# 2013 Feb 15;177(4):292-8. doi: 10.1093/aje/kws412.
# Epub 2013 Jan 30. PMID: 23371353; PMCID: PMC3626058.

# https://pubmed.ncbi.nlm.nih.gov/23371353/

#--------------
rm(list = ls())
gc()

set.seed(5232021)
n = 1e5
age <- rnorm(n = n, mean = 0, sd = 1) #standardize age so mean of 40y/o is zero
smoking <- 2*age + rnorm(n = n, mean = 0, sd = 1)
hiv <- 3*age + 2*smoking + rnorm(n = n, mean = 0, sd = 1)
stroke <- 2*age + 2*smoking + 3*hiv + rnorm(n = n, mean = 0, sd = 0.3)

# Model 1:
mod1 <- lm(stroke ~  age + smoking + hiv)
summary(mod1)
coef(mod1)
# g2g, Recovers generating paramers

# Include U in the generation of smoking, rest are same:
U <- rnorm(n = n, mean = 0, sd = 1)
age <- rnorm(n = n, mean = 0, sd = 1)
smoking <- 4*U + 2*age + rnorm(n = n, mean = 0, sd = 1)
hiv <- 3*age + 4*smoking + rnorm(n = n, mean = 0, sd = 1)
stroke <- 2*age + 2*smoking + 3*hiv + rnorm(n = n, mean = 0, sd = 1)

# Model 2:
mod.check <- lm(stroke ~ age + smoking + hiv)
coef(mod.check)
# Also g2g, note that model not in the manuscript! This is to check.

# Figure 2: U also contributes to stroke
stroke <- 2*U + 2*age + 2*smoking + 3*hiv + rnorm(n = n, mean = 0, sd = 1)
mod3 <- lm(stroke ~ age + smoking + hiv)
coef(mod3)
# hiv is correct, age and smoking are inaccurate
# Corresponds with the manuscript


# Include the variable U
mod4 <- lm(stroke ~ U + age + smoking + hiv)
coef(mod4)
# all estimates correct again - age and smoking, not just hiv

# heterogeneity
# Model 3:
stroke <-
  2*age + 2*smoking + 3*hiv + # main effects
  1*hiv*smoking + 1*hiv*age + 1*smoking*age + #interactions
  rnorm(n = n, mean = 0, sd = 1) # error term
mod5 <- lm(stroke ~ age + smoking + hiv)
coef(mod5)
# This estimate seems to be at odds with the following on page 5:
# "Thus, the total effect of HIV on the cohort would still be given by
# smoking-and-age standardization, and smoking-and-age–
# specific ffects of HIV on the log odds of stroke would still
# equal β1 + β4 × Smoking + β5×Age."

# try standardizing variables:
age <- as.vector(scale(age, center = T, scale = T))
smoking <- as.vector(scale(smoking, center = T, scale = T))
hiv <- as.vector(scale(hiv, center = T, scale = T))
#U <- as.vector(scale(U, center = T, scale = T))

stroke <-
  2*age + 2*smoking + 3*hiv +
  1*hiv*smoking + 1*hiv*age +
  1*smoking*age +
  rnorm(n = n, mean = 0, sd = 1)

#stroke <- as.vector(scale(stroke, center = T, scale = T))

mod6 <- lm(stroke ~ age + smoking + hiv)
coef(mod6)
# still not correct, what is going on here?
mod7 <- lm(stroke ~ age + smoking + hiv + hiv*smoking + hiv*age + smoking*age)
coef(mod7)

