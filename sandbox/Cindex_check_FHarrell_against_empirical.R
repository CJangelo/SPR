
# 2.7.2022: Using F Harrell's example of a violated PO assumption to check t
# the empiircal C-index. It looks fine!
# So our empirical C-index is not tracking with the model estimates
# Unclear why, exactly
# Suspect it's just too strong of a violation of PO
#


#https://www.fharrell.com/post/wpo/



rm(list = ls())
gc()

w <- expand.grid(group=1:2, y=0:2)
n <- c(100, 110, 50, 10, 30, 60)
u <- w[rep(1:6, n),]
with(u, table(group, y))

library(rms)

or2 <- exp(coef(lrm(y == 2 ~ group, data=u))['group'])
or1 <- exp(coef(lrm(y >= 1 ~ group, data=u))['group'])
or12 <- exp(coef(lrm(y ~ group, data=u))['group'])
ors <- c(or2, or1, or12)
names(ors) <- c('y=2', 'y>=1', 'y')
ors

wilcox.test(y ~ group, u, correct=FALSE)


sumr1 <- with(u, sum(rank(y)[group == 1]))
sumr2 <- with(u, sum(rank(y)[group == 2]))
n1 <- sum(u$group == 1)
n2 <- sum(u$group == 2)
# wilcox.test uses sum of ranks in group 1
W <- sumr1 - n1 * (n1 + 1) / 2   # equals wilcox.test
W
W <- sumr2 - n2 * (n2 + 1) / 2
W / (n1 * n2)

g1 <- u[u$group == 1, 'y']
g2 <- u[u$group == 2, 'y']

g1.re <- sample(x = g1, size = 1e4, replace = T)
g2.re <- sample(x = g2, size = 1e4, replace = T)

C_emp <- ifelse(g2.re > g1.re, 1,
                   ifelse(g2.re < g1.re, 0, 0.5))
C_emp <- mean(C_emp, na.rm = T)

C_emp

library(ordinal)
mod.clm <- ordinal::clm(as.factor(y) ~ group, nAGQ = 5, data= u)
beta.hat <- coef(mod.clm)['group']
#emm.clm <- emmeans::emmeans(mod.clm, pairwise ~ group)
#emm.clm

  C_clm <- 1/(1 + exp(-1*beta.hat/1.52))

  C_clm
  C_emp
  W / (n1 * n2)
  OR <- exp(beta.hat)
  OR ^ 0.66/(1 + OR ^ 0.66)
