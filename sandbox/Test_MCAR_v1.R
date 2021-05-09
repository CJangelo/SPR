
# MCAR test


# Structure the Data:
score
dat1 <- dat[dat$Time == 'Time_3', score]
dat2 <- dat[dat$Time == 'Time_4', score]
dropped.out <- !is.na(dat1) & is.na(dat2) # these subjects dropped out at timepoint 4
table(dropped.out)

# T-test that compares the mean score at timepoint 3 of 
# Group 1: subjects who dropped out at timepoint 4
# Group 2: subjects who did NOT drop out at timepoint 4
# in other words, does score at timepoint 3 diff between these groups?
# t.test(x = dat1[dropped.out], y = dat1[!dropped.out])


# Does the score at the previous timepoint predict drop-out?
mod <- glm(dropped.out ~ dat1, family = 'binomial')
summary(mod)

# Does the score at the previous timepoint predict drop-out AFTER adjusting for biomarker?
biomarker <- dat$Bio[dat$Time == 'Time_4']
mod.mcar <- glm(dropped.out ~ dat1 + biomarker, family = 'binomial')
summary(mod.mcar)

