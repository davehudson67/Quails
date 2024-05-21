library(tidyverse)
library(survminer)
library(survival)

rm(list=ls())

## load data
load("quails.RData")
quails <- read_csv("Survival.FINAL.csv")
quails$cens <- if_else(censored == 2, 0, 1)

## create survival object
fit <- survfit(Surv(quails$time, quails$cens) ~ quails$sex)

## survival curves
ggsurvplot(fit, data = quails, risk.table = FALSE,
           conf.int = TRUE, pval = TRUE, title = "Kaplan-Meier plot of survival probability",
           legend = "right", legend.title = "Sex", legend.labs = c("Female","Male"),
           pval.coord = c(200, 0.25), surv.median.line = "v",
           censor.shape=124, ggtheme = theme_bw(), font.main = c(16, "darkblue"), 
           subtitle = "Censored individuals shown by vertical lines, dashed line indicates median age.
           (p-value corresponds to log-rank comparison of survival curves)")


#Comparison of the 2 survival curves
fit2 <- Surv(mong$dur, mong$dead)
survdiff(fit2 ~ mong$sex)
survdiff(fit2 ~ mong$sex, rho = 1)

#Check smoothed Hazard rate curve
library(muhaz)
plot(muhaz(mong$dur, mong$dead))

#Plot life table hazard
library(KMsurv)
library(biostat3)

lt <- lifetab2(fit2 ~ 1, mong)
plot (lt$hazard, type = 'b')

#Cox ph
mong.fit <- coxph(fit2 ~ mong$sex)
summary(mong.fit)






