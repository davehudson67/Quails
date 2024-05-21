library(tidyverse)
library(survminer)
library(survival)
rm(list=ls())

## load data
quails <- read_csv("Survival.FINAL.csv")
sex <- as.numeric(as.factor(quails$sex))

## set up data
tU <- quails$time
tB <- rep(0, nrow(quails))
tL <- if_else(quails$time == 0, quails$time, quails$time - 1)

## check last seen against death
table(tL >= tU)
#summary(quails$birth_day[tL >= tU])
#summary(mong$death_day[tL < tU])

## set up censoring vector (1 = interval, 2 = right)
quails$censored <- if_else(quails$status == 1, 1, 2)
#quails$censored <- if_else(quails$status == 1, 1, 0)
censored <- quails$censored

## set death time for censored individuals
tU[quails$censored == 2] <- NA

## summaries
stopifnot(all(tL[!is.na(tU)] <= tU[!is.na(tU)]))
stopifnot(all(tB <= tL)) #some individuals died on the same day they were born

## define censoring matrices
cint <- cbind(tL, tU)
colnames(cint) <- NULL
cint[censored == 2, 2] <- cint[censored == 2, 1] 
cint[censored == 2, 1] <- 0

## check censoring times
summary(apply(cbind(tL, tU), 1, diff)[censored == 1])
summary(apply(cbind(tL, tU), 1, diff)[censored == 2])

## set up latent death times
tD <- rep(NA, length(tU))

## set up nind
nind <- nrow(cint)

## set up zL/zU input matrix
zL <- tL
zU <- tU

## produce K-M plot (assuming that interval-censored individuals
## die at the end time of the interval)
quails_km_dat <- select(quails, tD = time, tL = time, cens = censored) %>%
  mutate(tB = 0) %>%
  mutate(tD = ifelse(cens == 2, tL, tD)) %>%
  mutate(tD = tD - tB) %>%
  select(-tL, -tB) %>%
  mutate(cens = 1 - (cens - 1))

## produce K-M plot (assuming that interval-censored individuals
## plot non-sex specific KM
temp_fit <- surv_fit(Surv(tD, cens) ~ sex, data = quails_km_dat)
km_plot <- ggsurvplot(temp_fit, data = quails_km_dat, conf.int = TRUE)
km_plot

## set and sample some random seeds for different models
set.seed(42)
seeds <- round(runif(5, 0, 100000000))

## create output folder
#dir.create("outputs")

## save data
rm(temp_fit, tL, tU, quails, tB)
save.image("quails.RData")