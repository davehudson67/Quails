library(tidyverse)
library(patchwork)
library(boot)

## source additional R functions
source("ModelComparison_FUNCTIONS.R")

## source necessary R distributions
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")
source("Distributions/Dist_GompertzMakeham.R")
source("Distributions/Dist_GompertzMakehamNim.R")
source("Distributions/Dist_Siler.R")
source("Distributions/Dist_SilerNim.R")
source("Distributions/Dist_Expo.R")

## read in data required
logimpweight_s <- readRDS("Outputs/logimpweight_s.rds")
logimpweight_gm <- readRDS("Outputs/logimpweight_gm.rds")
logimpweight_g <- readRDS("Outputs/logimpweight_g.rds")
logimpweight_e <- readRDS("Outputs/logimpweight_e.rds")

## generate log marginal likelihoods
logmarg_s <- log_sum_exp_marg(logimpweight_s)
logmarg_gm <- log_sum_exp_marg(logimpweight_gm)
logmarg_g <- log_sum_exp_marg(logimpweight_g)
logmarg_e <- log_sum_exp_marg(logimpweight_e)

## bootstrap samples
imp_boot_s <- BootsPlot(logimpweight_s, 5000)
imp_boot_gm <- BootsPlot(logimpweight_gm, 5000)
imp_boot_g <- BootsPlot(logimpweight_g, 5000)
imp_boot_e <- BootsPlot(logimpweight_e, 5000)

## add prior model weights
priorp <- 1/4
ps <- logmarg_s + log(priorp)
pgm <- logmarg_gm + log(priorp)
pg <- logmarg_g + log(priorp)
pe <- logmarg_e + log(priorp)

p <- c(ps, pgm, pg, pe)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  S = imp_boot_s, 
  GM = imp_boot_gm, 
  G = imp_boot_g,
  E = imp_boot_e
)
m <- MargLike.plot(mods) +
  theme_bw() +
  scale_x_discrete(labels = c("Exponential", "Gompertz", "Gompertz-Makeham", "Siler")) +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  labs(title = "a.")

m
## plot survival curve of best fitting model with KM plot 
## load samples
samples <- readRDS("Outputs/samples_g.rds")

## set up required info
## load data
quails <- read_csv("Survival.FINAL.csv")

## set up data
tU <- quails$time
tB <- rep(0, nrow(quails))
tL <- if_else(quails$time == 0, quails$time, quails$time - 1)

## check last seen against death
#table(tL >= tU)
#summary(mong$birth_day[tL >= tU])
#summary(mong$death_day[tL < tU])

## set death time for censored individuals
#tU[tU == 0] <- NA

## set up censoring vector (1 = interval, 2 = right)
quails$censored <- if_else(quails$status == 1, 1, 2)
#quails$censored <- if_else(quails$status == 1, 1, 0)
censored <- quails$censored

## summaries
stopifnot(all(tL[!is.na(tU)] <= tU[!is.na(tU)]))
stopifnot(all(tB <= tL)) #some individuals died on the same day they were born

## normalise to survival times
tU <- tU - tB
tL <- tL - tB

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

## set up zL/zU input matrix
zL <- tL
zU <- tU

## predictive plot against data
t_pred <- seq(0, max(cint), length.out = 100)
preds <- as.matrix(samples)[1:2000, ] %>%
  apply(1, function(pars, t) {
    pGompertz(t, a = pars[1], b = pars[2], lower.tail = FALSE)
  }, t = t_pred) %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, )

## plot ignoring sex differences
pred_plot <- km_plot[[1]] +
  geom_line(aes(x = t, y = Median), data = preds) +
  geom_ribbon(aes(x = t, y = Median, ymin = LCI, ymax = UCI), data = preds, alpha = 0.3)
p <- pred_plot +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  labs(title = "b.", x = "Time (days)", y = "Survivorship")
p

## Mortality curve
pred_surv <- as.matrix(samples)[1:2000, ] %>%
  apply(1, function(pars, t) {
    pGompertz(t, a = pars[1], b = pars[2], lower.tail = FALSE)
  }, t = t_pred)
  
pred_dens <- as.matrix(samples)[1:2000, ] %>%
  apply(1, function(pars, t) {
    dGompertz(t, a = pars[1], b = pars[2])
  }, t = t_pred)

pred_mort <- pred_dens/pred_surv  
  
pred_m <- pred_mort %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, )

predm_plot <- ggplot() +
  geom_line(aes(x = t, y = Median), data = pred_m) +
  geom_ribbon(aes(x = t, y = Median, ymin = LCI, ymax = UCI), data = pred_m, alpha = 0.3)
pm <- predm_plot +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  labs(title = "c.", x = "Time (days)", y = "Mortality")
pm

p/pm
