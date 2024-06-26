## load libraries
library(tidyverse)
library(boot)
library(lamW)
library(nimble)
rm(list = ls())

## source necessary R distributions
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")
source("Distributions/Dist_GompertzMakeham.R")
source("Distributions/Dist_GompertzMakehamNim.R")
source("Distributions/Dist_Siler.R")
source("Distributions/Dist_SilerNim.R")
source("Distributions/Dist_Expo.R")

## source additional R functions
source("ModelComparison_FUNCTIONS.R")

## load data
load("quails.RData")

## set seed according to model
#set.seed(seeds[4])

## set up plot output file
#pdf("outputs/ModelComparisons.pdf")

###########################################################
##                                                      ###
##          Now conduct model comparisons               ###
##                                                      ###
###########################################################

## load IS samples
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
MargLike.plot(mods)

## which models within log(20) of best
logmarg <- map_dbl(mods, "logmarg")
bestind <- which(logmarg == max(logmarg))
logmargLCI <- mods[[bestind]]$LCI
logmarg <- map_dbl(mods, "UCI")
logmarg <- logmarg[map_lgl(logmarg, ~ . >= logmargLCI - log(20))]
logmarg

## turn graphics device off
dev.off()
