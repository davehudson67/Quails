## load libraries
library(nimble)
library(tidyverse)
library(mvtnorm)
library(boot)
library(lamW)
library(GGally)
library(coda)
library(mclust)
library(parallel)
library(survminer)
library(survival)
library(MCMCvis)
library(data.table)
rm(list = ls())

## source necessary R distributions
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

## source additional R functions
source("ModelComparison_FUNCTIONS.R")

## load data
load("quails.RData")
quails <- read_csv("Survival.FINAL.csv")
sex <- as.numeric(as.factor(quails$sex)) - 1

## set seed according to model
#set.seed(seeds[3])

## set up plot output file
#pdf("outputs/SilerSexSpecificA1.pdf")

###########################################################
##                                                      ###
##   Now fit gompertz model with sex differences        ###
##                                                      ###
###########################################################

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated Siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dgompzNim(amult[i], bmult[i])
    
    log(amult[i]) <- log(a) + betaSEX[1] * sex[i] * zsex[1]
    log(bmult[i]) <- log(b) + betaSEX[2] * sex[i] * zsex[2]
  }
  
  ## priors
  for (j in 1:2){
    betaSEX[j] ~ dnorm(0, 1)
    zsex[j] ~ dbern(0.5)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  
})

## set up other components of model
consts <- list(nind = nrow(cint))

data <- list(cint = cint, censored = censored, tD = tD, sex = sex)

## find overdispersed initial values
tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + 1
    } else {
      y <- runif(1, x[1], x[2])
    }
    y
  })
}
initFn <- function(cint, censored, sex) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
    if(any(pars < 0)) {
      return(NA)
    }
    ll <- sum(dGompertz(t, a = pars[1], b = pars[2], log = TRUE))
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)

    ## optimise to interval-censored only
    pars <- optim(rexp(2, 100), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  ## check log-likelihoods
  ll <- sum(dGompertz(tD, a = pars[1], b = pars[2], log = TRUE))
  stopifnot(is.finite(ll))
  ## output initial values
  list(
    tD = tD,
    a = pars[1],
    b = pars[2],
    betaSEX = rnorm(2, 0, 1),
    zsex = rep(0, times = 2)
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored, sex))
#model$initializeInfo()

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = c("a", "b"), thin = 1, enableWAIC = TRUE)
config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
config$addSampler(target = c("a1", "a2", "b1", "c1"), type = 'AF_slice', control = 50)
#config$addSampler(target = c("a1", "c1"), type = 'AF_slice', control = 20)
config$addSampler(target = c("b2"), type = 'slice')
#help(samplers)

## load in custom RJ-MCMC samplers
#source("MCMC_RJ_multi.R")

## Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("betaSEX"),
            indicatorNodes = c("zsex"),
            control = list(mean = 0, scale = 0.5))

config$addMonitors("betaSEX", "zsex")
config

#Check monitors and samplers
config$printMonitors()
config$printSamplers(c("a", "b"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 500000, 
                           nburnin = 50000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1,
                           WAIC = TRUE))
run$summary
MCMCsummary(run$samples)

saveRDS(run, "IS_RJMCMC_NewModelCodeComparisons/samples_RJMCMC_MongSexLONGRUNnew.rds")
#samples <- readRDS("outputs/samples_RJMCMC_MongSex.rds")
run$WAIC

#Plot mcmcm
samples <- as.matrix(run$samples)
plot(run$samples)


zNames <- model$expandNodeNames('zsex')
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]

res

## Plot BMA predictions from model
t_pred <- seq(0, max(cint), length.out = 100)

predsM <- apply(samples, 1, function(pars, t, sex) {
  a <- exp(log(pars[1]) + pars[3] * sex * pars[5])
  b <- exp(log(pars[2]) + pars[4] * sex * pars[6])
  pGompertz(t, a = a, b = b, lower.tail = FALSE)
}, t = t_pred, sex = 1)

predsF <- apply(samples, 1, function(pars, t, sex) {
  a <- exp(log(pars[1]) + pars[3] * sex * pars[5])
  b <- exp(log(pars[2]) + pars[4] * sex * pars[6])
  pGompertz(t, a = a, b = b, lower.tail = FALSE)
}, t = t_pred, sex = 0)


SurvF <- predsF %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Female_BMA")

SurvM <- predsM %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Male_BMA")

Surv_BMA_log20models <- rbind(SurvF, SurvM)

## plot posterior predictive urvival curves
pred_plot <- ggplot(Surv_BMA_log20models, aes(x = t)) +
  geom_line(aes(y = Median, colour = category)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = category), alpha = 0.3) +
  geom_line(aes(x = t, y = Median), data = SurvF, alpha = 0.7) +
  geom_line(aes(x = t, y = LCI), data = SurvF, linetype = "dashed", alpha = 0.7) +
  geom_line(aes(x = t, y = UCI), data = SurvF, linetype = "dashed", alpha = 0.7) +
  geom_line(aes(x = t, y = Median), data = SurvM, alpha = 0.7) +
  geom_line(aes(x = t, y = LCI), data = SurvM, linetype = "dotted", alpha = 0.7) +
  geom_line(aes(x = t, y = UCI), data = SurvM, linetype = "dotted", alpha = 0.7) +
  xlab("Age (yrs)") + ylab("Survival probability") + labs(colour = "Sex", fill = "Sex") + 
  labs(title = "a.") +
  scale_fill_discrete(labels = c("Female", "Male")) +
  scale_colour_discrete(labels = c("Female", "Male")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) #+


pred_plot


