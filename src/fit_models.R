##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    Fit Dynamic occupancy model
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())

library(stringr)
library(foreach)
library(tidyverse)
library(raster)

library(LaplacesDemon)
library(runjags)
library(rjags)
library(coda)

### 1. PARAMETERS FOR THE SIMULATION -------------------------------------------

# Years of interest
YEARMIN = 2006
YEARMAX = 2021

# Number of weeks per season
T=7

# min number of days sampled to consider week valid
Nobs_min = 3 

# min number of weeks sampled a year to consider the year valid
Nweek_min = 3 

# min number of photo a day to consider the observation valid
Npic_min = 35 

### 2. GET DATA AND UTILITY FUNCTIONS ------------------------------------------

source("src/format_data.R")
source("src/model_init.R")

### 3. RUN MODEL ---------------------------------------------------------------

## MCMC parameters ##

n.chains <- 4
adapt <- 1000
burnin <- 10000
sample <- 2500
thin <- 1

## Model specification ##


psi_covs <- c("int", "CLG", "TFG", "rodents_fall", "supp_feeding", "treatment")
gam_covs <- c("int", "CLG", "TFG", "rodents_fall", "supp_feeding", "treatment")
eps_covs <- c("int", "CLG", "TFG", "rodents_fall", "supp_feeding", "treatment")
pi_covs  <- c("int", "rodents_fall", "supp_feeding", "treatment")
tau_covs <- c("int", "rodents_fall", "supp_feeding", "treatment")
rho_covs <- c("int")

param_covs <- list(psi = psi_covs,
                   gam = gam_covs,
                   eps = eps_covs,
                   pi  = pi_covs,
                   tau = tau_covs,
                   rho = rho_covs)


data_list <- list(psi_cov = covs[,psi_covs]%>%as.matrix(),
                  gam_cov = covs[,gam_covs]%>%as.matrix(),
                  eps_cov = covs[,eps_covs]%>%as.matrix(),
                  pi_cov = covs[,pi_covs]%>%as.matrix(),
                  tau_cov = covs[,tau_covs]%>%as.matrix(),
                  rho_cov = covs[,rho_covs]%>%as.matrix(),
                  year_cov = covs[,"year"]%>%as.matrix(),
                  ncov_psi = length(psi_covs), ncov_gam = length(gam_covs),
                  ncov_eps = length(eps_covs), ncov_pi = length(pi_covs),
                  ncov_tau = length(tau_covs), ncov_rho = length(rho_covs), 
                  nyear = length(years),
                  nspec = 2, nseason = T, nsite = M,
                  nsurvey = K, nout = 4,
                  y = ob_state,
                  bait=bait)


M <- run.jags(model = "src/dcom.R",
              monitor = c("a", "b", "d","f","g","h", "bait_",
                          "yr_rho", "yr_psi", "yr_gam", "yr_eps",
                          "sigma_rho", "sigma_psi", "sigma_gam","sigma_eps",
                          "z"),
              data = data_list,
              n.chains = n.chains,
              inits = inits,
              adapt = adapt,
              burnin = burnin,
              sample = sample,
              thin = thin,
              summarise = TRUE,
              plots = FALSE,
              method = "parallel")

### 4. SAVE OUTPUT ---------------------------------------------------------------

saveRDS(M, "outputs/M_wtreat.rds")
