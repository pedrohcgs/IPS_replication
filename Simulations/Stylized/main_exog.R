################################################################
#
# This file reproduces the Monte Carlo tables of 
# "Integrated Propensity Score"
#
# All simulations were run in a Windows 10 Workstation, with 
#  R 4.1
# Setup with exogeneous treatments
################################################################
#-----------------------------------------------------------------------------
# Startup - clear memory, load packages, and set parameters
# Clear memory
rm(list = ls())
#-----------------------------------------------------------------------------
# Basic parameters for the simulation - Doesn't change over setups
ncores  <- 40                  # Number of cores to use in parallel
seed1   <- 1234            # Set initial seed (guaranteed reproducibility)
nrep <- 1000                   # Monte Carlo replications
#----------------------------------------------------------------------------
# load the necessary libraries
library(here)
library(MASS)
library(doSNOW)
library(doRNG)
library(CBPS)
library(matrixStats)
#library(quantreg)
library(Matrix)
library(stats)
# Install IPS library
#devtools::install_github("pedrohcgs/IPS")
library(IPS)
#-----------------------------------------------------------------------------
# Set seed
set.seed(seed1)
#-----------------------------------------------------------------------------
# Source Auxiliary functions
source(here("Codes", "Inflc_glm.R"))
source(here("Codes", "Inflc_CBPS.R"))
source(here("Codes", "dgps_ips.R"))
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Run the simulations
  for (nn in 1:3){
    for (dgp in 1:2){
      
      #set sample size
      if(nn==1) n <- 200
      if(nn==2) n <- 500
      if(nn==3) n <- 1000
      
      #do the Monte Carlo and produce all tables
      source(here::here("Simulations/Stylized/sims_exog.R"))

    }
  }
#}
#----------------------------------------------------------------------------