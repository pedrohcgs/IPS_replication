################################################################
#
# This file reproduces the Monte Carlo tables of 
# "Integrated Propensity Score"
#
# All simulations were run in a Windows 10 Workstation, with 
#  R 4.1
# Setup with endogeneity
################################################################
################################################################
#---------------------------------------------------
#      Main file to replicate Empirically-Driven Simulations
#      Sant'Anna, Song and Xu (2021),
# Simulations based on R 4.1
# Updated: 02/06/2021
#---------------------------------------------------
#-----------------------------------------------------------------------------
# Startup - clear memory, load packages, and set parameters
# Clear memory
rm(list = ls())
#-----------------------------------------------------------------------------
# Basic parameters for the simulation - Doesn't change over setups
ncores  <- 40                  # Number of cores to use in parallel
seed1   <- 07232021            # Set initial seed (guaranteed reproducibility)
nrep <- 1000                   # Monte Carlo replications
trim <- FALSE                # Should trim extreme pscores?
trim.at <- 0.025             # If trim is TRUE, what is threshold - keep (trim.at, 1-trim.at)
set.seed(seed1)
kappa_out <- 1000
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# load the necessary libraries
library(here)
library(foreach)
library(doSNOW)
library(doRNG)
library(dplyr)
# MAKE SURE YOU HAVE INSTALLED THE IPS library
#devtools::install_github("pedrohcgs/IPS")
library(IPS)
library(CBPS)
library(estimatr)  
library(foreign)
#devtools::install_github("echasnovski/pdqr")
library(pdqr)
#-----------------------------------------------------------------------------
# Set seed
set.seed(seed1)
#-----------------------------------------------------------------------------
# Source Auxiliary functions
# Influence functions for CBPS and GLM-based Pscores
source(here("Codes", "Inflc_glm.R"))
source(here("Codes", "Inflc_CBPS.R"))
#----------------------------------------------------------------------------
###########################################################################
# Load data from Chernozhukov and Hansen (2004)
load(here("Applications/401k/data/data401k.RData"))
data401k <- data401k[data401k$inc>0,]
#----------------------------------------------------------------------------
# Rescale some variables
data401k$age <- (data401k$age -25)/(64-25)
data401k$inc <- (data401k$inc+2652)/(242124+2652)
data401k$fsize <- data401k$fsize/13
data401k$educ <- data401k$educ/18
#----------------------------------------------------------------------------
# Generate some covariate transformations
data401k$incsq <- data401k$inc^2
data401k$loginc <- log(data401k$inc)
data401k$logincsq <- data401k$loginc^2
data401k$agesq <- data401k$age^2
data401k$educsq <- data401k$educ^2 
data401k$fsizesq <- data401k$fsize^2 
#----------------------------------------------------------------------------
# Keep only variables we will use
cts_variables_names <- c("inc", "loginc", "age", "fsize", "educ")

binary_variables_names <- c("hown","marr", "twoearn", "db", "pira")
covariates_names <- c(cts_variables_names,binary_variables_names )
outcomes_names <- c("tw", "net_tfa")
treat_names <- c("e401", "p401")

all_variables <- c(outcomes_names,
                   treat_names,
                   covariates_names)

df <- data401k[,colnames(data401k) %in% all_variables]
covariates <- df[,covariates_names]
#----------------------------------------------------------------------------
# Sample size of original data
n_total <- dim(df)[1]
#----------------------------------------------------------------------------
# Instrument propensity score Specifications
instrument_ps_formula <- as.formula("e401 ~ (inc + loginc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira) +
                                      I(inc^2) + I(loginc^2) + I(age^2) +
                                      I(fsize^2) + I(educ^2)")

instrument_ps <- glm(formula = instrument_ps_formula,
                     family  = binomial(link = "logit"),
                     data    = df,
                     x = TRUE)

X_prime_Beta_inst_ps <- instrument_ps$x %*% ( instrument_ps$coefficients)
df$instrument_ps <- 1/(1+exp(- X_prime_Beta_inst_ps))
#----------------------------------------------------------------------------
# Generate potential outcomes
OR_itt_formula <- as.formula("net_tfa ~ (inc + loginc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira)^2 +
                                      I(inc^2) + I(loginc^2) + I(age^2) +
                                      I(fsize^2) + I(educ^2)")

out_reg_d0 <- lm(OR_itt_formula, data = df[df$e401==0,])
out_reg_d1 <- lm(OR_itt_formula, data = df[df$e401==1,])

X_or <- model.matrix(OR_itt_formula, data = df)
# get conditional means
df$out_reg_d0 <- X_or %*% ( out_reg_d0$coefficients)
df$out_reg_d1 <- X_or %*% ( out_reg_d1$coefficients)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Generate potential outcomes 
Y0 <- df$out_reg_d0/kappa_out + rnorm(n_total)
Y1 <- df$out_reg_d1/kappa_out + rnorm(n_total)

df$sorting <- (Y1 - Y0)

#----------------------------------------------------------------------------
# propensity score Specifications
ps_formula <- as.formula("p401 ~ sorting")


ps_treat <- glm(formula = ps_formula,
                family  = binomial(link = "logit"),
                data    = df[df$e401==1,],
                x = TRUE)

ps_treat <- ps_treat$coefficients
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Keep only necessary portions of data
df <- data.frame(
  out_reg_d0 = df$out_reg_d0,
  out_reg_d1 = df$out_reg_d1,
  instrument_ps_baseline = df$instrument_ps,
  covariates
)
#----------------------------------------------------------------------------
# True :ATE and LQTE
if(kappa_out == 1000){
  late_true <- 9.042082
  lqte_true <- c(1.369598, 4.182532 , 7.760703, 13.350003, 15.232624)
}

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
    source(here::here("Simulations/401k/401k_sims_lte.R"))
  }
}
#}
#----------------------------------------------------------------------------