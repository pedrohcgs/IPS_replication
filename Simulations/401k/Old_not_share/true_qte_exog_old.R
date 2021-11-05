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
ncores  <- 8                  # Number of cores to use in parallel
seed1   <- 07232021            # Set initial seed (guaranteed reproducibility)
nrep <- 1000                  # Monte Carlo replications
n <- 1e6
set.seed(seed1)
kappa_out <- 2500
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
library(glm2)
library(np)
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
# File to run all ATE/ATTs dending on the PS estimation method.
source(here("Codes", "average_itt.R"))
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
cts_variables_names <- c("inc", "loginc", "age", "fsize", "educ",
                         "incsq", "logincsq", "agesq", "fsizesq", "educsq")

binary_variables_names <- c("hown","marr", "twoearn", "db", "pira")
covariates_names <- c(cts_variables_names,binary_variables_names )
outcomes_names <- c("tw", "net_tfa")
treat_names <- c("e401", "p401")

all_variables <- c(outcomes_names,
                   treat_names,
                   covariates_names)

df <- data401k[,colnames(data401k) %in% all_variables]
#----------------------------------------------------------------------------
# Sample size of original data
n_total <- dim(df)[1]
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Generate potential outcomes
# Outcome regressions for intention to treat

OR_itt_formula <- as.formula("net_tfa ~ (inc + loginc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira)")

# OR_itt_formula <- as.formula("net_tfa ~ (inc + loginc + age +fsize + educ+ hown +
#                                       marr+ twoearn+ db+ pira)^2 +
#                                       incsq + logincsq + agesq + fsizesq +
#                                       educsq")

out_reg_z0 <- lm(OR_itt_formula, data = df[df$e401==0,])
out_reg_z1 <- lm(OR_itt_formula, data = df[df$e401==1,])

X_or <- model.matrix(OR_itt_formula, data = df)
# get conditional means
df$out_reg_z0 <- X_or %*% ( out_reg_z0$coefficients)
df$out_reg_z1 <- X_or %*% ( out_reg_z1$coefficients)


# Compute single index conditional density of errors to serve as dgp
het_error_z0 <- glm(out_reg_z0$residuals^2 ~ -1 + subset(X_or,df$e401==0), 
                    family = "poisson")

het_error_z1 <- glm(out_reg_z1$residuals^2 ~ -1 + subset(X_or,df$e401==1), 
                    family = "poisson")

# get conditional standard deviations
df$het_error_z0 <- exp(X_or %*% ( het_error_z0$coefficients)/2)
df$het_error_z1 <- exp(X_or %*% ( het_error_z1$coefficients)/2)


r_error_het_z0 <- new_r(out_reg_z0$residuals/sqrt(het_error_z0$fitted.values),
                        type = "continuous")
r_error_het_z1 <- new_r(out_reg_z1$residuals/sqrt(het_error_z1$fitted.values),
                        type = "continuous")
#----------------------------------------------------------------------------
# Keep only necessary portions of data
df <- data.frame(
  out_reg_z0 = df$out_reg_z0,
  out_reg_z1 = df$out_reg_z1,
  het_error_z0 = df$het_error_z0,
  het_error_z1 = df$het_error_z1
)

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Now, we have all parameters for the simulation
# Let's start
cl <- parallel::makeCluster(ncores)
registerDoSNOW(cl)

# Progress bar
pb <- txtProgressBar(max = nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb,n)
opts <- list(progress=progress)

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#Start the MONTE CARLO loop
MC_sims <- foreach(nn = 1:nrep, .options.snow = opts) %dorng%
  {
    # Draw covariates (which implies we draw instrument pscores, and OR, too)
    id_sim <- sample(1:n_total, n, replace = TRUE )
    df_sim <- df[id_sim,]
    
    
    Y0 <- df_sim$out_reg_z0 + df_sim$het_error_z0 * r_error_het_z0(n)
    Y1 <- df_sim$out_reg_z1 + df_sim$het_error_z1 * r_error_het_z1(n)
    
    
    #----------------------------------------------------------------------------
    # compute QTE using different estimators
    tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
   #
    true_qte <- stats::quantile(Y1, probs = tau) - stats::quantile(Y0, probs = tau)
    
    #----------------------------------------------------------------------------
    #Return point estimates
    out <- true_qte
    
    return(out)
  }
#-----------------------------------------------------------------------------
#Stop the cluster
stopCluster(cl)
#-----------------------------------------------------------------------------
#Put the Monte Carlo Results in an matrix
mc <- do.call(rbind, MC_sims) 
mc_all <- mc
mc <- na.omit(mc)
#-----------------------------------------------------------------------------
# Mean in the Monte Carlo
mean.mc <- base::colMeans(mc)

base::colMeans(mc^2 ) - mean.mc^2
