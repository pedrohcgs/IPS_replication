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
ncores  <- 48                  # Number of cores to use in parallel
seed1   <- 07232021            # Set initial seed (guaranteed reproducibility)
nrep <- 200                  # Monte Carlo replications
n <- 200
set.seed(seed1)
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
# instrument_ps_formula <- as.formula("e401 ~ (inc + loginc + age +fsize + educ+ hown +
#                                       marr+ twoearn+ db+ pira)")

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
# Outcome regressions for intention to treat
# OR_itt_formula <- as.formula("net_tfa ~ (inc + loginc + age +fsize + educ+ hown +
#                                       marr+ twoearn+ db+ pira)")

OR_itt_formula <- as.formula("net_tfa ~ (inc + loginc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira) +
                                      I(inc^2) + I(loginc^2) + I(age^2) +
                                      I(fsize^2) + I(educ^2)")

out_reg_z0 <- lm(OR_itt_formula, data = df[df$e401==0,])
out_reg_z1 <- lm(OR_itt_formula, data = df[df$e401==1,])

X_or <- model.matrix(OR_itt_formula, data = df)
# get conditional means
df$out_reg_z0 <- X_or %*% ( out_reg_z0$coefficients)
df$out_reg_z1 <- X_or %*% ( out_reg_z1$coefficients)
#----------------------------------------------------------------------------
# True ATE
ate_true <- mean(df$out_reg_z1 - df$out_reg_z0)/1000

# True QTE (based on simulations) - linear OR
# qte_true <- c(-474.4316 , 3420.1549 , 9453.3591,  9544.2945, 15063.9459 )
# NOnlinear OR with interactions
#qte_true <- c( 1015.525,  4096.162,  7554.571, 12798.330, 14560.076 )
# Nonlinear OR without interactions
qte_true <- c( -4.1862382,  0.6219272,  8.7278857, 11.3987574, 20.6641778 )

#----------------------------------------------------------------------------
# Keep only necessary portions of data
df <- data.frame(
  out_reg_z0 = df$out_reg_z0,
  out_reg_z1 = df$out_reg_z1,
  instrument_ps_baseline = df$instrument_ps,
  covariates
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
MC_sims <- foreach(nn = 1:nrep, 
                   .options.snow = opts, 
                   .errorhandling='remove') %dorng%
  {
    # Draw covariates (which implies we draw instrument pscores, and OR, too)
    id_sim <- sample(1:n_total, n, replace = TRUE )
    df_sim <- df[id_sim,]
    
    #xcov <- df_sim[,-(1:3)]
    xcov  <- model.matrix(as.formula("~ (inc + loginc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira) +
                                      I(inc^2) + I(loginc^2) + I(age^2) +
                                      I(fsize^2) + I(educ^2)"), data = df_sim)
    
    
    Y0 <- df_sim$out_reg_z0/1000 +  stats::rnorm(n)
    Y1 <- df_sim$out_reg_z1/1000 +  stats::rnorm(n)
    
    treat <- rbinom(n, 1, df_sim$instrument_ps)
    
    Y <- Y0 * (1 - treat) + Y1 * treat
    
    #----------------------------------------------------------------------------
    # fit different estimators for the propensity score
    #Logit GLM estimator
    fit.glm <- (glm.fit(x = xcov, y = treat,
                       family  = binomial(link = "logit")))
    
    #Logit CBPS estimator - just identified
    fit.cbps <- (CBPS::CBPS(treat ~ -1 + as.matrix(xcov), ATT = 0,
                           twostep = T, standardize = F,
                           method = "exact"))
    
    #Logit CBPS estimator - over identified
    fit.cbps2 <- (CBPS::CBPS(treat ~ -1 + as.matrix(xcov),
                            ATT = 0,
                            twostep = T, standardize = F, data = data401k,
                            method = "over"))
    
    start_beta <-  fit.glm$coefficients
    start_beta[is.nan(start_beta)] <- 0
    start_beta[is.na(start_beta)] <- 0

    # IPS with indicator
    fit.ind <- (IPS::IPS_ind(treat, xcov, beta.initial = start_beta, maxit = 5000))
    
    #IPS with exponential
    fit.exp <- (IPS::IPS_exp(treat, xcov, beta.initial = start_beta, maxit = 5000))
    
    #IPS with projection
    fit.proj <- (IPS::IPS_proj(treat, xcov, beta.initial =start_beta, maxit = 5000))
    
    
    #----------------------------------------------------------------------------
    # compute ATE using different estimators
    # IPS with indicator
    ate_IPS_ind <- (IPS::ATE(Y, treat, xcov, fit.ind$fitted.values, fit.ind$lin.rep))
    
    # IPS with exponential
    ate_IPS_exp <- (IPS::ATE(Y, treat, xcov, fit.exp$fitted.values, fit.exp$lin.rep))
    
    # IPS with exponential
    ate_IPS_proj <- (IPS::ATE(Y, treat, xcov, fit.proj$fitted.values, fit.proj$lin.rep))
    
    # GLM
    glm.lin.rep <- (inflc_glm(treat, xcov, fit.glm$fitted.values))
    ate_glm <- (IPS::ATE(Y, treat, xcov, fit.glm$fitted.values, glm.lin.rep))
    
    # CBPS -- just-identification
    cbps.lin.rep <- (inflc_cbps(treat, xcov, fit.cbps$fitted.values, method = "exact"))
    ate_cbps <- (IPS::ATE(Y, treat, xcov, fit.cbps$fitted.values, cbps.lin.rep))
    
    # CBPS -- over-identification
    cbps.lin.rep2<- (inflc_cbps(treat, xcov, fit.cbps2$fitted.values, method = "over"))
    ate_cbps2 <- (IPS::ATE(Y, treat, xcov, fit.cbps2$fitted.values, cbps.lin.rep2))
    #----------------------------------------------------------------------------
    #----------------------------------------------------------------------------
    # compute QTE using different estimators
    tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    # Bandwidth 
    bw = "nrd0" #"nrd"
    # IPS with indicator
    qte_IPS_ind <- (IPS::QTE(Y, treat, xcov, fit.ind$fitted.values,
                            fit.ind$lin.rep, tau, bw = bw))
    
    # IPS with exponential
    qte_IPS_exp <- (IPS::QTE(Y, treat, xcov, fit.exp$fitted.values,
                            fit.exp$lin.rep, tau, bw = bw))
    
    # IPS with exponential
    qte_IPS_proj <- (IPS::QTE(Y, treat, xcov, fit.proj$fitted.values,
                             fit.proj$lin.rep, tau, bw = bw))
    
    # GLM
    qte_glm <- (IPS::QTE(Y, treat, xcov, fit.glm$fitted.values, 
                        glm.lin.rep, tau, bw = bw))
    
    # CBPS -- just-identification
    qte_cbps <- (IPS::QTE(Y, treat, xcov, fit.cbps$fitted.values, 
                         cbps.lin.rep, tau, bw = bw))
    
    # CBPS -- over-identification
    qte_cbps2 <- (IPS::QTE(Y, treat, xcov, fit.cbps2$fitted.values, 
                          cbps.lin.rep2, tau, bw = bw))
    
    #----------------------------------------------------------------------------
    #Return point estimates and standard errors 
    out <- (cbind(
      # Point estimates
      ate_IPS_exp$ate, ate_IPS_ind$ate, ate_IPS_proj$ate, ate_cbps$ate, ate_cbps2$ate, ate_glm$ate,
      
      qte_IPS_exp$qte[1], qte_IPS_ind$qte[1], qte_IPS_proj$qte[1], qte_cbps$qte[1], qte_cbps2$qte[1], qte_glm$qte[1],
      qte_IPS_exp$qte[2], qte_IPS_ind$qte[2], qte_IPS_proj$qte[2], qte_cbps$qte[2], qte_cbps2$qte[2], qte_glm$qte[2],
      qte_IPS_exp$qte[3], qte_IPS_ind$qte[3], qte_IPS_proj$qte[3], qte_cbps$qte[3], qte_cbps2$qte[3], qte_glm$qte[3],
      qte_IPS_exp$qte[4], qte_IPS_ind$qte[4], qte_IPS_proj$qte[4], qte_cbps$qte[4], qte_cbps2$qte[4], qte_glm$qte[4],
      qte_IPS_exp$qte[5], qte_IPS_ind$qte[5], qte_IPS_proj$qte[5], qte_cbps$qte[5], qte_cbps2$qte[5], qte_glm$qte[5],
      
      # std. errors
      ate_IPS_exp$ate.se, ate_IPS_ind$ate.se, ate_IPS_proj$ate.se, ate_cbps$ate.se, ate_cbps2$ate.se, ate_glm$ate.se,
      
      qte_IPS_exp$qte.se[1], qte_IPS_ind$qte.se[1], qte_IPS_proj$qte.se[1], qte_cbps$qte.se[1], qte_cbps2$qte.se[1], qte_glm$qte.se[1],
      qte_IPS_exp$qte.se[2], qte_IPS_ind$qte.se[2], qte_IPS_proj$qte.se[2], qte_cbps$qte.se[2], qte_cbps2$qte.se[2], qte_glm$qte.se[2],
      qte_IPS_exp$qte.se[3], qte_IPS_ind$qte.se[3], qte_IPS_proj$qte.se[3], qte_cbps$qte.se[3], qte_cbps2$qte.se[3], qte_glm$qte.se[3],
      qte_IPS_exp$qte.se[4], qte_IPS_ind$qte.se[4], qte_IPS_proj$qte.se[4], qte_cbps$qte.se[4], qte_cbps2$qte.se[4], qte_glm$qte.se[4],
      qte_IPS_exp$qte.se[5], qte_IPS_ind$qte.se[5], qte_IPS_proj$qte.se[5], qte_cbps$qte.se[5], qte_cbps2$qte.se[5], qte_glm$qte.se[5]
    ))

    
    return(out)
  }
#-----------------------------------------------------------------------------
#Stop the cluster
stopCluster(cl)
#-----------------------------------------------------------------------------
#Put the Monte Carlo Results in an matrix
mc <- do.call(rbind, MC_sims) 
mc_all <- as.matrix(mc)
# mode(mc) = "numeric"
# mc2 <- as.matrix(MC_sims, nrow = nrep, ncol = 72)

#-----------------------------------------------------------------------------
# Mean in the Monte Carlo
mean.mc <- base::colMeans(mc, na.rm = TRUE)

true_effects  <- c(rep(ate_true, 6),
                   rep(qte_true[1], 6),
                   rep(qte_true[2], 6),
                   rep(qte_true[3], 6),
                   rep(qte_true[4], 6),
                   rep(qte_true[5], 6)
)
true_effects <- matrix(true_effects, nrow = nrow(mc), ncol = 36, byrow = TRUE)

bias_mc <- mean.mc[1:36] - colMeans(true_effects)

RMSE_mc <- sqrt(base::colMeans((mc[,1:36] - true_effects)^2, na.rm = TRUE ))

abs_bias_mc <- (base::colMeans(abs(mc[,1:36] - true_effects), na.rm = TRUE))


coverage <-  colMeans(((mc[,1:36] - 1.96*mc[,37:72]) <= true_effects) * 
                        ((mc[,1:36] + 1.96*mc[,37:72]) >= true_effects), na.rm = TRUE)

#----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#Export output: Inference
source(here("Codes/MC_tables_inf_exog.R"))
write.csv(out_sim_inf, file = here("Simulations/401k/", 
                                   paste0("MC_exog3_inf_n", n, ".csv")))
#Export output: point estimate
source(here("Codes", "MC_tables_point_exog.R"))
write.csv(out_sim_point, file = here("Simulations/401k/", 
                                     paste0("MC_exog3_point_n", n, ".csv")))
#----------------------------------------------------------------------------
# Save Dataset
save.image(here("Simulations/401k/",
                paste0("401k_sim_exog3_n", n, ".RData"))
)

