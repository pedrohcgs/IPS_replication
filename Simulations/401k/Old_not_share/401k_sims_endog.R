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
nrep <- 1000                   # Monte Carlo replications
n <- 1000
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
# Instrument propensity score Specifications
instrument_ps_formula <- as.formula("e401 ~ (inc + loginc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira)")

# instrument_ps_formula <- as.formula("e401 ~ (inc + loginc + age +fsize + educ+ hown +
#                                       marr+ twoearn+ db+ pira)^2 +
#                                       incsq + logincsq + agesq + fsizesq +
#                                       educsq")
instrument_ps <- glm(formula = instrument_ps_formula,
                     family  = binomial(link = "logit"),
                     data    = df,
                     x = TRUE)

X_prime_Beta_inst_ps <- instrument_ps$x %*% ( instrument_ps$coefficients)
df$instrument_ps <- 1/(1+exp(- X_prime_Beta_inst_ps))
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Generate potential outcomes
# Outcome regressions
OR_ate_formula <- as.formula("net_tfa ~ (inc + loginc + age +fsize + educ+ 
                                      hown + marr+ twoearn + db + pira)")

out_reg_d0 <- lm(OR_ate_formula, data = df[df$p401==0,])
out_reg_d1 <- lm(OR_ate_formula, data = df[df$p401==1,])

X_or <- model.matrix(OR_ate_formula, data = df)
# get conditional means
df$out_reg_d0 <- X_or %*% ( out_reg_d0$coefficients)
df$out_reg_d1 <- X_or %*% ( out_reg_d1$coefficients)
#----------------------------------------------------------------------------
# Compute single index conditional density of errors to serve as dgp
het_error_d0 <- suppressWarnings(
  glm(out_reg_d0$residuals^2 ~ -1 + subset(X_or,df$p401==0), 
      family = "poisson")
)

het_error_d1 <- suppressWarnings(
  glm(out_reg_d1$residuals^2 ~ -1 + subset(X_or,df$p401==1), 
      family = "poisson")
)


# get conditional standard deviations
df$het_error_d0 <- exp(X_or %*% ( het_error_d0$coefficients)/2)
df$het_error_d1 <- exp(X_or %*% ( het_error_d1$coefficients)/2)


r_error_het_d0 <- new_r(het_error_d0$residuals/sqrt(het_error_d0$fitted.values),
                        type = "continuous")
r_error_het_d1 <- new_r(het_error_d1$residuals/sqrt(het_error_d1$fitted.values),
                           type = "continuous")
#----------------------------------------------------------------------------
# Generate potential outcomes 
Y0 <- df$out_reg_d0 +  df$het_error_d0 * r_error_het_d0(n_total)
Y1 <- df$out_reg_d1 +  df$het_error_d1 * r_error_het_d1(n_total)

df$sorting <- (Y1 - Y0)/100000

#----------------------------------------------------------------------------
# propensity score Specifications
ps_formula <- as.formula("p401 ~ sorting")


ps_treat <- glm(formula = ps_formula,
                     family  = binomial(link = "logit"),
                     data    = df[df$e401==1,],
                     x = TRUE)

ps_treat <- ps_treat$coefficients
#----------------------------------------------------------------------------
# True ATT/LATE
ate_true <- 12762.719

# True QTE (based on simulations)
qte_true <- c(1806.972, 6971.817, 14855.657, 16438.111, 23080.030)
  
  
#----------------------------------------------------------------------------
# Keep only necessary portions of data
df <- data.frame(
  out_reg_d0 = df$out_reg_d0,
  out_reg_d1 = df$out_reg_d1,
  het_error_d0 = df$het_error_d0,
  het_error_d1 = df$het_error_d1,
  instrument_ps_baseline = df$instrument_ps,
  instrument_ps$x
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
MC_sims <- foreach(nn = 1:nrep, .options.snow = opts,
                   .errorhandling = c("remove")) %dorng%
  {
    # Draw covariates (which implies we draw instrument pscores, and OR, too)
    id_sim <- sample(1:n_total, n, replace = TRUE )
    df_sim <- df[id_sim,]
    
    xcov <- df_sim[,-(1:5)]
    
    Y0 <- df_sim$out_reg_d0 + df_sim$het_error_d0 * r_error_het_d0(n)
    Y1 <- df_sim$out_reg_d1 + df_sim$het_error_d1 * r_error_het_d1(n)
    
    Z <- stats::rbinom(n, 1, df_sim$instrument_ps)
    
    # Generate pscore and Treatment Status
    endog_index <- cbind(1,(Y1-Y0)/10000) %*% ps_treat
    ps <- 1/(1+exp(- endog_index))
    D1 <- stats::rbinom(n, size=1, prob = ps)
    D0 <- rep(0,n)
    treat <- Z*D1 + (1-Z) * D0

    Y <- Y0 * (1 - treat) + Y1 * treat
    
    #----------------------------------------------------------------------------
    # fit different estimators for the propensity score
    #Logit GLM estimator
    fit.glm <- glm.fit(x = xcov, y = Z,
                       family  = binomial(link = "logit"))
    
    #Logit CBPS estimator - just identified
    fit.cbps <- CBPS::CBPS(Z ~ -1 + as.matrix(xcov), ATT = 0,
                           twostep = T, standardize = F,
                           method = "exact")
    
    #Logit CBPS estimator - over identified
    fit.cbps2 <- CBPS::CBPS(Z ~ -1 + as.matrix(xcov),
                            ATT = 0,
                            twostep = T, standardize = F, data = data401k,
                            method = "over")
    # IPS with indicator
    fit.ind <- IPS::LIPS_ind(Z, treat, xcov, beta.initial = fit.cbps2$coefficients, maxit = 5000)
    
    #IPS with exponential
    fit.exp <- IPS::LIPS_exp(Z, treat, xcov, beta.initial = fit.cbps2$coefficients, maxit = 5000)
    
    #IPS with projection
    fit.proj <- IPS::LIPS_proj(Z, treat, xcov, beta.initial = fit.cbps2$coefficients, maxit = 5000)
    
    
    #----------------------------------------------------------------------------
    # compute ATE using different estimators
    # IPS with indicator
    late_IPS_ind <- IPS::LATE(Y, Z, treat, xcov, fit.ind$fitted.values, fit.ind$lin.rep)
    
    # IPS with exponential
    late_IPS_exp <- IPS::LATE(Y,  Z, treat, xcov, fit.exp$fitted.values, fit.exp$lin.rep)
    
    # IPS with exponential
    late_IPS_proj <- IPS::LATE(Y, Z,  treat, xcov, fit.proj$fitted.values, fit.proj$lin.rep)
    
    # GLM
    glm.lin.rep <- inflc_glm(Z, xcov, fit.glm$fitted.values)
    late_glm <- IPS::LATE(Y,  Z, treat, xcov, fit.glm$fitted.values, glm.lin.rep)
    
    # CBPS -- just-identification
    cbps.lin.rep <- inflc_cbps(Z, xcov, fit.cbps$fitted.values, method = "exact")
    late_cbps <- IPS::LATE(Y,  Z, treat, xcov, fit.cbps$fitted.values, cbps.lin.rep)
    
    # CBPS -- over-identification
    cbps.lin.rep2<- inflc_cbps(Z, xcov, fit.cbps2$fitted.values, method = "over")
    late_cbps2 <- IPS::LATE(Y,  Z, treat, xcov, fit.cbps2$fitted.values, cbps.lin.rep2)
    #----------------------------------------------------------------------------
    #----------------------------------------------------------------------------
    # compute QTE using different estimators
    tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    # Bandwidth 
    bw = "nrd0"
    # IPS with indicator
    lqte_IPS_ind <- IPS::LQTE(Y,  Z, treat, xcov, fit.ind$fitted.values,
                            fit.ind$lin.rep, tau, bw = bw)
    
    # IPS with exponential
    lqte_IPS_exp <- IPS::LQTE(Y,  Z, treat, xcov, fit.exp$fitted.values,
                            fit.exp$lin.rep, tau, bw = bw)
    
    # IPS with exponential
    lqte_IPS_proj <- IPS::LQTE(Y, Z,  treat, xcov, fit.proj$fitted.values,
                             fit.proj$lin.rep, tau, bw = bw)
    
    # GLM
    lqte_glm <- IPS::LQTE(Y,  Z, treat, xcov, fit.glm$fitted.values, 
                        glm.lin.rep, tau, bw = bw)
    
    # CBPS -- just-identification
    lqte_cbps <- IPS::LQTE(Y, Z,  treat, xcov, fit.cbps$fitted.values, 
                         cbps.lin.rep, tau, bw = bw)
    
    # CBPS -- over-identification
    lqte_cbps2 <- IPS::LQTE(Y, Z,  treat, xcov, fit.cbps2$fitted.values, 
                          cbps.lin.rep2, tau, bw = bw)
    
    #----------------------------------------------------------------------------
    #Return point estimates and standard errors 
    out <- cbind(
      # Point estimates
      late_IPS_exp$late, late_IPS_ind$late, late_IPS_proj$late, late_cbps$late, late_cbps2$late, late_glm$late,
      
      lqte_IPS_exp$lqte[1], lqte_IPS_ind$lqte[1], lqte_IPS_proj$lqte[1], lqte_cbps$lqte[1], lqte_cbps2$lqte[1], lqte_glm$lqte[1],
      lqte_IPS_exp$lqte[2], lqte_IPS_ind$lqte[2], lqte_IPS_proj$lqte[2], lqte_cbps$lqte[2], lqte_cbps2$lqte[2], lqte_glm$lqte[2],
      lqte_IPS_exp$lqte[3], lqte_IPS_ind$lqte[3], lqte_IPS_proj$lqte[3], lqte_cbps$lqte[3], lqte_cbps2$lqte[3], lqte_glm$lqte[3],
      lqte_IPS_exp$lqte[4], lqte_IPS_ind$lqte[4], lqte_IPS_proj$lqte[4], lqte_cbps$lqte[4], lqte_cbps2$lqte[4], lqte_glm$lqte[4],
      lqte_IPS_exp$lqte[5], lqte_IPS_ind$lqte[5], lqte_IPS_proj$lqte[5], lqte_cbps$lqte[5], lqte_cbps2$lqte[5], lqte_glm$lqte[5],
      
      # std. errors
      late_IPS_exp$late.se, late_IPS_ind$late.se, late_IPS_proj$late.se, late_cbps$late.se, late_cbps2$late.se, late_glm$late.se,
      
      lqte_IPS_exp$lqte.se[1], lqte_IPS_ind$lqte.se[1], lqte_IPS_proj$lqte.se[1], lqte_cbps$lqte.se[1], lqte_cbps2$lqte.se[1], lqte_glm$lqte.se[1],
      lqte_IPS_exp$lqte.se[2], lqte_IPS_ind$lqte.se[2], lqte_IPS_proj$lqte.se[2], lqte_cbps$lqte.se[2], lqte_cbps2$lqte.se[2], lqte_glm$lqte.se[2],
      lqte_IPS_exp$lqte.se[3], lqte_IPS_ind$lqte.se[3], lqte_IPS_proj$lqte.se[3], lqte_cbps$lqte.se[3], lqte_cbps2$lqte.se[3], lqte_glm$lqte.se[3],
      lqte_IPS_exp$lqte.se[4], lqte_IPS_ind$lqte.se[4], lqte_IPS_proj$lqte.se[4], lqte_cbps$lqte.se[4], lqte_cbps2$lqte.se[4], lqte_glm$lqte.se[4],
      lqte_IPS_exp$lqte.se[5], lqte_IPS_ind$lqte.se[5], lqte_IPS_proj$lqte.se[5], lqte_cbps$lqte.se[5], lqte_cbps2$lqte.se[5], lqte_glm$lqte.se[5]
    )
    
    return(out)
  }
#-----------------------------------------------------------------------------
#Stop the cluster
stopCluster(cl)
#-----------------------------------------------------------------------------
#Put the Monte Carlo Results in an matrix
mc <- do.call(rbind, MC_sims) 
#mc_all <- mc
#mc <- na.omit(mc)
#-----------------------------------------------------------------------------
# Mean in the Monte Carlo
mean.mc <- base::colMeans(mc)

true_effects  <- c(rep(ate_true, 6),
                   rep(qte_true[1], 6),
                   rep(qte_true[2], 6),
                   rep(qte_true[3], 6),
                   rep(qte_true[4], 6),
                   rep(qte_true[5], 6)
)
true_effects <- matrix(true_effects, nrow = nrep, ncol = 36, byrow = TRUE)

bias_mc <- mean.mc[1:36] - colMeans(true_effects)

RMSE_mc <- sqrt(base::colMeans((mc[,1:36]-true_effects)^2) )



abs_bias_mc <- (base::colMeans(abs(mc[,1:36] - true_effects)))


coverage <-  colMeans(((mc[,1:36] - 1.96*mc[,37:72]) <= true_effects) * 
                        ((mc[,1:36] + 1.96*mc[,37:72]) >= true_effects))

#----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#Export output: Inference
source(here("Codes/MC_tables_inf_exog.R"))
write.csv(out_sim_inf, file = here("Simulations/401k/", 
                                   paste0("MC_endog_inf_n", n, ".csv")))
#Export output: point estimate
source(here("Codes", "MC_tables_point_exog.R"))
write.csv(out_sim_point, file = here("Simulations/401k/", 
                                     paste0("MC_endog_point_n", n, ".csv")))
#----------------------------------------------------------------------------
# Save Dataset
save.image(here("Simulations/401k/",
                paste0("401k_sim_endog_n", n, ".RData"))
)

