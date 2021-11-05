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
ncores  <- 36                  # Number of cores to use in parallel
seed1   <- 07232021            # Set initial seed (guaranteed reproducibility)
nrep <- 1000              # Monte Carlo replications
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

out_reg_d0 <- lm(OR_itt_formula, data = df[df$e401==0,])
out_reg_d1 <- lm(OR_itt_formula, data = df[df$e401==1,])

X_or <- model.matrix(OR_itt_formula, data = df)
# get conditional means
df$out_reg_d0 <- X_or %*% ( out_reg_d0$coefficients)
df$out_reg_d1 <- X_or %*% ( out_reg_d1$coefficients)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Generate potential outcomes 
Y0 <- df$out_reg_d0/1000 + rnorm(n_total)
Y1 <- df$out_reg_d1/1000 + rnorm(n_total)

df$sorting <- (Y1 - Y0)

#----------------------------------------------------------------------------
# propensity score Specifications
ps_formula <- as.formula("p401 ~ sorting")


ps_treat <- glm(formula = ps_formula,
                family  = binomial(link = "logit"),
                data    = df[df$e401==1,],
                x = TRUE)

ps_treat <- ps_treat$coefficients
#ps_treat <- c(0.811, 0.055)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Keep only necessary portions of data
df <- data.frame(
  out_reg_d0 = df$out_reg_d0,
  out_reg_d1 = df$out_reg_d1,
  instrument_ps_baseline = df$instrument_ps,
  covariates
)

# df <- df[order(df$inc),]
#----------------------------------------------------------------------------
# True ATE
ate_true <- 0.8106185
# True QTE (based on simulations) - linear OR
# qte_true <- c(-474.4316 , 3420.1549 , 9453.3591,  9544.2945, 15063.9459 )
# NOnlinear OR
qte_true <- c( -0.1476849 , 0.2006515 , 0.7719787 , 1.3601936,  2.0244570  )



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
    df_sim <- df[(id_sim),]
    
    #xcov <- df_sim[,-(1:3)]
    xcov  <- model.matrix(as.formula("~ (inc + loginc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira) +
                                      I(inc^2) + I(loginc^2) + I(age^2) +
                                      I(fsize^2) + I(educ^2)"), data = df_sim)
    
    
    # Y0 <- df_sim$out_reg_d0 +  stats::rnorm(n)
    # Y1 <- df_sim$out_reg_d1 +  stats::rnorm(n)
    # 
    Y0 <-   df_sim$out_reg_d0/1000 + stats::rnorm(n)
    Y1 <-   df_sim$out_reg_d1/1000 + stats::rnorm(n)
    Z <- stats::rbinom(n, 1, df_sim$instrument_ps)
    
    # Generate pscore and Treatment Status
    endog_index <- cbind(1, (Y1-Y0)) %*% ps_treat
    ps <- 1/(1+exp(- endog_index))
    D1 <- stats::rbinom(n, size=1, prob = ps)
    D0 <- rep(0,n)
    treat <- Z*D1 + (1-Z) * D0
    
    Y <- Y0 * (1 - treat) + Y1 * treat
    
    if(sum(treat) < 20) next
    if(sum(Z) < 50) next
    #----------------------------------------------------------------------------
    # fit different estimators for the propensity score
    #Logit GLM estimator
    fit.glm <- glm.fit(x = xcov, y = Z,
                       family  = binomial(link = "logit"))
    
   
    # GLM
    glm.lin.rep <- inflc_glm(Z, xcov, fit.glm$fitted.values)
    late_glm <- IPS::LATE(Y,Z, treat, xcov, fit.glm$fitted.values, glm.lin.rep)
    
    
    #----------------------------------------------------------------------------
    #----------------------------------------------------------------------------
    # compute QTE using different estimators
    tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    # Bandwidth 
    bw = "nrd0"
    
    # GLM
    lqte_glm <- IPS::LQTE(Y,  Z, treat, xcov, fit.glm$fitted.values, 
                          glm.lin.rep, tau, bw = bw)
    
    
    
    #----------------------------------------------------------------------------
    #Return point estimates and standard errors 
    out <- c(
      # Point estimates
      late_glm$late,
      
      lqte_glm$lqte,
      
      # std. errors
      late_glm$late.se,
      
      lqte_glm$lqte.se
    )
    
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

true_effects  <- c(rep(ate_true, 1),
                   rep(qte_true[1], 1),
                   rep(qte_true[2], 1),
                   rep(qte_true[3], 1),
                   rep(qte_true[4], 1),
                   rep(qte_true[5], 1)
)
true_effects <- matrix(true_effects, nrow = nrow(mc), ncol = 6, byrow = TRUE)

bias_mc <- mean.mc[1:6] - colMeans(true_effects)

RMSE_mc <- sqrt(base::colMeans((mc[,1:6] - true_effects)^2, na.rm = TRUE ))

abs_bias_mc <- (base::colMeans(abs(mc[,1:6] - true_effects), na.rm = TRUE))


coverage <-  colMeans(((mc[,1:6] - 1.96*mc[,7:12]) <= true_effects) * 
                        ((mc[,1:6] + 1.96*mc[,7:12]) >= true_effects), na.rm = TRUE)



