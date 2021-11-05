###########################################################################
# Application: 401k - Endog setup
###########################################################################
# Startup - clear memory, load packages, and set parameters 
# Clear memory
rm(list = ls(all = T))
#----------------------------------------------------------------------------
# Set seed
set.seed(1234)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Install and load package from github
library(foreign)
library(MASS)
library(CBPS)
library(Rcpp)
library(Matrix)
library(stats)
# MAKE SURE YOU HAVE INSTALLED THE IPS library
#devtools::install_github("pedrohcgs/IPS")
library(IPS)
library(lmtest)
library(sandwich)
library(here)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Source the files with the Auxiliary R-functions
source(here("Codes", "Inflc_glm.R"))
source(here("Codes", "Inflc_CBPS.R"))
source(here("Codes", "average_late.R"))
source(here("Codes", "quantile_lte.R"))

#----------------------------------------------------------------------------
###########################################################################
###########################################################################
# Application: 401k
# Based on Benjamin (2003), Chernozhukov and Hansen (2004), and Wuthrich (2019)
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
#----------------------------------------------------------------------------
# propensity score Specifications
# Spec following Benjamin (2003)
inst.ps.baseline <- as.formula("e401 ~ (inc + loginc + age +fsize + educ+
                                          hown +
                                          marr+ twoearn+ db+ pira)^2 +
                                          incsq + logincsq + agesq + fsizesq + educsq")

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# LTE with baseline specification
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Instrument Propensity scores
#----------------------------------------------------------------------------
# Logit with GLM - baseline
ps.glm <- glm(formula = inst.ps.baseline,
              family  = binomial(link = "logit"),
              data    = data401k,
              x = T)
#----------------------------------------------------------------------------
# Matrix of covariate transformations and binary variable "treatment" indicator
baselineX <- stats::model.matrix(ps.glm)
treat <- data401k$p401
z <- data401k$e401

xbal_names <- c("inc", "loginc" , "age" ,"fsize" , "educ",
                "hown" ,
                "marr", "twoearn", "db", "pira")
xbal <- baselineX[,xbal_names]
#----------------------------------------------------------------------------
# Now, fit the balancing pscores
#----------------------------------------------------------------------------
#Logit CBPS estimator - just identified
ps.cbps <- CBPS::CBPS(inst.ps.baseline, ATT = 0,
                      twostep = T, standardize = F, data = data401k,
                      method = "exact")

#Logit CBPS estimator - over identified
ps.cbps.over <- CBPS::CBPS(inst.ps.baseline, ATT = 0,
                           twostep = T, standardize = F, data = data401k,
                           method = "over")
# IPS with indicator
s.ind <- Sys.time() 
ps.ind <- IPS::LIPS_ind(z = z, 
                        d = treat, 
                        x = baselineX, 
                        xbal = xbal)
e.ind <- Sys.time() 
e.ind - s.ind

#IPS with exponential
s.exp <- Sys.time() 
ps.exp <- IPS::LIPS_exp(z = z,
                        d = treat, 
                        x = baselineX)
e.exp <- Sys.time() 
e.exp - s.exp

#IPS with projection
s.proj <- Sys.time() 
ps.proj <- IPS::LIPS_proj(z = z, 
                          d = treat,
                          x = baselineX)
e.proj <- Sys.time() 
e.proj - s.proj
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# compute LATE using different PS estimators
#----------------------------------------------------------------------------
#outcome is total wealth or Net total Financial Assets
tw.late <- average_late(data401k$tw)
write.csv(tw.late,
          file =  here("Applications/401k/Results/", 
                       paste0("out_TotalWealth_late.csv"))
)
#outcome is Net total Financial Assets
net_tfa.late <- average_late(data401k$net_tfa)

write.csv(net_tfa.late,
          file =  here("Applications/401k/Results/", 
                       paste0("out_NetTotalFA_late.csv"))
)

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Estimate LQTE using different PS estimators
#----------------------------------------------------------------------------
# The tau for quantile effects
# CODE BELOW IS FOR these 5 quantile indexes only!!!
tau <- c(0.10, 0.25, 0.5, 0.75, 0.90)
#outcome is total wealth or Net total Financial Assets
tw.lqte <- quantile_lte(data401k$tw, tau)

write.csv(tw.lqte,
          file =  here("Applications/401k/Results/", 
                       paste0("out_TotalWealth_lqte.csv"))
)

#outcome is Net total Financial Assets
net_tfa.lqte <- quantile_lte(data401k$net_tfa, tau)
write.csv(net_tfa.lqte,
          file =  here("Applications/401k/Results/", 
                       paste0("out_NetTotalFA_lqte.csv"))
)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Distributional Balance
#----------------------------------------------------------------------------
xbal_names <- c("inc", "loginc" , "age" ,"fsize" , "educ",
                "hown" ,
                "marr", "twoearn", "db", "pira")
xbal <- baselineX[,xbal_names]
n <- dim(xbal)[1]
#----------------------------------------------------------------------------
# Balance checks
#GLM
w_glm1 <- treat * ((z/ps.glm$fitted.values) - ((1 - z)/(1 - ps.glm$fitted.values)))
w_glm0 <- (1 - treat) * ((z/ps.glm$fitted.values) - ((1 - z)/(1 - ps.glm$fitted.values)))
w_glm1 <- w_glm1/mean(w_glm1)
w_glm0 <- w_glm0/mean(w_glm0)
w_glm_ate <- w_glm1 - w_glm0
w_glm_ov <- 1 - (1 - treat) * (z/ps.glm$fitted.values) - 
  treat * (1 - z)/(1 - ps.glm$fitted.values)
w_glm_ov = w_glm_ov/mean(w_glm_ov)

#CBps - just identified
w_cbps1 <- treat * ((z/ps.cbps$fitted.values) - ((1 - z)/(1 - ps.cbps$fitted.values)))
w_cbps0 <- (1 - treat) * ((z/ps.cbps$fitted.values) - ((1 - z)/(1 - ps.cbps$fitted.values)))
w_cbps1 <- w_cbps1/mean(w_cbps1)
w_cbps0 <- w_cbps0/mean(w_cbps0)
w_cbps_ate <- w_cbps1 - w_cbps0
w_cbps_ov <- 1 - (1 - treat) * (z/ps.cbps$fitted.values) - 
  treat * (1 - z)/(1 - ps.cbps$fitted.values)
w_cbps_ov = w_cbps_ov/mean(w_cbps_ov)

#CBPS-over identified
w_cbps21 <- treat * ((z/ps.cbps.over$fitted.values) - ((1 - z)/(1 - ps.cbps.over$fitted.values)))
w_cbps20 <- (1 - treat) * ((z/ps.cbps.over$fitted.values) - ((1 - z)/(1 - ps.cbps.over$fitted.values)))
w_cbps21 <- w_cbps21/mean(w_cbps21)
w_cbps20 <- w_cbps20/mean(w_cbps20)
w_cbps2_ate <- w_cbps21 - w_cbps20
w_cbps2_ov <- 1 - (1 - treat) * (z/ps.cbps.over$fitted.values) - 
  treat * (1 - z)/(1 - ps.cbps.over$fitted.values)
w_cbps2_ov = w_cbps2_ov/mean(w_cbps2_ov)

#IPS-exp
w_exp1 <- treat * ((z/ps.exp$fitted.values) - ((1 - z)/(1 - ps.exp$fitted.values)))
w_exp0 <- (1 - treat) * ((z/ps.exp$fitted.values) - ((1 - z)/(1 - ps.exp$fitted.values)))
w_exp1 <- w_exp1/mean(w_exp1)
w_exp0 <- w_exp0/mean(w_exp0)
w_exp_ate <- w_exp1 - w_exp0
w_exp_ov <- 1 - (1 - treat) * (z/ps.exp$fitted.values) - 
  treat * (1 - z)/(1 - ps.exp$fitted.values)
w_exp_ov = w_exp_ov/mean(w_exp_ov)

#IPS-Proj
w_proj1 <- treat * ((z/ps.proj$fitted.values) - ((1 - z)/(1 - ps.proj$fitted.values)))
w_proj0 <- (1 - treat) * ((z/ps.proj$fitted.values) - ((1 - z)/(1 - ps.proj$fitted.values)))
w_proj1 <- w_proj1/mean(w_proj1)
w_proj0 <- w_proj0/mean(w_proj0)
w_proj_ate <- w_proj1 - w_proj0
w_proj_ov <- 1 - (1 - treat) * (z/ps.proj$fitted.values) - 
  treat * (1 - z)/(1 - ps.proj$fitted.values)
w_proj_ov = w_proj_ov/mean(w_proj_ov)
#IPS-Ind
w_ind1 <- treat * ((z/ps.ind$fitted.values) - ((1 - z)/(1 - ps.ind$fitted.values)))
w_ind0 <- (1 - treat) * ((z/ps.ind$fitted.values) - ((1 - z)/(1 - ps.ind$fitted.values)))
w_ind1 <- w_ind1/mean(w_ind1)
w_ind0 <- w_ind0/mean(w_ind0)
w_ind_ate <- w_ind1 - w_ind0
w_ind_ov <- 1 - (1 - treat) * (z/ps.ind$fitted.values) - 
  treat * (1 - z)/(1 - ps.ind$fitted.values)
w_ind_ov = w_ind_ov/mean(w_ind_ov)


k_dim = dim(xbal)[2]

indicators <- base::outer(xbal[,1], xbal[,1], "<=")
for ( bb in 2:k_dim){
  indicators <- indicators * 
    base::outer(xbal[, bb], xbal[, bb], "<=")
}

cdf_balance_exp <- abs(base::colMeans(w_exp_ate * indicators))
cdf_balance_ind <- abs(base::colMeans(w_ind_ate * indicators))
cdf_balance_proj <- abs(base::colMeans(w_proj_ate * indicators))
cdf_balance_cbps <- abs(base::colMeans(w_cbps_ate * indicators))
cdf_balance_cbps2 <- abs(base::colMeans(w_cbps2_ate * indicators))
cdf_balance_glm <- abs(base::colMeans(w_glm_ate * indicators))


cdf_balance_exp_1 <- abs(base::colMeans((w_exp1 - w_exp_ov) * indicators))
cdf_balance_ind_1 <- abs(base::colMeans((w_ind1 - w_ind_ov) * indicators))
cdf_balance_proj_1 <- abs(base::colMeans((w_proj1 - w_proj_ov) * indicators))
cdf_balance_cbps_1 <- abs(base::colMeans((w_cbps1 - w_cbps_ov) * indicators))
cdf_balance_cbps2_1 <- abs(base::colMeans((w_cbps21 - w_cbps2_ov) * indicators))
cdf_balance_glm_1 <- abs(base::colMeans((w_glm1 - w_glm_ov) * indicators))



cdf_balance_exp_0 <- abs(base::colMeans((w_exp0 - w_exp_ov) * indicators))
cdf_balance_ind_0 <- abs(base::colMeans((w_ind0 - w_ind_ov) * indicators))
cdf_balance_proj_0 <- abs(base::colMeans((w_proj0 - w_proj_ov) * indicators))
cdf_balance_cbps_0 <- abs(base::colMeans((w_cbps0 - w_cbps_ov) * indicators))
cdf_balance_cbps2_0 <- abs(base::colMeans((w_cbps20 - w_cbps2_ov) * indicators))
cdf_balance_glm_0 <- abs(base::colMeans((w_glm0 - w_glm_ov) * indicators))
#ks

ks_exp <-  sqrt(n) * c(
  max(cdf_balance_exp), 
  max(cdf_balance_exp_1), 
  max(cdf_balance_exp_0)
)

ks_ind <-  sqrt(n) * c(
  max(cdf_balance_ind), 
  max(cdf_balance_ind_1), 
  max(cdf_balance_ind_0)
)

ks_proj <-  sqrt(n) * c(
  max(cdf_balance_proj), 
  max(cdf_balance_proj_1), 
  max(cdf_balance_proj_0)
)

ks_cbps <-  sqrt(n) * c(
  max(cdf_balance_cbps), 
  max(cdf_balance_cbps_1), 
  max(cdf_balance_cbps_0)
)

ks_cbps2 <-  sqrt(n) * c(
  max(cdf_balance_cbps2), 
  max(cdf_balance_cbps2_1), 
  max(cdf_balance_cbps2_0)
)

ks_glm <-  sqrt(n) * c(
  max(cdf_balance_glm), 
  max(cdf_balance_glm_1), 
  max(cdf_balance_glm_0)
)


cvm_exp <-  c(
  sum(cdf_balance_exp^2), 
  sum(cdf_balance_exp_1^2), 
  sum(cdf_balance_exp_0^2)
)

cvm_ind <-  c(
  sum(cdf_balance_ind^2), 
  sum(cdf_balance_ind_1^2), 
  sum(cdf_balance_ind_0^2)
)

cvm_proj <-  c(
  sum(cdf_balance_proj^2), 
  sum(cdf_balance_proj_1^2), 
  sum(cdf_balance_proj_0^2)
)

cvm_cbps <-  c(
  sum(cdf_balance_cbps^2), 
  sum(cdf_balance_cbps_1^2), 
  sum(cdf_balance_cbps_0^2)
)

cvm_cbps2 <-  c(
  sum(cdf_balance_cbps2^2), 
  sum(cdf_balance_cbps2_1^2), 
  sum(cdf_balance_cbps2_0^2)
)

cvm_glm <-  c(
  sum(cdf_balance_glm^2), 
  sum(cdf_balance_glm_1^2), 
  sum(cdf_balance_glm_0^2)
)




#Generate output table

out_p_bal <- matrix(0, ncol = 6, nrow = 6)
out_p_bal <- data.frame(out_p_bal)

out_p_bal[1,] <- 100 * c(ks_exp/sqrt(n), sqrt(cvm_exp/n))
out_p_bal[2,] <- 100 * c(ks_ind/sqrt(n), sqrt(cvm_ind/n))
out_p_bal[3,] <- 100 * c(ks_proj/sqrt(n), sqrt(cvm_proj/n))
out_p_bal[4,] <- 100 * c(ks_cbps/sqrt(n), sqrt(cvm_cbps/n))
out_p_bal[5,] <- 100 * c(ks_cbps2/sqrt(n), sqrt(cvm_cbps2/n))
out_p_bal[6,] <- 100 * c(ks_glm/sqrt(n), sqrt(cvm_glm/n))

out_p_bal <- t(out_p_bal)

rownames(out_p_bal) <- c("ks", "ks_1", "ks_0", "cvm", "cvm_1", "cvm_0")
colnames(out_p_bal) <- c("Exp", "Ind", "Proj", "CBPS-just", "CBPS-over", "GLM") 


write.csv(out_p_bal,
          file =  here("Applications/401k/Results/", 
                       paste0("dist_balance_lte.csv"))
)

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Save Dataset
save.image( here("Applications/401k/Results/401k_lte.RData"))
#----------------------------------------------------------------------------
