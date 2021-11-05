# Run the simulations for IPS paper
# Unified
####################################################################
####################################################################
#Make cluster
cl <- parallel::makeCluster(ncores) 
doSNOW::registerDoSNOW(cl)


# Progress bar
pb <- txtProgressBar(max = nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb,n)
opts <- list(progress=progress)
####################################################################
#Set seed
base::set.seed(seed1)
iseed <- floor(seed1 + n*2598 + dgp*711)
#base::set.seed(iseed)
####################################################################
#Start the MONTE CARLO loop
MC_sims <- foreach::foreach(nn = 1:nrep, .options.snow = opts, 
                            .errorhandling = c("remove")) %dorng% {
  #----------------------------------------------------------------------------
  #Generate the data
  datamatrix <- dgps_ips_late(n, dgp)
  treat <- datamatrix$Treat
  instru <- datamatrix$Instrument
  xpscore1 <- base::as.matrix(datamatrix[,(4:7)])
  trueX <- base::as.matrix(datamatrix[,-(1:7)])
  y <- datamatrix$Y
  
  if(sum(treat) < 20) next
  if(sum(instru) < 50) next
  #----------------------------------------------------------------------------
  # fit different estimators for the propensity score
  #Logit GLM estimator
  fit.glm <- stats::glm(instru ~ xpscore1 , family = binomial, x=TRUE)
  xcov <- fit.glm$x

   #Logit CBPS estimator -- just-identification
  fit.cbps <- CBPS::CBPS(instru ~ xpscore1, ATT = 0,
                         twostep = T, standardize = F, method = "exact")
  
  #Logit CBPS estimator -- over-identification
  fit.cbps2 <- CBPS::CBPS(instru ~ xpscore1, ATT = 0,
                          twostep = T, standardize = F, method = "over")
  
  # IPS with indicator
  fit.ind <-IPS::LIPS_ind(z = instru, d = treat, x = xcov,
                              beta.initial = fit.cbps2$coefficients,
                              maxit = 25000)
  
  #IPS with exponential
  fit.exp <- IPS::LIPS_exp(z = instru, d = treat, x = xcov,
                               beta.initial = fit.cbps2$coefficients,
                              maxit = 25000)
  #IPS with projection
  fit.proj <- IPS::LIPS_proj(z = instru, d = treat, x = xcov,
                                 beta.initial = fit.cbps2$coefficients,
                                  maxit = 25000)
  
  #----------------------------------------------------------------------------
  # compute ATE using different estimators
  # IPS with indicator
  late_IPS_ind <- IPS::LATE(y, instru, treat, xcov, 
                            fit.ind$fitted.values,
                            fit.ind$lin.rep,
                            trim = FALSE, trim.at = 0.005)
  
  # IPS with exponential
  late_IPS_exp <- IPS::LATE(y, instru, treat, xcov, fit.exp$fitted.values,
                            fit.exp$lin.rep,
                            trim = FALSE, trim.at = 0.005)
  
  # IPS with exponential
  late_IPS_proj <- IPS::LATE(y, instru, treat, xcov, fit.proj$fitted.values, 
                             fit.proj$lin.rep,
                             trim = FALSE, trim.at = 0.005)
  
  # GLM
  glm.lin.rep <- inflc_glm(instru, xcov, fit.glm$fitted.values)
  late_glm <- IPS::LATE(y, instru, treat, xcov, fit.glm$fitted.values, glm.lin.rep,
                        trim = F, trim.at = 0.005)
  
  
  # CBPS -- just-identification
  cbps.lin.rep <- inflc_cbps(instru, xcov, fit.cbps$fitted.values, method = "exact")
  late_cbps <- IPS::LATE(y, instru, treat, xcov, fit.cbps$fitted.values, cbps.lin.rep,
                         trim = F, trim.at = 0.005)
  
  # CBPS -- over-identification
  cbps.lin.rep2<- inflc_cbps(instru, xcov, fit.cbps2$fitted.values, method = "over")
  late_cbps2 <- IPS::LATE(y, instru, treat, xcov, fit.cbps2$fitted.values, cbps.lin.rep2,
                          trim = F, trim.at = 0.005)
  
  #----------------------------------------------------------------------------
  # compute QTE using different estimators
  tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  # Bandwidth 
  bw = "nrd0"
  # IPS with indicator
  lqte_IPS_ind <- IPS::LQTE(y, instru, treat, xcov, fit.ind$fitted.values,
                            fit.ind$lin.rep, tau,
                            trim = F, trim.at = 0.005,
                            bw = bw)
  
  # IPS with exponential
  lqte_IPS_exp <- IPS::LQTE(y, instru, treat, xcov, fit.exp$fitted.values,
                            fit.exp$lin.rep, tau,  
                            trim = F, trim.at = 0.005,
                            bw = bw)
  
  
  # IPS with exponential
  lqte_IPS_proj <- IPS::LQTE(y, instru, treat, xcov, fit.proj$fitted.values,
                             fit.proj$lin.rep, tau,  
                             trim = F, trim.at = 0.005,
                             bw = bw)
  
  
  # GLM
  lqte_glm <- IPS::LQTE(y, instru, treat, xcov, fit.glm$fitted.values, 
                        glm.lin.rep, tau, 
                        trim = F, trim.at = 0.005,
                        bw = bw)
  
  
  # CBPS -- just-identification
  lqte_cbps <- IPS::LQTE(y, instru, treat, xcov, fit.cbps$fitted.values, 
                         cbps.lin.rep, tau,  
                         trim = F, trim.at = 0.005,
                         bw = bw)
  
  
  # CBPS -- over-identification
  lqte_cbps2 <- IPS::LQTE(y, instru, treat, xcov, fit.cbps2$fitted.values, 
                          cbps.lin.rep2, tau, 
                          trim = F, trim.at = 0.005,
                          bw = bw)
  
  #----------------------------------------------------------------------------
  # Balance checks
  #GLM
  w_glm1 <- treat * ((instru/fit.glm$fitted.values) - ((1 - instru)/(1 - fit.glm$fitted.values)))
  w_glm0 <- (1 - treat) * ((instru/fit.glm$fitted.values) - ((1 - instru)/(1 - fit.glm$fitted.values)))
  w_glm1 <- w_glm1/mean(w_glm1)
  w_glm0 <- w_glm0/mean(w_glm0)
  w_glm_ate <- w_glm1 - w_glm0
  w_glm_ov <- 1 - (1 - treat) * (instru/fit.glm$fitted.values) - 
    treat * (1 - instru)/(1 - fit.glm$fitted.values)
  w_glm_ov = w_glm_ov/mean(w_glm_ov)
  
  #CBps - just identified
  w_cbps1 <- treat * ((instru/fit.cbps$fitted.values) - ((1 - instru)/(1 - fit.cbps$fitted.values)))
  w_cbps0 <- (1 - treat) * ((instru/fit.cbps$fitted.values) - ((1 - instru)/(1 - fit.cbps$fitted.values)))
  w_cbps1 <- w_cbps1/mean(w_cbps1)
  w_cbps0 <- w_cbps0/mean(w_cbps0)
  w_cbps_ate <- w_cbps1 - w_cbps0
  w_cbps_ov <- 1 - (1 - treat) * (instru/fit.cbps$fitted.values) - 
    treat * (1 - instru)/(1 - fit.cbps$fitted.values)
  w_cbps_ov = w_cbps_ov/mean(w_cbps_ov)
  
  #CBPS-over identified
  w_cbps21 <- treat * ((instru/fit.cbps2$fitted.values) - ((1 - instru)/(1 - fit.cbps2$fitted.values)))
  w_cbps20 <- (1 - treat) * ((instru/fit.cbps2$fitted.values) - ((1 - instru)/(1 - fit.cbps2$fitted.values)))
  w_cbps21 <- w_cbps21/mean(w_cbps21)
  w_cbps20 <- w_cbps20/mean(w_cbps20)
  w_cbps2_ate <- w_cbps21 - w_cbps20
  w_cbps2_ov <- 1 - (1 - treat) * (instru/fit.cbps2$fitted.values) - 
    treat * (1 - instru)/(1 - fit.cbps2$fitted.values)
  w_cbps2_ov = w_cbps2_ov/mean(w_cbps2_ov)
  
  #IPS-exp
  w_exp1 <- treat * ((instru/fit.exp$fitted.values) - ((1 - instru)/(1 - fit.exp$fitted.values)))
  w_exp0 <- (1 - treat) * ((instru/fit.exp$fitted.values) - ((1 - instru)/(1 - fit.exp$fitted.values)))
  w_exp1 <- w_exp1/mean(w_exp1)
  w_exp0 <- w_exp0/mean(w_exp0)
  w_exp_ate <- w_exp1 - w_exp0
  w_exp_ov <- 1 - (1 - treat) * (instru/fit.exp$fitted.values) - 
    treat * (1 - instru)/(1 - fit.exp$fitted.values)
  w_exp_ov = w_exp_ov/mean(w_exp_ov)
  
  #IPS-Proj
  w_proj1 <- treat * ((instru/fit.proj$fitted.values) - ((1 - instru)/(1 - fit.proj$fitted.values)))
  w_proj0 <- (1 - treat) * ((instru/fit.proj$fitted.values) - ((1 - instru)/(1 - fit.proj$fitted.values)))
  w_proj1 <- w_proj1/mean(w_proj1)
  w_proj0 <- w_proj0/mean(w_proj0)
  w_proj_ate <- w_proj1 - w_proj0
  w_proj_ov <- 1 - (1 - treat) * (instru/fit.proj$fitted.values) - 
    treat * (1 - instru)/(1 - fit.proj$fitted.values)
  w_proj_ov = w_proj_ov/mean(w_proj_ov)
  #IPS-Ind
  w_ind1 <- treat * ((instru/fit.ind$fitted.values) - ((1 - instru)/(1 - fit.ind$fitted.values)))
  w_ind0 <- (1 - treat) * ((instru/fit.ind$fitted.values) - ((1 - instru)/(1 - fit.ind$fitted.values)))
  w_ind1 <- w_ind1/mean(w_ind1)
  w_ind0 <- w_ind0/mean(w_ind0)
  w_ind_ate <- w_ind1 - w_ind0
  w_ind_ov <- 1 - (1 - treat) * (instru/fit.ind$fitted.values) - 
    treat * (1 - instru)/(1 - fit.ind$fitted.values)
  w_ind_ov = w_ind_ov/mean(w_ind_ov)
  
  #xcov = trueX
  k_dim = dim(xcov)[2]
  
  indicators <- base::outer(xcov[,1], xcov[,1], "<=")
  for ( bb in 2:k_dim){
    indicators <- indicators * 
      base::outer(xcov[, bb], xcov[, bb], "<=")
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
  ks_tests <- sqrt(n) * c(
    max(cdf_balance_exp), max(cdf_balance_ind), max(cdf_balance_proj),
    max(cdf_balance_cbps), max(cdf_balance_cbps2), max(cdf_balance_glm),
    
    
    max(cdf_balance_exp_1), max(cdf_balance_ind_1), max(cdf_balance_proj_1),
    max(cdf_balance_cbps_1), max(cdf_balance_cbps2_1), max(cdf_balance_glm_1),
    
    max(cdf_balance_exp_0), max(cdf_balance_ind_0), max(cdf_balance_proj_0),
    max(cdf_balance_cbps_0), max(cdf_balance_cbps2_0), max(cdf_balance_glm_0)
  )
  
  cvm_tests <- c(
    sum(cdf_balance_exp^2), sum(cdf_balance_ind^2), sum(cdf_balance_proj^2),
    sum(cdf_balance_cbps^2), sum(cdf_balance_cbps2^2), sum(cdf_balance_glm^2),
    
    sum(cdf_balance_exp_1^2), sum(cdf_balance_ind_1^2), sum(cdf_balance_proj_1^2),
    sum(cdf_balance_cbps_1^2), sum(cdf_balance_cbps2_1^2), sum(cdf_balance_glm_1^2),
    
    sum(cdf_balance_exp_0^2), sum(cdf_balance_ind_0^2), sum(cdf_balance_proj_0^2),
    sum(cdf_balance_cbps_0^2), sum(cdf_balance_cbps2_0^2), sum(cdf_balance_glm_0^2)
  )
  
  
  
  
  
  #----------------------------------------------------------------------------
  #Return point estimates and standard errors 
  out <- c(
    # Point estimates
    late_IPS_exp$late, late_IPS_ind$late, late_IPS_proj$late, late_cbps$late, late_cbps2$late, late_glm$late,
    
    lqte_IPS_exp$lqte[1], lqte_IPS_ind$lqte[1], lqte_IPS_proj$lqte[1], lqte_cbps$lqte[1], lqte_cbps2$lqte[1], lqte_glm$lqte[1],
    lqte_IPS_exp$lqte[2], lqte_IPS_ind$lqte[2], lqte_IPS_proj$lqte[2], lqte_cbps$lqte[2], lqte_cbps2$lqte[2], lqte_glm$lqte[2],
    lqte_IPS_exp$lqte[3], lqte_IPS_ind$lqte[3], lqte_IPS_proj$lqte[3], lqte_cbps$lqte[3], lqte_cbps2$lqte[3], lqte_glm$lqte[3],
    lqte_IPS_exp$lqte[4], lqte_IPS_ind$lqte[4], lqte_IPS_proj$lqte[4], lqte_cbps$lqte[4], lqte_cbps2$lqte[4], lqte_glm$lqte[4],
    lqte_IPS_exp$lqte[5], lqte_IPS_ind$lqte[5], lqte_IPS_proj$lqte[5], lqte_cbps$lqte[5], lqte_cbps2$lqte[5], lqte_glm$lqte[5],
    
    # Stadandard errors
    late_IPS_exp$late.se, late_IPS_ind$late.se, late_IPS_proj$late.se, late_cbps$late.se, late_cbps2$late.se, late_glm$late.se,
    
    lqte_IPS_exp$lqte.se[1], lqte_IPS_ind$lqte.se[1], lqte_IPS_proj$lqte.se[1], lqte_cbps$lqte.se[1], lqte_cbps2$lqte.se[1], lqte_glm$lqte.se[1],
    lqte_IPS_exp$lqte.se[2], lqte_IPS_ind$lqte.se[2], lqte_IPS_proj$lqte.se[2], lqte_cbps$lqte.se[2], lqte_cbps2$lqte.se[2], lqte_glm$lqte.se[2],
    lqte_IPS_exp$lqte.se[3], lqte_IPS_ind$lqte.se[3], lqte_IPS_proj$lqte.se[3], lqte_cbps$lqte.se[3], lqte_cbps2$lqte.se[3], lqte_glm$lqte.se[3],
    lqte_IPS_exp$lqte.se[4], lqte_IPS_ind$lqte.se[4], lqte_IPS_proj$lqte.se[4], lqte_cbps$lqte.se[4], lqte_cbps2$lqte.se[4], lqte_glm$lqte.se[4],
    lqte_IPS_exp$lqte.se[5], lqte_IPS_ind$lqte.se[5], lqte_IPS_proj$lqte.se[5], lqte_cbps$lqte.se[5], lqte_cbps2$lqte.se[5], lqte_glm$lqte.se[5],
    
    # Balancing test
    ks_tests, cvm_tests
  )
  out <- matrix(out, ncol = length(out), nrow = 1)
  return(out)
}

#----------------------------------------------------------------------------
####################################################################
#Stop the cluster
stopCluster(cl)


#-----------------------------------------------------------------------------
#Put the Monte Carlo Results in an matrix
mc <- do.call(rbind, MC_sims) 

#-----------------------------------------------------------------------------
# Mean in the Monte Carlo
mean.mc <- base::colMeans(mc, na.rm = TRUE)


true1.mean = 39.25974
true1.q10 = 42.93326
true1.q25 = 36.93488 
true1.q50 = 34.39781 
true1.q75 = 36.93470 
true1.q90 = 42.93315

true_effects <- c(rep(true1.mean,6),
            rep(true1.q10,6),
            rep(true1.q25,6),
            rep(true1.q50,6), 
            rep(true1.q75,6),
            rep(true1.q90,6))

true_effects <- matrix(true_effects, nrow = nrow(mc), ncol = 36, byrow = TRUE)

bias_mc <- mean.mc[1:36] - base::colMeans(true_effects)

RMSE_mc <- sqrt(base::colMeans((mc[,1:36]-true_effects)^2, na.rm = TRUE) )

abs_bias_mc <- (base::colMeans(abs(mc[,1:36] - true_effects), na.rm = TRUE))

critical_value <- stats::qnorm(0.975)
coverage <-  colMeans(((mc[,1:36] - critical_value*mc[,37:72]) <= true_effects) * 
                        ((mc[,1:36] + critical_value*mc[,37:72]) >= true_effects),
                      na.rm = TRUE)


ks_bal <- mean.mc[73:90]
cvm_bal <- mean.mc[91:108]

#-----------------------------------------------------------------------------
#Export output: Inference
source(here("Codes/MC_tables_inf.R"))
write.csv(out_sim_inf, file = here("Simulations/Stylized/Results/Endog/", 
                                   paste0("MC_endog_inf.csv")))
#Export output: point estimate
source(here("Codes", "MC_tables_point.R"))
write.csv(out_sim_point, file = here("Simulations/Stylized/Results/Endog/", 
                                     paste0("MC_endog_point.csv")))

#Export output: balance
source(here("Codes", "MC_tables_bal.R"))
write.csv(out_sim_bal, file = here("Simulations/Stylized/Results/Endog/", 
                                   paste0("MC_endog_bal.csv")))
#----------------------------------------------------------------------------
# Save Dataset
save.image(here("Simulations/Stylized/Results/Endog/",
                paste0("stylized_sim_endog_n", n,"_dgp_",dgp,".RData"))
)

