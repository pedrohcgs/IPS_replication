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
MC_sims <- foreach::foreach(nn = 1:nrep, 
                            .errorhandling = c("remove"),
                            .options.snow = opts
                      ) %dorng% {
  #----------------------------------------------------------------------------
  #Generate the data
  datamatrix <- dgps_ips(n, dgp)
  treat <- datamatrix$Treat
  xpscore1 <- base::as.matrix(datamatrix[,(3:6)])
  trueX <- base::as.matrix(datamatrix[,(7:10)])
  truePS <- datamatrix[,11]
  y <- datamatrix$Y
  #----------------------------------------------------------------------------
  # fit different estimators for the propensity score
  #Logit GLM estimator
  fit.glm <- stats::glm(treat ~ xpscore1 , family = binomial, x=T)
  xcov <- fit.glm$x
  #Logit CBPS estimator -- just-identification
  fit.cbps <- CBPS::CBPS(treat ~ xpscore1, ATT = 0,
                         twostep = T, standardize = F, method = "exact")
  
  #Logit CBPS estimator -- over-identification
  fit.cbps2 <- CBPS::CBPS(treat ~ xpscore1, ATT = 0,
                          twostep = T, standardize = F, method = "over")
  
  # IPS with indicator
  fit.ind <- IPS::IPS_ind(d = treat, x = xcov,
                          beta.initial = fit.cbps2$coefficients, maxit = 25000)
  
  #IPS with exponential
  fit.exp <- IPS::IPS_exp(d = treat, x = xcov,
                          beta.initial = fit.cbps2$coefficients, maxit = 25000)
  
  #IPS with projection
  fit.proj <- IPS::IPS_proj(d = treat, x = xcov,
                            beta.initial = fit.cbps2$coefficients, maxit = 25000)
  #----------------------------------------------------------------------------
  # compute ATE using different estimators
  # IPS with indicator
  ate_IPS_ind <- IPS::ATE(y, treat, xcov, 
                          fit.ind$fitted.values, fit.ind$lin.rep)

  # IPS with exponential
  ate_IPS_exp <- IPS::ATE(y, treat, xcov, 
                          fit.exp$fitted.values, fit.exp$lin.rep)
  
  # IPS with exponential
  ate_IPS_proj <- IPS::ATE(y, treat, xcov,
                           fit.proj$fitted.values, fit.proj$lin.rep)
  
    # GLM
  glm.lin.rep <- inflc_glm(treat, xcov, fit.glm$fitted.values)
  ate_glm <- IPS::ATE(y, treat, xcov, fit.glm$fitted.values, glm.lin.rep)
  
  
  # CBPS -- just-identification
  cbps.lin.rep <- inflc_cbps(treat, xcov,
                             fit.cbps$fitted.values, method = "exact")
  ate_cbps <- IPS::ATE(y, treat, xcov, 
                       fit.cbps$fitted.values, cbps.lin.rep)
  
  # CBPS -- over-identification
  cbps.lin.rep2<- inflc_cbps(treat, xcov, 
                             fit.cbps2$fitted.values, method = "over")
  ate_cbps2 <- IPS::ATE(y, treat, xcov, 
                        fit.cbps2$fitted.values, cbps.lin.rep2)
  
  #----------------------------------------------------------------------------
  # compute QTE using different estimators
  tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  # Bandwidth 
  bw = "nrd0" 
  # IPS with indicator
  qte_IPS_ind <- IPS::QTE(y, treat, xcov, fit.ind$fitted.values,
                          fit.ind$lin.rep, tau, bw = bw)
  
  # IPS with exponential
  qte_IPS_exp <- IPS::QTE(y, treat, xcov, fit.exp$fitted.values,
                          fit.exp$lin.rep, tau, bw = bw)
  
  # IPS with exponential
  qte_IPS_proj <- IPS::QTE(y, treat, xcov, fit.proj$fitted.values,
                           fit.proj$lin.rep, tau, bw = bw)
  
  # GLM
  qte_glm <- IPS::QTE(y, treat, xcov, fit.glm$fitted.values, 
                      glm.lin.rep, tau, bw = bw)
  # CBPS -- just-identification
  qte_cbps <- IPS::QTE(y, treat, xcov, fit.cbps$fitted.values, 
                       cbps.lin.rep, tau, bw = bw)
  
  # CBPS -- over-identification
  qte_cbps2 <- IPS::QTE(y, treat, xcov, fit.cbps2$fitted.values, 
                        cbps.lin.rep2, tau, bw = bw)
  
  #----------------------------------------------------------------------------
  # Balance checks
  #GLM
  w_glm1 <- treat/fit.glm$fitted.values
  w_glm0 <- (1-treat)/(1 - fit.glm$fitted.values)
  w_glm1 <- w_glm1/mean(w_glm1)
  w_glm0 <- w_glm0/mean(w_glm0)
  w_glm_ate <- w_glm1 - w_glm0
  
  #CBps - just identified
  w_cbps1 <- treat/fit.cbps$fitted.values
  w_cbps0 <- (1-treat)/(1 - fit.cbps$fitted.values)
  w_cbps1 <- w_cbps1/mean(w_cbps1)
  w_cbps0 <- w_cbps0/mean(w_cbps0)
  w_cbps_ate <- w_cbps1 - w_cbps0
  
  
  #CBPS-over identified
  w_cbps21 <- treat/fit.cbps2$fitted.values
  w_cbps20 <- (1-treat)/(1 - fit.cbps2$fitted.values)
  w_cbps21 <- w_cbps21/mean(w_cbps21)
  w_cbps20 <- w_cbps20/mean(w_cbps20)
  w_cbps2_ate <- w_cbps21 - w_cbps20
  
  
  #IPS-exp
  w_exp1 <- treat/fit.exp$fitted.values
  w_exp0 <- (1-treat)/(1 - fit.exp$fitted.values)
  w_exp1 <- w_exp1/mean(w_exp1)
  w_exp0 <- w_exp0/mean(w_exp0)
  w_exp_ate <- w_exp1 - w_exp0
  
  
  #IPS-Proj
  w_proj1 <- treat/fit.proj$fitted.values
  w_proj0 <- (1-treat)/(1 - fit.proj$fitted.values)
  w_proj1 <- w_proj1/mean(w_proj1)
  w_proj0 <- w_proj0/mean(w_proj0)
  w_proj_ate <- w_proj1 - w_proj0
  
  #IPS-Ind
  w_ind1 <- treat/fit.ind$fitted.values
  w_ind0 <- (1-treat)/(1 - fit.ind$fitted.values)
  w_ind1 <- w_ind1/mean(w_ind1)
  w_ind0 <- w_ind0/mean(w_ind0)
  w_ind_ate <- w_ind1 - w_ind0
  
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
  
  
  cdf_balance_exp_1 <- abs(base::colMeans((w_exp1 - 1) * indicators))
  cdf_balance_ind_1 <- abs(base::colMeans((w_ind1 - 1) * indicators))
  cdf_balance_proj_1 <- abs(base::colMeans((w_proj1 - 1) * indicators))
  cdf_balance_cbps_1 <- abs(base::colMeans((w_cbps1 - 1) * indicators))
  cdf_balance_cbps2_1 <- abs(base::colMeans((w_cbps21 - 1) * indicators))
  cdf_balance_glm_1 <- abs(base::colMeans((w_glm1 - 1) * indicators))
  
  
  
  cdf_balance_exp_0 <- abs(base::colMeans((w_exp0 - 1) * indicators))
  cdf_balance_ind_0 <- abs(base::colMeans((w_ind0 - 1) * indicators))
  cdf_balance_proj_0 <- abs(base::colMeans((w_proj0 - 1) * indicators))
  cdf_balance_cbps_0 <- abs(base::colMeans((w_cbps0 - 1) * indicators))
  cdf_balance_cbps2_0 <- abs(base::colMeans((w_cbps20 - 1) * indicators))
  cdf_balance_glm_0 <- abs(base::colMeans((w_glm0 - 1) * indicators))
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
  
  error_glm <- mean((fit.glm$fitted.values - truePS)^2)
  error_cbps <- mean((fit.cbps$fitted.values - truePS)^2)
  error_cbps2 <- mean((fit.cbps2$fitted.values - truePS)^2)
  error_IPS_ind <- mean((fit.ind$fitted.values - truePS)^2)
  error_IPS_exp <- mean((fit.exp$fitted.values - truePS)^2)
  error_IPS_proj <- mean((fit.proj$fitted.values - truePS)^2)
  
  error_ps <- c(error_IPS_exp, error_IPS_ind, error_IPS_proj,
                error_cbps, error_cbps2, error_glm
                )
  #----------------------------------------------------------------------------
  #Return point estimates and standard errors and balance cehcks
  out <- c(
    # Point estimates
    ate_IPS_exp$ate, ate_IPS_ind$ate, ate_IPS_proj$ate, ate_cbps$ate, ate_cbps2$ate, ate_glm$ate,
    
    qte_IPS_exp$qte[1], qte_IPS_ind$qte[1], qte_IPS_proj$qte[1], qte_cbps$qte[1], qte_cbps2$qte[1], qte_glm$qte[1],
    qte_IPS_exp$qte[2], qte_IPS_ind$qte[2], qte_IPS_proj$qte[2], qte_cbps$qte[2], qte_cbps2$qte[2], qte_glm$qte[2],
    qte_IPS_exp$qte[3], qte_IPS_ind$qte[3], qte_IPS_proj$qte[3], qte_cbps$qte[3], qte_cbps2$qte[3], qte_glm$qte[3],
    qte_IPS_exp$qte[4], qte_IPS_ind$qte[4], qte_IPS_proj$qte[4], qte_cbps$qte[4], qte_cbps2$qte[4], qte_glm$qte[4],
    qte_IPS_exp$qte[5], qte_IPS_ind$qte[5], qte_IPS_proj$qte[5], qte_cbps$qte[5], qte_cbps2$qte[5], qte_glm$qte[5],

    # Stadandard errors
    ate_IPS_exp$ate.se, ate_IPS_ind$ate.se, ate_IPS_proj$ate.se, ate_cbps$ate.se, ate_cbps2$ate.se, ate_glm$ate.se,
    
    qte_IPS_exp$qte.se[1], qte_IPS_ind$qte.se[1], qte_IPS_proj$qte.se[1], qte_cbps$qte.se[1], qte_cbps2$qte.se[1], qte_glm$qte.se[1],
    qte_IPS_exp$qte.se[2], qte_IPS_ind$qte.se[2], qte_IPS_proj$qte.se[2], qte_cbps$qte.se[2], qte_cbps2$qte.se[2], qte_glm$qte.se[2],
    qte_IPS_exp$qte.se[3], qte_IPS_ind$qte.se[3], qte_IPS_proj$qte.se[3], qte_cbps$qte.se[3], qte_cbps2$qte.se[3], qte_glm$qte.se[3],
    qte_IPS_exp$qte.se[4], qte_IPS_ind$qte.se[4], qte_IPS_proj$qte.se[4], qte_cbps$qte.se[4], qte_cbps2$qte.se[4], qte_glm$qte.se[4],
    qte_IPS_exp$qte.se[5], qte_IPS_ind$qte.se[5], qte_IPS_proj$qte.se[5], qte_cbps$qte.se[5], qte_cbps2$qte.se[5], qte_glm$qte.se[5],
    
    # Balancing test
    ks_tests, cvm_tests, error_ps
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

true_effects  <- c(rep(10, 36))

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

ps_error <- mean.mc[109:114]
#-----------------------------------------------------------------------------
#Export output: Inference
source(here("Codes/MC_tables_inf.R"))
write.csv(out_sim_inf, file = here("Simulations/Stylized/Results/Exog/", 
                                   paste0("MC_exog_inf.csv")))
#Export output: point estimate
source(here("Codes", "MC_tables_point.R"))
write.csv(out_sim_point, file = here("Simulations/Stylized/Results/Exog/", 
                                     paste0("MC_exog_point.csv")))

#Export output: balance
source(here("Codes", "MC_tables_bal.R"))
write.csv(out_sim_bal, file = here("Simulations/Stylized/Results/Exog/", 
                                   paste0("MC_exog_bal.csv")))
#----------------------------------------------------------------------------
# Save Dataset
save.image(here("Simulations/Stylized/Results/Exog/",
                paste0("stylized_sim_exog_n", n,"_dgp_",dgp,".RData"))
)

  