
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
    
    if(dgp == 1) {
      xcov  <- model.matrix(as.formula("~ (inc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira) +
                                      loginc + 
                                      I(inc^2) + I(loginc^2) + I(age^2) +
                                      I(fsize^2) + I(educ^2)"), data = df_sim)
    } 
    if(dgp == 2){
      xcov  <- model.matrix(as.formula("~ (inc + age +fsize + educ+ hown +
                                      marr+ twoearn+ db+ pira)"), data = df_sim)
    }
    
    
    Y0 <- df_sim$out_reg_d0/kappa_out +  stats::rnorm(n)
    Y1 <- df_sim$out_reg_d1/kappa_out +  stats::rnorm(n)
    
    Z <- stats::rbinom(n, 1, df_sim$instrument_ps)
    
    # Generate pscore and Treatment Status
    endog_index <- cbind(1,(Y1-Y0)) %*% ps_treat
    ps <- 1/(1+exp(- endog_index))
    D1 <- stats::rbinom(n, size=1, prob = ps)
    D0 <- rep(0,n)
    treat <- Z*D1 + (1-Z) * D0
    
    Y <- Y0 * (1 - treat) + Y1 * treat
    
    if(sum(treat) < 20) next
    if(sum(Z) < 50) next
    
    xbal <- as.data.frame(df_sim[,-c(1,2,3,5)])
    #----------------------------------------------------------------------------
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
    fit.ind <- IPS::LIPS_ind(z = Z, d = treat, x = xcov, xbal = xbal,
                             beta.initial = fit.cbps2$coefficients, maxit = 25000)
    
    #IPS with exponential
    fit.exp <- IPS::LIPS_exp(z = Z, d = treat, x = xcov,
                             beta.initial = fit.cbps2$coefficients, maxit = 25000)
    
    #IPS with projection
    fit.proj <- IPS::LIPS_proj(z = Z, d = treat, x = xcov,
                               beta.initial = fit.cbps2$coefficients, maxit = 25000)
    
    
    #----------------------------------------------------------------------------
    # compute ATE using different estimators
    # IPS with indicator
    late_IPS_ind <- IPS::LATE(Y, Z, treat, xcov, fit.ind$fitted.values, fit.ind$lin.rep, 
                              trim = trim, trim.at = trim.at)
    
    # IPS with exponential
    late_IPS_exp <- IPS::LATE(Y,  Z, treat, xcov, fit.exp$fitted.values, fit.exp$lin.rep, 
                              trim = trim, trim.at = trim.at)
    
    # IPS with exponential
    late_IPS_proj <- IPS::LATE(Y, Z,  treat, xcov, fit.proj$fitted.values, fit.proj$lin.rep,
                               trim = trim, trim.at = trim.at)
    
    # GLM
    glm.lin.rep <- inflc_glm(Z, xcov, fit.glm$fitted.values)
    late_glm <- IPS::LATE(Y,  Z, treat, xcov, fit.glm$fitted.values, glm.lin.rep, 
                          trim = trim, trim.at = trim.at)
    
    # CBPS -- just-identification
    cbps.lin.rep <- inflc_cbps(Z, xcov, fit.cbps$fitted.values, method = "exact")
    late_cbps <- IPS::LATE(Y,  Z, treat, xcov, fit.cbps$fitted.values, cbps.lin.rep,
                           trim = trim, trim.at = trim.at)
    
    # CBPS -- over-identification
    cbps.lin.rep2<- inflc_cbps(Z, xcov, fit.cbps2$fitted.values, method = "over")
    late_cbps2 <- IPS::LATE(Y,  Z, treat, xcov, fit.cbps2$fitted.values, cbps.lin.rep2,
                            trim = trim, trim.at = trim.at)
    #----------------------------------------------------------------------------
    #----------------------------------------------------------------------------
    # compute QTE using different estimators
    tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    # Bandwidth 
    bw = "nrd0"
    # IPS with indicator
    lqte_IPS_ind <- IPS::LQTE(Y,  Z, treat, xcov, fit.ind$fitted.values,
                              fit.ind$lin.rep, tau, bw = bw, trim = trim, trim.at = trim.at)
    
    # IPS with exponential
    lqte_IPS_exp <- IPS::LQTE(Y,  Z, treat, xcov, fit.exp$fitted.values,
                              fit.exp$lin.rep, tau, bw = bw,
                              trim = trim, trim.at = trim.at)
    
    # IPS with exponential
    lqte_IPS_proj <- IPS::LQTE(Y, Z,  treat, xcov, fit.proj$fitted.values,
                               fit.proj$lin.rep, tau, bw = bw,
                               trim = trim, trim.at = trim.at)
    
    # GLM
    lqte_glm <- IPS::LQTE(Y,  Z, treat, xcov, fit.glm$fitted.values, 
                          glm.lin.rep, tau, bw = bw, 
                          trim = trim, trim.at = trim.at)
    
    # CBPS -- just-identification
    lqte_cbps <- IPS::LQTE(Y, Z,  treat, xcov, fit.cbps$fitted.values, 
                           cbps.lin.rep, tau, bw = bw,
                           trim = trim, trim.at = trim.at)
    
    # CBPS -- over-identification
    lqte_cbps2 <- IPS::LQTE(Y, Z,  treat, xcov, fit.cbps2$fitted.values, 
                            cbps.lin.rep2, tau, bw = bw, 
                            trim = trim, trim.at = trim.at)
    #----------------------------------------------------------------------------
    # Balance checks
    #GLM
    w_glm1 <- treat * ((Z/fit.glm$fitted.values) - ((1 - Z)/(1 - fit.glm$fitted.values)))
    w_glm0 <- (1 - treat) * ((Z/fit.glm$fitted.values) - ((1 - Z)/(1 - fit.glm$fitted.values)))
    w_glm1 <- w_glm1/mean(w_glm1)
    w_glm0 <- w_glm0/mean(w_glm0)
    w_glm_ate <- w_glm1 - w_glm0
    w_glm_ov <- 1 - (1 - treat) * (Z/fit.glm$fitted.values) - 
      treat * (1 - Z)/(1 - fit.glm$fitted.values)
    w_glm_ov = w_glm_ov/mean(w_glm_ov)
    
    #CBps - just identified
    w_cbps1 <- treat * ((Z/fit.cbps$fitted.values) - ((1 - Z)/(1 - fit.cbps$fitted.values)))
    w_cbps0 <- (1 - treat) * ((Z/fit.cbps$fitted.values) - ((1 - Z)/(1 - fit.cbps$fitted.values)))
    w_cbps1 <- w_cbps1/mean(w_cbps1)
    w_cbps0 <- w_cbps0/mean(w_cbps0)
    w_cbps_ate <- w_cbps1 - w_cbps0
    w_cbps_ov <- 1 - (1 - treat) * (Z/fit.cbps$fitted.values) - 
      treat * (1 - Z)/(1 - fit.cbps$fitted.values)
    w_cbps_ov = w_cbps_ov/mean(w_cbps_ov)
    
    #CBPS-over identified
    w_cbps21 <- treat * ((Z/fit.cbps2$fitted.values) - ((1 - Z)/(1 - fit.cbps2$fitted.values)))
    w_cbps20 <- (1 - treat) * ((Z/fit.cbps2$fitted.values) - ((1 - Z)/(1 - fit.cbps2$fitted.values)))
    w_cbps21 <- w_cbps21/mean(w_cbps21)
    w_cbps20 <- w_cbps20/mean(w_cbps20)
    w_cbps2_ate <- w_cbps21 - w_cbps20
    w_cbps2_ov <- 1 - (1 - treat) * (Z/fit.cbps2$fitted.values) - 
      treat * (1 - Z)/(1 - fit.cbps2$fitted.values)
    w_cbps2_ov = w_cbps2_ov/mean(w_cbps2_ov)
    
    #IPS-exp
    w_exp1 <- treat * ((Z/fit.exp$fitted.values) - ((1 - Z)/(1 - fit.exp$fitted.values)))
    w_exp0 <- (1 - treat) * ((Z/fit.exp$fitted.values) - ((1 - Z)/(1 - fit.exp$fitted.values)))
    w_exp1 <- w_exp1/mean(w_exp1)
    w_exp0 <- w_exp0/mean(w_exp0)
    w_exp_ate <- w_exp1 - w_exp0
    w_exp_ov <- 1 - (1 - treat) * (Z/fit.exp$fitted.values) - 
      treat * (1 - Z)/(1 - fit.exp$fitted.values)
    w_exp_ov = w_exp_ov/mean(w_exp_ov)
    
    #IPS-Proj
    w_proj1 <- treat * ((Z/fit.proj$fitted.values) - ((1 - Z)/(1 - fit.proj$fitted.values)))
    w_proj0 <- (1 - treat) * ((Z/fit.proj$fitted.values) - ((1 - Z)/(1 - fit.proj$fitted.values)))
    w_proj1 <- w_proj1/mean(w_proj1)
    w_proj0 <- w_proj0/mean(w_proj0)
    w_proj_ate <- w_proj1 - w_proj0
    w_proj_ov <- 1 - (1 - treat) * (Z/fit.proj$fitted.values) - 
      treat * (1 - Z)/(1 - fit.proj$fitted.values)
    w_proj_ov = w_proj_ov/mean(w_proj_ov)
    #IPS-Ind
    w_ind1 <- treat * ((Z/fit.ind$fitted.values) - ((1 - Z)/(1 - fit.ind$fitted.values)))
    w_ind0 <- (1 - treat) * ((Z/fit.ind$fitted.values) - ((1 - Z)/(1 - fit.ind$fitted.values)))
    w_ind1 <- w_ind1/mean(w_ind1)
    w_ind0 <- w_ind0/mean(w_ind0)
    w_ind_ate <- w_ind1 - w_ind0
    w_ind_ov <- 1 - (1 - treat) * (Z/fit.ind$fitted.values) - 
      treat * (1 - Z)/(1 - fit.ind$fitted.values)
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
      
      # std. errors
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
mean.mc <- base::colMeans(mc, na.rm = TRUE)

true_effects  <- c(rep(late_true, 6),
                   rep(lqte_true[1], 6),
                   rep(lqte_true[2], 6),
                   rep(lqte_true[3], 6),
                   rep(lqte_true[4], 6),
                   rep(lqte_true[5], 6)
)
true_effects <- matrix(true_effects, nrow = nrow(mc), ncol = 36, byrow = TRUE)

bias_mc <- mean.mc[1:36] - colMeans(true_effects)

RMSE_mc <- sqrt(base::colMeans((mc[,1:36]-true_effects)^2, na.rm = TRUE) )



abs_bias_mc <- (base::colMeans(abs(mc[,1:36] - true_effects), na.rm = TRUE))

critical_value <- stats::qnorm(0.975)
coverage <-  colMeans(((mc[,1:36] - critical_value*mc[,37:72]) <= true_effects) * 
                        ((mc[,1:36] + critical_value*mc[,37:72]) >= true_effects),
                      na.rm = TRUE)


ks_bal <- mean.mc[73:90]
cvm_bal <- mean.mc[91:108]

#----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
if(trim == FALSE){
  #Export output: Inference
  source(here("Codes/MC_tables_inf.R"))
  write.csv(out_sim_inf, file = here("Simulations/401k/Results/LTE/", 
                                     paste0("MC_lte_inf.csv")))
  #Export output: point estimate
  source(here("Codes", "MC_tables_point.R"))
  write.csv(out_sim_point, file = here("Simulations/401k/Results/LTE/", 
                                       paste0("MC_lte_point.csv")))
  
  #Export output: balance
  source(here("Codes", "MC_tables_bal.R"))
  write.csv(out_sim_bal, file = here("Simulations/401k/Results/LTE/", 
                                     paste0("MC_lte_bal.csv")))
  #----------------------------------------------------------------------------
  # Save Dataset
  save.image(here("Simulations/401k/Results/LTE/",
                  paste0("401k_sim_lte_n", n,"_dgp_",dgp,".RData"))
  )
}
if(trim == TRUE){
  #Export output: Inference
  source(here("Codes/MC_tables_inf.R"))
  write.csv(out_sim_inf, file = here("Simulations/401k/Results/LTE/", 
                                     paste0("MC_lte_inf_trimmed.csv")))
  #Export output: point estimate
  source(here("Codes", "MC_tables_point.R"))
  write.csv(out_sim_point, file = here("Simulations/401k/Results/LTE/", 
                                       paste0("MC_lte_point_trimmed.csv")))
  
  #Export output: balance
  source(here("Codes", "MC_tables_bal.R"))
  write.csv(out_sim_bal, file = here("Simulations/401k/Results/LTE/", 
                                     paste0("MC_lte_bal_trimmed.csv")))
  #----------------------------------------------------------------------------
  # Save Dataset
  save.image(here("Simulations/401k/Results/LTE/",
                  paste0("401k_sim_lte_n", n,"_dgp_",dgp,'trimmed' ,".RData"))
  )
}


