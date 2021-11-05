# Generate Excel files related to Balance
#----------------------------------------------------------------------------
# Summarize the simulation results
#----------------------------------------------------------------------------
# Compute average of KS and CvM balance tests
#----------------------------------------------------------------------------
#Generate output table
out_dim <- length(ks_bal)*2 + 3
out_p_bal <- rep(0, out_dim)
dim(out_p_bal) <- c(1, out_dim)
out_p_bal <- data.frame(out_p_bal)

out_p_bal[1,1] <- nrep
out_p_bal[1,2] <- n
out_p_bal[1,3] <- dgp
#----------------------------------------------------------------------------
# Empirical coverage 
out_p_bal[1,4:(length(ks_bal)+3)] <- ks_bal
# Average standard error 
out_p_bal[1,(length(ks_bal)+4):(out_dim)] <- cvm_bal
#----------------------------------------------------------------------------
colnames(out_p_bal) <- c("nrep", "n", "dgp (1: correct, 2:misspec)",
                         #KS
                         "ks-IPS-exp","ks-IPS-ind", "ks-IPS-proj",
                         "ks-CBPS-Just","ks-CBPS-over","ks-GLM",
                         
                         "ks-IPS-exp_1","ks-IPS-ind_1", "ks-IPS-proj_1",
                         "ks-CBPS-Just_1","ks-CBPS-over_1","ks-GLM_1",
                         
                         "ks-IPS-exp_0","ks-IPS-ind_0", "ks-IPS-proj_0",
                         "ks-CBPS-Just_0","ks-CBPS-over_0","ks-GLM_0",
                         
                         #CvM
                         "cvm-IPS-exp","cvm-IPS-ind", "cvm-IPS-proj",
                         "cvm-CBPS-Just","cvm-CBPS-over","cvm-GLM",
                         
                         "cvm-IPS-exp_1","cvm-IPS-ind_1", "cvm-IPS-proj_1",
                         "cvm-CBPS-Just_1","cvm-CBPS-over_1","cvm-GLM_1",
                         
                         "cvm-IPS-exp_0","cvm-IPS-ind_0", "cvm-IPS-proj_0",
                         "cvm-CBPS-Just_0","cvm-CBPS-over_0","cvm-GLM_0"
                         
)
####################################################################
if(exists("out_sim_bal") == FALSE){out_sim_bal <- out_p_bal}
if(identical(out_sim_bal,out_p_bal)==FALSE) {out_sim_bal <- rbind(out_sim_bal, out_p_bal)}

