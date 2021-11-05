# Generate Bias and RMSE
#----------------------------------------------------------------------------
# Summarize the simulation results 
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Compute Median bias and RMSE
bias_hat = base::colMeans(mc[,1:36], na.rm = T) - true_effects[1,]
RMSE_hat = sqrt(base::colMeans((mc[,1:36] - true_effects)^2,
                                na.rm = T))

#----------------------------------------------------------------------------
#Generate output table
out_p_point <- rep(0, 75)
dim(out_p_point) <- c(1, 75)
out_p_point <- data.frame(out_p_point)

out_p_point[1,1] <- nrep
out_p_point[1,2] <- n
out_p_point[1,3] <- dgp
#----------------------------------------------------------------------------
# Bias
out_p_point[1,4:39] <- bias_hat
# RMSE
out_p_point[1,40:75] <- RMSE_hat
#----------------------------------------------------------------------------
colnames(out_p_point) <- c("nrep", "n", "dgp (1: correct, 2:misspec)",
                           #ATE
                           "ATE-bias-IPS-exp","ATE-bias-IPS-ind", "ATE-bias-IPS-proj",
                           "ATE-bias-CBPS-just","ATE-bias-CBPS-over","ATE-bias-GLM",
                           #Qte-0.1                     
                           "QTE-0.10-bias-IPS-exp","QTE-0.10-bias-IPS-ind", "QTE-0.10-bias-IPS-proj",
                           "QTE-0.10-bias-CBPS-Just","QTE-0.10-bias-CBPS-over", "QTE-0.10-bias-GLM",
                           #Qte-0.25
                           "QTE-0.25-bias-IPS-exp","QTE-0.25-bias-IPS-ind", "QTE-0.25-bias-IPS-proj",
                           "QTE-0.25-bias-CBPS-Just","QTE-0.25-bias-CBPS-over","QTE-0.25-bias-GLM",
                           #Qte-0.5 
                           "QTE-0.5-bias-IPS-exp","QTE-0.5-bias-IPS-ind", "QTE-0.5-bias-IPS-proj",
                           "QTE-0.5-bias-CBPS-Just","QTE-0.5-bias-CBPS-over","QTE-0.5-bias-GLM",
                           #Qte-0.75 
                           "QTE-0.75-bias-IPS-exp","QTE-0.75-bias-IPS-ind", "QTE-0.75-bias-IPS-proj",
                           "QTE-0.75-bias-CBPS-Just","QTE-0.75-bias-CBPS-over","QTE-0.75-bias-GLM",
                           #Qte-0.9
                           "QTE-0.9-bias-IPS-exp","QTE-0.9-bias-IPS-ind", "QTE-0.9-bias-IPS-proj",
                           "QTE-0.9-bias-CBPS-Just","QTE-0.9-bias-CBPS-over","QTE-0.9-bias-GLM",
                           
                           #ATE
                           "ATE-RMSE-IPS-exp","ATE-RMSE-IPS-ind", "ATE-RMSE-IPS-proj",
                           "ATE-RMSE-CBPS-Just","ATE-RMSE-CBPS-over","ATE-RMSE-GLM",
                           #Qte-0.1                     
                           "QTE-0.10-RMSE-IPS-exp","QTE-0.10-RMSE-IPS-ind", "QTE-0.10-RMSE-IPS-proj",
                           "QTE-0.10-RMSE-CBPS-Just","QTE-0.10-RMSE-CBPS-over","QTE-0.10-RMSE-GLM",
                           #Qte-0.25
                           "QTE-0.25-RMSE-IPS-exp","QTE-0.25-RMSE-IPS-ind", "QTE-0.25-RMSE-IPS-proj",
                           "QTE-0.25-RMSE-CBPS-Just","QTE-0.25-RMSE-CBPS-over","QTE-0.25-RMSE-GLM",
                           #Qte-0.5 
                           "QTE-0.5-RMSE-IPS-exp","QTE-0.5-RMSE-IPS-ind", "QTE-0.5-RMSE-IPS-proj",
                           "QTE-0.5-RMSE-CBPS-Just","QTE-0.5-RMSE-CBPS-over","QTE-0.5-RMSE-GLM",
                           #Qte-0.75 
                           "QTE-0.75-RMSE-IPS-exp","QTE-0.75-RMSE-IPS-ind", "QTE-0.75-RMSE-IPS-proj",
                           "QTE-0.75-RMSE-CBPS-Just","QTE-0.75-RMSE-CBPS-over","QTE-0.75-RMSE-GLM",
                           #Qte-0.9
                           "QTE-0.9-RMSE-IPS-exp","QTE-0.9-RMSE-IPS-ind", "QTE-0.9-RMSE-IPS-proj",
                           "QTE-0.9-RMSE-CBPS-Just","QTE-0.9-RMSE-CBPS-over","QTE-0.9-RMSE-GLM"
                           
)
####################################################################
if(exists("out_sim_point")==FALSE){out_sim_point <- out_p_point}
if(identical(out_sim_point, out_p_point)==FALSE) {out_sim_point <- rbind(out_sim_point,out_p_point)}
