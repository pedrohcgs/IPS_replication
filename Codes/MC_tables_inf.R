# Generate Excel files related to Inference
#----------------------------------------------------------------------------
# Summarize the simulation results
#----------------------------------------------------------------------------
# Compute empirical coverage and  asy. std deviation
emp.cov <- base::colMeans((abs(mc[,1:36] - true_effects) / 
                             mc[,37:72]) <= 1.96, na.rm = TRUE)

ASSD <- base::colMeans(sqrt(n) * mc[ ,37:72], na.rm = TRUE)
#----------------------------------------------------------------------------
#Generate output table
out_p_inf <- rep(0, 75)
dim(out_p_inf) <- c(1, 75)
out_p_inf <- data.frame(out_p_inf)

out_p_inf[1,1] <- nrep
out_p_inf[1,2] <- n
out_p_inf[1,3] <- dgp
#----------------------------------------------------------------------------
# Empirical coverage 
out_p_inf[1,4:39] <- emp.cov
# Average standard error 
out_p_inf[1,40:75] <- ASSD
#----------------------------------------------------------------------------
colnames(out_p_inf) <- c("nrep", "n", "dgp (1: correct, 2:misspec)",
                     #ATE
                     "ATE-Empcov-IPS-exp","ATE-Empcov-IPS-ind", "ATE-Empcov-IPS-proj",
                     "ATE-Empcov-CBPS-Just","ATE-Empcov-CBPS-over","ATE-Empcov-GLM",
                     #Qte-0.1
                     "QTE-0.1-Empcov-IPS-exp","QTE-0.1-Empcov-IPS-ind", "QTE-0.1-Empcov-IPS-proj",
                     "QTE-0.1-Empcov-CBPS-Just","QTE-0.1-Empcov-CBPS-over","QTE-0.1-Empcov-GLM",
                     #Qte-0.25
                     "QTE-0.25-Empcov-IPS-exp","QTE-0.25-Empcov-IPS-ind", "QTE-0.25-Empcov-IPS-proj",
                     "QTE-0.25-Empcov-CBPS-Just","QTE-0.25-Empcov-CBPS-over","QTE-0.25-Empcov-GLM",
                     #Qte-0.5 
                     "QTE-0.5-Empcov-IPS-exp","QTE-0.5-Empcov-IPS-ind", "QTE-0.5-Empcov-IPS-proj",
                     "QTE-0.5-Empcov-CBPS-Just","QTE-0.5-Empcov-CBPS-over","QTE-0.5-Empcov-GLM",
                     #Qte-0.75 
                     "QTE-0.75-Empcov-IPS-exp","QTE-0.75-Empcov-IPS-ind", "QTE-0.75-Empcov-IPS-proj",
                     "QTE-0.75-Empcov-CBPS-Just","QTE-0.75-Empcov-CBPS-over","QTE-0.75-Empcov-GLM",
                     #Qte-0.9 
                     "QTE-0.9-Empcov-IPS-exp","QTE-0.9-Empcov-IPS-ind", "QTE-0.9-Empcov-IPS-proj",
                     "QTE-0.9-Empcov-CBPS-Just","QTE-0.9-Empcov-CBPS-over","QTE-0.9-Empcov-GLM",
                     
                     #ATE
                     "ATE-ASSD-IPS-exp","ATE-ASSD-IPS-ind", "ATE-ASSD-IPS-proj",
                     "ATE-ASSD-CBPS-Just","ATE-ASSD-CBPS-over","ATE-ASSD-GLM",
                     #Qte-0.1
                     "QTE-0.1-ASSD-IPS-exp","QTE-0.1-ASSD-IPS-ind", "QTE-0.1-ASSD-IPS-proj",
                     "QTE-0.1-ASSD-CBPS-Just","QTE-0.1-ASSD-CBPS-over","QTE-0.1-ASSD-GLM",
                     #Qte-0.25
                     "QTE-0.25-ASSD-IPS-exp","QTE-0.25-ASSD-IPS-ind", "QTE-0.25-ASSD-IPS-proj",
                     "QTE-0.25-ASSD-CBPS-Just","QTE-0.25-ASSD-CBPS-over","QTE-0.25-ASSD-GLM",
                     #Qte-0.5 
                     "QTE-0.5-ASSD-IPS-exp","QTE-0.5-ASSD-IPS-ind", "QTE-0.5-ASSD-IPS-proj",
                     "QTE-0.5-ASSD-CBPS-Just","QTE-0.5-ASSD-CBPS-over","QTE-0.5-ASSD-GLM",
                     #Qte-0.75 
                     "QTE-0.75-ASSD-IPS-exp","QTE-0.75-ASSD-IPS-ind", "QTE-0.75-ASSD-IPS-proj",
                     "QTE-0.75-ASSD-CBPS-Just","QTE-0.75-ASSD-CBPS-over","QTE-0.75-ASSD-GLM",
                     #Qte-0.9
                     "QTE-0.9-ASSD-IPS-exp","QTE-0.9-ASSD-IPS-ind", "QTE-0.9-ASSD-IPS-proj",
                     "QTE-0.9-ASSD-CBPS-Just","QTE-0.9-ASSD-CBPS-over","QTE-0.9-ASSD-GLM"
      
)
####################################################################
if(exists("out_sim_inf") == FALSE){out_sim_inf <- out_p_inf}
if(identical(out_sim_inf,out_p_inf)==FALSE) {out_sim_inf <- rbind(out_sim_inf, out_p_inf)}

