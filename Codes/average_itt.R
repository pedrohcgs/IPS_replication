average_itt <- function(outcomeY){
  
  #----------------------------------------------------------------------------
  # GLM  
  glm.lin.rep <- inflc_glm(treat, baselineX, 
                           ps.glm$fitted.values)
  
  ate.ps.glm <- IPS::ATE(outcomeY, 
                         treat, 
                         baselineX, 
                         ps.glm$fitted.values, 
                         glm.lin.rep)
  
  # CBPS
  cbps.lin.rep <- inflc_cbps(treat, baselineX, 
                             ps.cbps$fitted.values, method = "exact")
  
  ate.ps.cbps <- IPS::ATE(outcomeY, 
                          treat,
                          baselineX, 
                          ps.cbps$fitted.values, 
                          cbps.lin.rep)
  
  # CBPS - overidentified
  cbps.over.lin.rep <- inflc_cbps(treat, baselineX, 
                                  ps.cbps.over$fitted.values, method = "over")
  
  ate.ps.cbps.over <- IPS::ATE(outcomeY, 
                               treat,
                               baselineX, 
                               ps.cbps.over$fitted.values, 
                               cbps.over.lin.rep)
  
  # IPS with indicator
  ate.ps.ind <- IPS::ATE(outcomeY, 
                         treat,
                         baselineX, 
                         ps.ind$fitted.values, 
                         ps.ind$lin.rep)
  
  # IPS with exponential
  ate.ps.exp <- IPS::ATE(outcomeY, 
                         treat,
                         baselineX, 
                         ps.exp$fitted.values, 
                         ps.exp$lin.rep)
  
  # IPS with projection
  ate.ps.proj <- IPS::ATE(outcomeY, 
                          treat,
                          baselineX, 
                          ps.proj$fitted.values, 
                          ps.proj$lin.rep)
  
  
  ps.ate.treat <-  rbind(cbind(ate.ps.glm$ate, 
                               ate.ps.cbps$ate,
                               ate.ps.cbps.over$ate,
                               ate.ps.ind$ate, 
                               ate.ps.exp$ate,
                               ate.ps.proj$ate
  ), 
  cbind(ate.ps.glm$ate.se, 
        ate.ps.cbps$ate.se,
        ate.ps.cbps.over$ate.se,
        ate.ps.ind$ate.se, 
        ate.ps.exp$ate.se, 
        ate.ps.proj$ate.se))
  #----------------------------------------------------------------------------
  # compute ITT on the treated using different estimators
  #----------------------------------------------------------------------------
  # GLM  
  att.ps.glm <- IPS::ATT(outcomeY, 
                         treat, 
                         baselineX, 
                         ps.glm$fitted.values, 
                         glm.lin.rep)
  
  # CBPS
  att.ps.cbps <- IPS::ATT(outcomeY, 
                          treat,
                          baselineX, 
                          ps.cbps$fitted.values, 
                          cbps.lin.rep)
  
  # CBPS - over identified
  att.ps.cbps.over <- IPS::ATT(outcomeY, 
                               treat,
                               baselineX, 
                               ps.cbps.over$fitted.values, 
                               cbps.over.lin.rep)
  
  # IPS with indicator
  att.ps.ind <- IPS::ATT(outcomeY, 
                         treat,
                         baselineX, 
                         ps.ind$fitted.values, 
                         ps.ind$lin.rep)
  
  # IPS with exponential
  att.ps.exp <- IPS::ATT(outcomeY, 
                         treat,
                         baselineX, 
                         ps.exp$fitted.values, 
                         ps.exp$lin.rep)
  
  # IPS with projection
  att.ps.proj <- IPS::ATT(outcomeY, 
                          treat,
                          baselineX, 
                          ps.proj$fitted.values, 
                          ps.proj$lin.rep)
  
  
  ps.att.treat <-  rbind(cbind(att.ps.glm$att, 
                               att.ps.cbps$att,
                               att.ps.cbps.over$att,
                               att.ps.ind$att, 
                               att.ps.exp$att,
                               att.ps.proj$att
  ), 
  cbind(att.ps.glm$att.se, 
        att.ps.cbps$att.se,
        att.ps.cbps.over$att.se,
        att.ps.ind$att.se, 
        att.ps.exp$att.se, 
        att.ps.proj$att.se))
  
  # Combine ATE and ATT
  ps.treat <- rbind(ps.ate.treat, ps.att.treat)
  colnames(ps.treat) <- c("glm","cbps-just","cbps-over","ind","exp","proj")
  rownames(ps.treat) <- c("ate-z","ate-z.se", "att-z","att-z.se")
  
  return(ps.treat)
  
}