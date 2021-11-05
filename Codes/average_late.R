average_late <- function(outcomeY){
  
  #----------------------------------------------------------------------------
  # GLM  
  glm.lin.rep <- inflc_glm(z, baselineX, 
                           ps.glm$fitted.values)
  
  late.ps.glm <- IPS::LATE(outcomeY, 
                          z,
                          treat, 
                         baselineX, 
                         ps.glm$fitted.values, 
                         glm.lin.rep)
  
  # CBPS
  cbps.lin.rep <- inflc_cbps(z, baselineX, 
                             ps.cbps$fitted.values, method = "exact")
  
  late.ps.cbps <- IPS::LATE(outcomeY, 
                          z, treat,
                          baselineX, 
                          ps.cbps$fitted.values, 
                          cbps.lin.rep)
  
  # CBPS - overidentified
  cbps.over.lin.rep <- inflc_cbps(z, baselineX, 
                                  ps.cbps.over$fitted.values, method = "over")
  
  late.ps.cbps.over <- IPS::LATE(outcomeY, 
                               z,treat,
                               baselineX, 
                               ps.cbps.over$fitted.values, 
                               cbps.over.lin.rep)
  
  # IPS with indicator
  late.ps.ind <- IPS::LATE(outcomeY, 
                         z,treat,
                         baselineX, 
                         ps.ind$fitted.values, 
                         ps.ind$lin.rep)
  
  # IPS with exponential
  late.ps.exp <- IPS::LATE(outcomeY, 
                          z,treat,
                         baselineX, 
                         ps.exp$fitted.values, 
                         ps.exp$lin.rep)
  
  # IPS with projection
  late.ps.proj <- IPS::LATE(outcomeY, 
                           z,treat,
                          baselineX, 
                          ps.proj$fitted.values, 
                          ps.proj$lin.rep)
  
  
  ps.late.treat <-  rbind(cbind(late.ps.glm$late, 
                               late.ps.cbps$late,
                               late.ps.cbps.over$late,
                               late.ps.ind$late, 
                               late.ps.exp$late,
                               late.ps.proj$late
  ), 
  cbind(late.ps.glm$late.se, 
        late.ps.cbps$late.se,
        late.ps.cbps.over$late.se,
        late.ps.ind$late.se, 
        late.ps.exp$late.se, 
        late.ps.proj$late.se))
  
  # Combine ATE and ATT
  ps.treat <- rbind(ps.late.treat)
  colnames(ps.treat) <- c("glm","cbps-just","cbps-over","ind","exp","proj")
  rownames(ps.treat) <- c("late-z","late-z.se")
  
  return(ps.treat)
  
}