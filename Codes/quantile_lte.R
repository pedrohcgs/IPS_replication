quantile_lte <- function(outcomeY, tauY){
  
  
  # Get some lin rep of other methods
  glm.lin.rep <- inflc_glm(z, baselineX, 
                           ps.glm$fitted.values)
  
  cbps.lin.rep <- inflc_cbps(z, baselineX, 
                             ps.cbps$fitted.values, method = "exact")
  # CBPS - overidentified
  cbps.over.lin.rep <- inflc_cbps(z, baselineX, 
                                  ps.cbps.over$fitted.values, method = "over")
  # Compute the QTE & QTT estimators
  # indicator
  lqte.ps.ind <- IPS::LQTE(outcomeY, 
                           z,treat, 
                           baselineX, 
                           ps.ind$fitted.values, 
                           ps.ind$lin.rep,
                           tauY, bw = "nrd0")
  
  # exponential
  lqte.ps.exp <- IPS::LQTE(outcomeY, 
                           z,treat, 
                           baselineX, 
                           ps.exp$fitted.values, 
                           ps.exp$lin.rep,
                           tauY, bw = "nrd0")
  
  # projection 
  lqte.ps.proj <- IPS::LQTE(outcomeY,
                            z,treat,
                            baselineX,
                            ps.proj$fitted.values,
                            ps.proj$lin.rep,
                            tauY, bw = "nrd0")
  
  # GLM
  lqte.ps.glm <- IPS::LQTE(outcomeY, 
                           z,treat, 
                           baselineX,
                           ps.glm$fitted.values, 
                           glm.lin.rep,
                           tauY, bw = "nrd0")
  
  # CBPS
  lqte.ps.cbps <- IPS::LQTE(outcomeY, 
                            z,treat, 
                            baselineX,
                            ps.cbps$fitted.values, 
                            cbps.lin.rep,
                            tauY, bw = "nrd0" )
  # CBPS - overidentified
  lqte.ps.cbps.over <- IPS::LQTE(outcomeY, 
                                 z,treat, 
                                 baselineX,
                                 ps.cbps.over$fitted.values, 
                                 cbps.over.lin.rep,
                                 tauY, bw = "nrd0")
  
  ret.qte <-  rbind(
    cbind(lqte.ps.glm$lqte[1], 
          lqte.ps.cbps$lqte[1], 
          lqte.ps.cbps.over$lqte[1], 
          lqte.ps.ind$lqte[1],
          lqte.ps.exp$lqte[1],
          lqte.ps.proj$lqte[1]
    ),
    cbind(lqte.ps.glm$lqte.se[1], 
          lqte.ps.cbps$lqte.se[1], 
          lqte.ps.cbps.over$lqte.se[1], 
          lqte.ps.ind$lqte.se[1], 
          lqte.ps.exp$lqte.se[1],
          lqte.ps.proj$lqte.se[1]
    ),
    cbind(lqte.ps.glm$lqte[2], 
          lqte.ps.cbps$lqte[2], 
          lqte.ps.cbps.over$lqte[2], 
          lqte.ps.ind$lqte[2],
          lqte.ps.exp$lqte[2],
          lqte.ps.proj$lqte[2]
    ),
    cbind(lqte.ps.glm$lqte.se[2], 
          lqte.ps.cbps$lqte.se[2], 
          lqte.ps.cbps.over$lqte.se[2], 
          lqte.ps.ind$lqte.se[2], 
          lqte.ps.exp$lqte.se[2],
          lqte.ps.proj$lqte.se[2]
    ),
    cbind(lqte.ps.glm$lqte[3], 
          lqte.ps.cbps$lqte[3], 
          lqte.ps.cbps.over$lqte[3], 
          lqte.ps.ind$lqte[3],
          lqte.ps.exp$lqte[3],
          lqte.ps.proj$lqte[3]
    ),
    cbind(lqte.ps.glm$lqte.se[3], 
          lqte.ps.cbps$lqte.se[3], 
          lqte.ps.cbps.over$lqte.se[3], 
          lqte.ps.ind$lqte.se[3], 
          lqte.ps.exp$lqte.se[3],
          lqte.ps.proj$lqte.se[3]
    ),
    cbind(lqte.ps.glm$lqte[4], 
          lqte.ps.cbps$lqte[4], 
          lqte.ps.cbps.over$lqte[4], 
          lqte.ps.ind$lqte[4],
          lqte.ps.exp$lqte[4],
          lqte.ps.proj$lqte[4]
    ),
    cbind(lqte.ps.glm$lqte.se[4], 
          lqte.ps.cbps$lqte.se[4], 
          lqte.ps.cbps.over$lqte.se[4], 
          lqte.ps.ind$lqte.se[4], 
          lqte.ps.exp$lqte.se[4],
          lqte.ps.proj$lqte.se[4]
    ),
    cbind(lqte.ps.glm$lqte[5], 
          lqte.ps.cbps$lqte[5], 
          lqte.ps.cbps.over$lqte[5], 
          lqte.ps.ind$lqte[5],
          lqte.ps.exp$lqte[5],
          lqte.ps.proj$lqte[5]
    ),
    cbind(lqte.ps.glm$lqte.se[5], 
          lqte.ps.cbps$lqte.se[5], 
          lqte.ps.cbps.over$lqte.se[5], 
          lqte.ps.ind$lqte.se[5], 
          lqte.ps.exp$lqte.se[5],
          lqte.ps.proj$lqte.se[5]
    )
    
  )
  
  
  
  
  ps.q.treat <- rbind(ret.qte)
  colnames(ps.q.treat) <- c("glm","cbps-just","cbps-over","ind","exp","proj")
  rownames(ps.q.treat) <- c("lqte-z(0.10)","lqte-z(0.10).se",
                            "lqte-z(0.25)","lqte-z(0.25).se",
                            "lqte-z(0.5)","lqte-z(0.5).se",
                            "lqte-z(0.75)","lqte-z(0.75).se",
                            "lqte-z(0.90)","lqte-z(0.90).se")
  
  return(ps.q.treat)
  
}
