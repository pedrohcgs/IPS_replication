quantile_itt <- function(outcomeY, tauY){
  
  
  # Get some lin rep of other methods
  glm.lin.rep <- inflc_glm(treat, baselineX, 
                           ps.glm$fitted.values)
  
  cbps.lin.rep <- inflc_cbps(treat, baselineX, 
                             ps.cbps$fitted.values, method = "exact")
  # CBPS - overidentified
  cbps.over.lin.rep <- inflc_cbps(treat, baselineX, 
                                  ps.cbps.over$fitted.values, method = "over")
  # Compute the QTE & QTT estimators
  # indicator
  qte.ps.ind <- IPS::QTE(outcomeY, 
                         treat, 
                         baselineX, 
                         ps.ind$fitted.values, 
                         ps.ind$lin.rep,
                         tauY,
                         bw = "nrd0",
                         whs = NULL )
  
  # exponential
  qte.ps.exp <- IPS::QTE(outcomeY, 
                         treat, 
                         baselineX, 
                         ps.exp$fitted.values, 
                         ps.exp$lin.rep,
                         tauY, bw = "nrd0",
                         whs = NULL )
  
  # projection 
  qte.ps.proj <- IPS::QTE(outcomeY,
                          treat,
                          baselineX,
                          ps.proj$fitted.values,
                          ps.proj$lin.rep,
                          tauY, bw = "nrd0",
                          whs = NULL )
  
  # GLM
  qte.ps.glm <- IPS::QTE(outcomeY, 
                         treat, 
                         baselineX,
                         ps.glm$fitted.values, 
                         glm.lin.rep,
                         tauY, bw = "nrd0",
                         whs = NULL )
  
  # CBPS
  qte.ps.cbps <- IPS::QTE(outcomeY, 
                          treat, 
                          baselineX,
                          ps.cbps$fitted.values, 
                          cbps.lin.rep,
                          tauY, bw = "nrd0",
                          whs = NULL )
  # CBPS - overidentified
  qte.ps.cbps.over <- IPS::QTE(outcomeY, 
                               treat, 
                               baselineX,
                               ps.cbps.over$fitted.values, 
                               cbps.over.lin.rep,
                               tauY, bw = "nrd0",
                               whs = NULL)
  
  ret.qte <-  rbind(
    cbind(qte.ps.glm$qte[1], 
          qte.ps.cbps$qte[1], 
          qte.ps.cbps.over$qte[1], 
          qte.ps.ind$qte[1],
          qte.ps.exp$qte[1],
          qte.ps.proj$qte[1]
    ),
    cbind(qte.ps.glm$qte.se[1], 
          qte.ps.cbps$qte.se[1], 
          qte.ps.cbps.over$qte.se[1], 
          qte.ps.ind$qte.se[1], 
          qte.ps.exp$qte.se[1],
          qte.ps.proj$qte.se[1]
    ),
    cbind(qte.ps.glm$qte[2], 
          qte.ps.cbps$qte[2], 
          qte.ps.cbps.over$qte[2], 
          qte.ps.ind$qte[2],
          qte.ps.exp$qte[2],
          qte.ps.proj$qte[2]
    ),
    cbind(qte.ps.glm$qte.se[2], 
          qte.ps.cbps$qte.se[2], 
          qte.ps.cbps.over$qte.se[2], 
          qte.ps.ind$qte.se[2], 
          qte.ps.exp$qte.se[2],
          qte.ps.proj$qte.se[2]
    ),
    cbind(qte.ps.glm$qte[3], 
          qte.ps.cbps$qte[3], 
          qte.ps.cbps.over$qte[3], 
          qte.ps.ind$qte[3],
          qte.ps.exp$qte[3],
          qte.ps.proj$qte[3]
    ),
    cbind(qte.ps.glm$qte.se[3], 
          qte.ps.cbps$qte.se[3], 
          qte.ps.cbps.over$qte.se[3], 
          qte.ps.ind$qte.se[3], 
          qte.ps.exp$qte.se[3],
          qte.ps.proj$qte.se[3]
    ),
    cbind(qte.ps.glm$qte[4], 
          qte.ps.cbps$qte[4], 
          qte.ps.cbps.over$qte[4], 
          qte.ps.ind$qte[4],
          qte.ps.exp$qte[4],
          qte.ps.proj$qte[4]
    ),
    cbind(qte.ps.glm$qte.se[4], 
          qte.ps.cbps$qte.se[4], 
          qte.ps.cbps.over$qte.se[4], 
          qte.ps.ind$qte.se[4], 
          qte.ps.exp$qte.se[4],
          qte.ps.proj$qte.se[4]
    ),
    cbind(qte.ps.glm$qte[5], 
          qte.ps.cbps$qte[5], 
          qte.ps.cbps.over$qte[5], 
          qte.ps.ind$qte[5],
          qte.ps.exp$qte[5],
          qte.ps.proj$qte[5]
    ),
    cbind(qte.ps.glm$qte.se[5], 
          qte.ps.cbps$qte.se[5], 
          qte.ps.cbps.over$qte.se[5], 
          qte.ps.ind$qte.se[5], 
          qte.ps.exp$qte.se[5],
          qte.ps.proj$qte.se[5]
    )
  )
  
  #----------------------------------------------------
  # QTT estimator
  # indicator
  qtt.ps.ind <- IPS::QTT(outcomeY, 
                         treat, 
                         baselineX, 
                         ps.ind$fitted.values, 
                         ps.ind$lin.rep,
                         tauY, bw = "nrd0",
                         whs = NULL )
  
  # exponential
  qtt.ps.exp <- IPS::QTT(outcomeY, 
                         treat, 
                         baselineX, 
                         ps.exp$fitted.values, 
                         ps.exp$lin.rep,
                         tauY, bw = "nrd0",
                         whs = NULL )
  
  # projection 
  qtt.ps.proj <- IPS::QTT(outcomeY,
                          treat,
                          baselineX,
                          ps.proj$fitted.values,
                          ps.proj$lin.rep,
                          tauY, bw = "nrd0",
                          whs = NULL )
  
  # GLM
  qtt.ps.glm <- IPS::QTT(outcomeY, 
                         treat, 
                         baselineX,
                         ps.glm$fitted.values, 
                         glm.lin.rep,
                         tauY, bw = "nrd0",
                         whs = NULL )
  
  # CBPS
  qtt.ps.cbps <- IPS::QTT(outcomeY, 
                          treat, 
                          baselineX,
                          ps.cbps$fitted.values, 
                          cbps.lin.rep,
                          tauY, bw = "nrd0",
                          whs = NULL )
  
  
  # CBPS - overidentified
  qtt.ps.cbps.over <- IPS::QTT(outcomeY, 
                               treat, 
                               baselineX,
                               ps.cbps.over$fitted.values, 
                               cbps.over.lin.rep,
                               tauY, bw = "nrd0",
                               whs = NULL)
  
  
  ret.qtt <-  rbind(
    cbind(qtt.ps.glm$qtt[1], 
          qtt.ps.cbps$qtt[1], 
          qtt.ps.cbps.over$qtt[1], 
          qtt.ps.ind$qtt[1],
          qtt.ps.exp$qtt[1],
          qtt.ps.proj$qtt[1]
    ),
    cbind(qtt.ps.glm$qtt.se[1], 
          qtt.ps.cbps$qtt.se[1], 
          qtt.ps.cbps.over$qtt.se[1], 
          qtt.ps.ind$qtt.se[1], 
          qtt.ps.exp$qtt.se[1],
          qtt.ps.proj$qtt.se[1]
    ),
    cbind(qtt.ps.glm$qtt[2], 
          qtt.ps.cbps$qtt[2], 
          qtt.ps.cbps.over$qtt[2], 
          qtt.ps.ind$qtt[2],
          qtt.ps.exp$qtt[2],
          qtt.ps.proj$qtt[2]
    ),
    cbind(qtt.ps.glm$qtt.se[2], 
          qtt.ps.cbps$qtt.se[2], 
          qtt.ps.cbps.over$qtt.se[2], 
          qtt.ps.ind$qtt.se[2], 
          qtt.ps.exp$qtt.se[2],
          qtt.ps.proj$qtt.se[2]
    ),
    cbind(qtt.ps.glm$qtt[3], 
          qtt.ps.cbps$qtt[3], 
          qtt.ps.cbps.over$qtt[3], 
          qtt.ps.ind$qtt[3],
          qtt.ps.exp$qtt[3],
          qtt.ps.proj$qtt[3]
    ),
    cbind(qtt.ps.glm$qtt.se[3], 
          qtt.ps.cbps$qtt.se[3], 
          qtt.ps.cbps.over$qtt.se[3], 
          qtt.ps.ind$qtt.se[3], 
          qtt.ps.exp$qtt.se[3],
          qtt.ps.proj$qtt.se[3]
    ),
    cbind(qtt.ps.glm$qtt[4], 
          qtt.ps.cbps$qtt[4], 
          qtt.ps.cbps.over$qtt[4], 
          qtt.ps.ind$qtt[4],
          qtt.ps.exp$qtt[4],
          qtt.ps.proj$qtt[4]
    ),
    cbind(qtt.ps.glm$qtt.se[4], 
          qtt.ps.cbps$qtt.se[4], 
          qtt.ps.cbps.over$qtt.se[4], 
          qtt.ps.ind$qtt.se[4], 
          qtt.ps.exp$qtt.se[4],
          qtt.ps.proj$qtt.se[4]
    ),
    cbind(qtt.ps.glm$qtt[5], 
          qtt.ps.cbps$qtt[5], 
          qtt.ps.cbps.over$qtt[5], 
          qtt.ps.ind$qtt[5],
          qtt.ps.exp$qtt[5],
          qtt.ps.proj$qtt[5]
    ),
    cbind(qtt.ps.glm$qtt.se[5], 
          qtt.ps.cbps$qtt.se[5], 
          qtt.ps.cbps.over$qtt.se[5], 
          qtt.ps.ind$qtt.se[5], 
          qtt.ps.exp$qtt.se[5],
          qtt.ps.proj$qtt.se[5]
    )
  )
  
  
  
  ps.q.treat <- rbind(ret.qte, ret.qtt)
  colnames(ps.q.treat) <- c("glm","cbps-just","cbps-over","ind","exp","proj")
  rownames(ps.q.treat) <- c("qte-z(0.10)","qte-z(0.10).se",
                            "qte-z(0.25)","qte-z(0.25).se",
                            "qte-z(0.5)","qte-z(0.5).se",
                            "qte-z(0.75)","qte-z(0.75).se",
                            "qte-z(0.90)","qte-z(0.90).se",
                            
                            "qtt-z(0.10)","qtt-z(0.10).se",
                            "qtt-z(0.25)","qtt-z(0.25).se",
                            "qtt-z(0.5)","qtt-z(0.5).se",
                            "qtt-z(0.75)","qtt-z(0.75).se",
                            "qtt-z(0.90)","qtt-z(0.90).se"
  )
  
  return(ps.q.treat)
  
}
