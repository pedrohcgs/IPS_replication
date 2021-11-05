inflc_cbps<- function(Treat, X, pscore, weights = NULL, method = "exact"){
  Treat <- as.vector(Treat)
  n <- length(Treat)
  X <- as.matrix(X)
  k <- ncol(X)
  PS <- as.vector(pscore)
  if(is.null(weights)) weights <- rep(1, n)
  if(!is.numeric(weights)) base::stop("weights must be a NULL or a numeric vector")
  #-----------------------------------------------------------------------------
  # Avoid diving by zero
  probs.min<- 1e-6
  PS <- base::pmin(1-probs.min, PS)
  PS <- base::pmax(probs.min, PS)
  PS <- base::as.vector(PS)
  #-----------------------------------------------------------------------------
  if (method == "exact"){
    
    aux <- as.vector(weights * (Treat - PS)/(PS*(1-PS))) 
    g <- aux * X  # n x k 
    aux.dot <- as.vector(weights * (((Treat - PS)^2))/(PS*(1-PS)))
    G <- base::crossprod( aux.dot * X, X) / n 
    infl  <- g %*% MASS::ginv(G)   # n x k  k x k
    
  } else if (method == "over"){
    
    aux1 <- as.vector(weights * (Treat - PS)) 
    aux2 <- as.vector(weights * (Treat - PS)/(PS*(1-PS))) 
    g <- cbind( aux1 * X, aux2 * X) # n x 2k
    aux1.dot <- as.vector(weights * PS*(1 - PS)) 
    aux2.dot <- as.vector( - weights * (((Treat - PS)^2))/(PS*(1-PS)))
    G <- t(cbind(base::crossprod( aux1.dot * X, X), base::crossprod( aux2.dot * X, X))) / n  # 2k x k
    
    Omega <-  rbind(cbind(base::crossprod(weights *PS*(1 - PS)* X, X), 
                          base::crossprod(weights *X, X)),
                    cbind(base::crossprod(weights *X, X), 
                          base::crossprod(weights *X/(PS*(1 - PS)), X))) /n   # 2k x 2k
    
    infl <- -( g %*% MASS::ginv(Omega) %*% G) %*% MASS::ginv( t(G) %*% MASS::ginv(Omega) %*% G) 
    
  }
  
  return(infl)
}