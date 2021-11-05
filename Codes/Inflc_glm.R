inflc_glm<- function(Treat, X, pscore, weights = NULL){
  n <- length(Treat)
  X <- as.matrix(X)
  PS <- as.vector(pscore) #1/(1+exp(-X%*%beta0))
  if(is.null(weights)) weights <- rep(1, n)
  if(!is.numeric(weights)) base::stop("weights must be a NULL or a numeric vector")
  #-----------------------------------------------------------------------------
  # Avoid diving by zero
  probs.min<- 1e-6
  PS <- base::pmin(1 - probs.min, PS)
  PS <- base::pmax(probs.min, PS)
  PS <- base::as.vector(PS)
  #-----------------------------------------------------------------------------
  #Fisher Information Matrix: k x k
  fisher <- (base::crossprod(weights * X, diag(as.vector(PS*(1-PS)))) %*% (weights * X)) / n
  # Score of Logit MLE
  ps.score <- as.vector(weights * (Treat - PS)) * X     
  #Asymptotic linear representation of GLM pscore: nxk
  lin.rep.ps.mle <- ps.score %*% MASS::ginv(fisher)
  
  return(lin.rep.ps.mle)
}