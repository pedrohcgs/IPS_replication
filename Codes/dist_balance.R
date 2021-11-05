dist_balance <- function(covs, weights, chunk = 1000){
  #-----------------------------------------------------------------------------
  # Ensure covs are stored as matrix
  covs <- as.matrix(covs)
  # Ensure weights are stored as numetic
  weights <- as.numeric(weights)
  #-----------------------------------------------------------------------------
  # Define variables to be used in the loop
  # Number of covariates
  k_dim = dim(covs)[2]
  # number of observations
  n_obs = dim(covs)[1]
  
  if(n_obs != length(weights)){
    stop("Weights must have same number of obs. as covs")
  }
  
  # Initialize `Rw` row vector (n_obs dimension)
  cdf_balance <- rep(0,n_obs)
  
  # We split n columns into l tiles, each with chunk columns
  l <- floor(n_obs/chunk) + 1
  
  #--------------------------------
  for (i in 1:l) {
    start <- min(chunk * (i - 1) + 1, n_obs)
    end <- min(chunk * i, n_obs)
    indicators <- base::outer(covs[,1], covs[start:end,1], "<=")
    if (k_dim>1) {
      for ( bb in 2:k_dim){
        indicators <- indicators * 
          base::outer(covs[, bb], covs[start:end, bb], "<=")
      }
    }
    
    cdf_balance[start:end] <- base::colMeans(weights * indicators)
  }
  
  ks <- sqrt(n_obs) * max(cdf_balance)
  cvm <- sum(cdf_balance^2)
  
  ret <- list(ks_bal = ks,
              cvm_bal = cvm)
}
