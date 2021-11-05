dgps_ips_late <- function(n, dgp, scale.ps=1, rho=1){
  #----------------------------------------------------------------------------
  
  #----------------------------------------------------------------------------
  #----------------------------------------------------------------------------
  # Kand and Schafer design - Correct model
  #----------------------------------------------------------------------------
  #----------------------------------------------------------------------------
  if (dgp==1){
    # Generate X
    X1 <- stats::rnorm(n)
    X2 <- stats::rnorm(n)
    X3 <- stats::rnorm(n)
    X4 <- stats::rnorm(n)
    X = cbind(X1, X2, X3, X4)
    
    # Beta of pscores
    b.ps <- c(-1, 0.5, -0.25, -0.1)
    # Generate pscore and Treatment Status
    b.ps.Times.x <- X %*% b.ps
    pr.z <- 1/(1+exp(- scale.ps * b.ps.Times.x))
    Z <- stats::rbinom(n, size=1, prob = pr.z)
    
    # the potential outcomes
    # Beta of regression
    b.reg <- c(27.4, 13.7, 13.7, 13.7)
    b.reg.Times.x <- X %*% b.reg
    
    u0 <- stats::rnorm(n)
    u1 <- stats::rnorm(n)
    
    Y0 <- 200 - b.reg.Times.x + u0
    Y1 <- 210 + b.reg.Times.x + u1
    
    # The potential treatment status
    d0 <- rep(0,n)
    d1.index <- 2 + rho * (Y1 - Y0)/20
    pr.d1 <- 1/(1+exp(- d1.index))
    d1 <- stats::rbinom(n, size=1, prob = pr.d1)
    
    # Observed D
    Treat <- d1*Z + d0 * (1-Z)
    
    # Observed outcome
    Y <- Y1*Treat + Y0*(1-Treat)
    
    # Dataset
    datamatrix <- data.frame(cbind(Y, Treat, Z, X, X))
    
    colnames(datamatrix) <- c("Y", "Treat", "Instrument",
                              paste("X", 1:dim(X)[2], sep = ""),
                              paste("trueX", 1:dim(X)[2], sep = ""))
  }
  
  #----------------------------------------------------------------------------
  # Kand and Schafer design - Misspecified model
  #----------------------------------------------------------------------------
  if (dgp==2){
    # Generate X
    X1 <- stats::rnorm(n)
    X2 <- stats::rnorm(n)
    X3 <- stats::rnorm(n)
    X4 <- stats::rnorm(n)
    X = cbind(X1,X2,X3,X4)
    
    # Beta of pscores
    b.ps <- c(-1, 0.5, -0.25, -0.1)
    # Generate pscore and Treatment Status
    b.ps.Times.x <- X %*% b.ps
    pr.z <- 1/(1+exp(- scale.ps * b.ps.Times.x))
    Z <- stats::rbinom(n, size=1, prob = pr.z)
    
    
    # the potential outcomes
    # Beta of regression
    b.reg <- c(27.4, 13.7, 13.7, 13.7)
    b.reg.Times.x <- X %*% b.reg
    
    u0 <- stats::rnorm(n)
    u1 <- stats::rnorm(n)
    
    Y0 <- 200 - b.reg.Times.x + u0
    Y1 <- 210 + b.reg.Times.x + u1
    
    # The potential treatment status
    d0 <- rep(0,n)
    d1.index <- 2 + rho * (Y1 - Y0)/20
    pr.d1 <- 1/(1+exp(- d1.index))
    d1 <- stats::rbinom(n, size=1, prob = pr.d1)
    
    
    # Observed D
    Treat <- d1*Z + d0 * (1-Z)
    
    # Observed outcome
    Y <- Y1*Treat + Y0*(1-Treat)
    
    # Dataset
    W1 <- exp(X1/2) 
    W2 <- X2/(1+exp(X1))
    W3 <- (X1*X3/25 + 0.6)^3 
    W4 <- (X1 + X4 + 20)^2 
    W <- cbind(W1, W2, W3, W4)
    datamatrix <- data.frame(cbind(Y, Treat, Z, W, X))
    
    colnames(datamatrix) <- c("Y", "Treat", "Instrument",
                              paste("X", 1:dim(X)[2], sep = ""),
                              paste("trueX", 1:dim(X)[2], sep = ""))
  }
  #----------------------------------------------------------------------------
  
  
  return(datamatrix)
}
