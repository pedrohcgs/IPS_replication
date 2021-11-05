#late = -3.461332
# late = 25.56

late <- function(n, scale.ps=1, rho = 1){
  
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
  z <- stats::rbinom(n, size=1, prob = pr.z)
   
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
  #d1.index <- 1 + rho * (Y1 - Y0 - 10)/30
  pr.d1 <- (1/(1+exp(- d1.index)))
  d1 <- stats::rbinom(n, size=1, prob = pr.d1)
  
  #d <- d1*z + d0 * (1-z)
  #y <- Y1*d + Y0*(1-d)
  
  ## First subindex is for d, second for z
  #w11.ps <-  d * (z/pr.z)
  #w10.ps <-  d * ((1 - z)/(1 - pr.z))
  #w01.ps <-  (1 - d) * (z/pr.z)
  #w00.ps <- (1 - d) * ((1 - z) / (1 - pr.z))
  
  # Complier weights
  #w1c <- w11.ps - w10.ps
  #w0c <- (w01.ps - w00.ps)
  #kappa1 <- base::mean(w1c)
  #kappa0 <- base::mean(w0c)
  
  # Normalized complier weights
  #w1c <- w1c/kappa1
  #w0c <- w0c/kappa0
  #-----------------------------------------------------------------------------
  # Estimate LATE
  #mu.summand.Y1 <- w1c * y
  #mu.summand.Y0 <- w0c * y
  
  #mu.Y1 <- base::mean(mu.summand.Y1)
  #mu.Y0 <- base::mean(mu.summand.Y0)
 
  #late.ipw <- mu.Y1 - mu.Y0
  

  late.unf.mean <- mean((Y1-Y0)*d1)/mean(d1)
  late.unf.10 <- quantile(Y1[d1==1], probs = .1) - quantile(Y0[d1==1], probs = .1)
  late.unf.25 <- quantile(Y1[d1==1], probs = 0.25) - quantile(Y0[d1==1], probs = 0.25)
  late.unf.50 <- quantile(Y1[d1==1], probs = 0.5) - quantile(Y0[d1==1], probs = 0.5)
  late.unf.75 <- quantile(Y1[d1==1], probs = 0.75) - quantile(Y0[d1==1], probs = 0.75)
  late.unf.90 <- quantile(Y1[d1==1], probs = 0.9) - quantile(Y0[d1==1], probs = 0.9)
  #late <- mean((Y1-Y0)[d1==1])
  
  out <- cbind(late.unf.mean,late.unf.10, late.unf.25, late.unf.50, late.unf.75,
               late.unf.90)
  #----------------------------------------------------------------------------
  return(out)
}


n=1e7
nrep = 1000
ncores = 10
iseed = 1234#floor(1000*runif(1))
scale.ps <- 1
library(doParallel)
library(doRNG)

cl <- parallel::makeCluster(ncores) 
doSNOW::registerDoSNOW(cl)

a <- foreach::foreach(nn = 1:nrep, .options.RNG = iseed) %dorng% {
  late.est <- late(n, scale.ps)
  return(late.est)
}
#Stop the cluster
stopCluster(cl)

#Put the Monte Carlo Results in an Array
simu <- array(unlist(a), dim = c(nrow(a[[1]]), ncol(a[[1]]), length(a)))
#simu <- t(matrix(simu, 3, nrow = nrep, byrow = T))
simu2 <- t(matrix(as.vector(simu), nrow = 6, ncol = nrep))

m.simu = base::colMeans(simu2)

matrixStats::colSds(simu2)
