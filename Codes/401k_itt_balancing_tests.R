library(here)
source(here("Codes", "dist_balance.R"))


treat = data401k$e401
#GLM
w_glm1 <- treat/ps.glm$fitted.values
w_glm0 <- (1-treat)/(1 - ps.glm$fitted.values)
w_glm1 <- w_glm1/mean(w_glm1)
w_glm0 <- w_glm0/mean(w_glm0)
w_glm <- w_glm1 - w_glm0

#CBps - just identified
w_cbps1 <- treat/ps.cbps$fitted.values
w_cbps0 <- (1-treat)/(1 - ps.cbps$fitted.values)
w_cbps1 <- w_cbps1/mean(w_cbps1)
w_cbps0 <- w_cbps0/mean(w_cbps0)
w_cbps <- w_cbps1 - w_cbps0


#CBPS-over identified
w_cbps21 <- treat/ps.cbps.over$fitted.values
w_cbps20 <- (1-treat)/(1 - ps.cbps.over$fitted.values)
w_cbps21 <- w_cbps21/mean(w_cbps21)
w_cbps20 <- w_cbps20/mean(w_cbps20)
w_cbps2 <- w_cbps21 - w_cbps20


#IPS-exp
w_exp1 <- treat/ps.exp$fitted.values
w_exp0 <- (1-treat)/(1 - ps.exp$fitted.values)
w_exp1 <- w_exp1/mean(w_exp1)
w_exp0 <- w_exp0/mean(w_exp0)
w_exp <- w_exp1 - w_exp0


#IPS-Proj
w_proj1 <- treat/ps.proj$fitted.values
w_proj0 <- (1-treat)/(1 - ps.proj$fitted.values)
w_proj1 <- w_proj1/mean(w_proj1)
w_proj0 <- w_proj0/mean(w_proj0)
w_proj <- w_proj1 - w_proj0

#IPS-Ind
w_ind1 <- treat/ps.ind$fitted.values
w_ind0 <- (1-treat)/(1 - ps.ind$fitted.values)
w_ind1 <- w_ind1/mean(w_ind1)
w_ind0 <- w_ind0/mean(w_ind0)
w_ind <- w_ind1 - w_ind0



#xcov = baselineX[,-1]

covariates_names = c("inc", "loginc" , "age" ,"fsize" , "educ",
                "hown" ,
                "marr", "twoearn", "db", "pira")
xcov = baselineX[,covariates_names]

db_glm <- dist_balance(covs = xcov, weights = w_glm, chunk = 10000)
db_cbps <- dist_balance(covs = xcov, weights = w_cbps, chunk = 10000)
db_cbps2 <- dist_balance(covs = xcov, weights = w_cbps2, chunk = 10000)
db_exp <- dist_balance(covs = xcov, weights = w_exp, chunk = 10000)
db_proj <- dist_balance(covs = xcov, weights = w_proj, chunk = 10000)
db_ind <- dist_balance(covs = xcov, weights = w_ind, chunk = 10000)
