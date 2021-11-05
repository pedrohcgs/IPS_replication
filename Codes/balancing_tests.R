source(here("Codes", "dist_balance.R"))


#GLM
w_glm1 <- treat/fit.glm$fitted.values
w_glm0 <- (1-treat)/(1 - fit.glm$fitted.values)
w_glm1 <- w_glm1/mean(w_glm1)
w_glm0 <- w_glm0/mean(w_glm0)
w_glm <- w_glm1 - w_glm0

#CBps - just identified
w_cbps1 <- treat/fit.cbps$fitted.values
w_cbps0 <- (1-treat)/(1 - fit.cbps$fitted.values)
w_cbps1 <- w_cbps1/mean(w_cbps1)
w_cbps0 <- w_cbps0/mean(w_cbps0)
w_cbps <- w_cbps1 - w_cbps0


#CBPS-over identified
w_cbps21 <- treat/fit.cbps2$fitted.values
w_cbps20 <- (1-treat)/(1 - fit.cbps2$fitted.values)
w_cbps21 <- w_cbps21/mean(w_cbps21)
w_cbps20 <- w_cbps20/mean(w_cbps20)
w_cbps2 <- w_cbps21 - w_cbps20


#IPS-exp
w_exp1 <- treat/fit.exp$fitted.values
w_exp0 <- (1-treat)/(1 - fit.exp$fitted.values)
w_exp1 <- w_exp1/mean(w_exp1)
w_exp0 <- w_exp0/mean(w_exp0)
w_exp <- w_exp1 - w_exp0


#IPS-Proj
w_proj1 <- treat/fit.proj$fitted.values
w_proj0 <- (1-treat)/(1 - fit.proj$fitted.values)
w_proj1 <- w_proj1/mean(w_proj1)
w_proj0 <- w_proj0/mean(w_proj0)
w_proj <- w_proj1 - w_proj0

#IPS-Ind
w_ind1 <- treat/fit.ind$fitted.values
w_ind0 <- (1-treat)/(1 - fit.ind$fitted.values)
w_ind1 <- w_ind1/mean(w_ind1)
w_ind0 <- w_ind0/mean(w_ind0)
w_ind <- w_ind1 - w_ind0


db_glm <- dist_balance(covs = xcov, weights = w_glm)
db_cbps <- dist_balance(covs = xcov, weights = w_cbps)
db_cbps2 <- dist_balance(covs = xcov, weights = w_cbps2)
db_exp <- dist_balance(covs = xcov, weights = w_exp)
db_proj <- dist_balance(covs = xcov, weights = w_proj)
db_ind <- dist_balance(covs = xcov, weights = w_ind)
