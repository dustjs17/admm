## glmnet lasso
rm(list=ls())
if(!require(glmnet)){install.packages("glmnet"); library(glmnet)}

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("lassofunction.R")

set.seed(2018)
p = 10 ;n = 100
x =  matrix(rnorm(n*p), nrow = n) 
beta_true = rep(0,p); beta_true[c(2, 4)] = c(2, -1)
y =  x %*% beta_true + rnorm(n)


# init value
rho = 0.5
z_init =  c(rep(0,p)) ; u_init = c(rep(1,p))
lam = 100
rhomatrix = diag(rho, p, p)
xtxrho <- (t(x) %*% x) + rhomatrix ; invxtxrho <- solve(xtxrho)
xty <- t(x) %*% y

#admm update
i = 1
while(i <= 500){
  
  tmp_beta <- func_beta(z = z_init,u = u_init,
                        rho = rho, xty, invxtxrho)
  
  tmp_z <- func_z(u = u_init,beta = tmp_beta, lam = lam,
                  rho = rho)
  
  tmp_u <- func_u(beta = tmp_beta,z = tmp_z,u = u_init)
  
  u_init = tmp_u
  z_init = tmp_z
  
  i <- i+1
}

solution <- cbind(tmp_beta, tmp_z, tmp_u)

# GLMNET REF
fit_lasso = glmnet(x, y, standardize = F, intercept = F, lambda = (lam/n))
coef(fit_lasso)[-1, 1]

# plot(coef(fit_lasso)[-1,1], solution[,1] )
# abline(a = 0 , b = 1) 

coef(fit_lasso)[-1, 1]
solution[,2]



####################################
glmnet_yeon = function(x, y, rho, z_init, u_init, lam, max_iter = 1000)
{
  p = ncol(x)
  n = nrow(x)
  
  rhomatrix = diag(rho, p, p)
  xtxrho = (t(x) %*% x) + rhomatrix
  invxtxrho = solve(xtxrho)
  xty = t(x) %*% y
  
  for(i in 1:max_iter)
  {
    # beta 구하기
    tmp_beta = invxtxrho %*% (xty + rho * (z_init - u_init)) 
    
    # z 구하기
    tmp_z = ifelse(abs(u_init + tmp_beta) > (lam/rho), u_init + tmp_beta - sign(u_init + tmp_beta) * (lam /rho), 0)
    
    # u
    tmp_u = u_init + tmp_beta - tmp_z
    
    u_init = tmp_u
    z_init = tmp_z
    
    cat(i, '\n')
  }
  return(list(beta = tmp_beta, z = tmp_z, u = tmp_u))
}

fit_yeon = glmnet_yeon(x = x, y = y, rho = 0.0005, z_init = z_init, u_init = u_init, lam = cvlasso$lambda.min)


