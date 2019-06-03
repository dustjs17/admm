rm(list = ls())

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("update_function.R")

set.seed(2018)
n = 100

#init value
sol_beta = seq(0,100,length.out = n) # beta : n*1
y = sol_beta + rnorm(n,0,1)  ## y : n*1
z_init = c(rep(0,n-2)); u_init = c(rep(1,n-2))
rho = 0.5; lam = 500

dmatrix = mat_func(n); imatrix = diag(1, n, n)
rhodd = rho * t(dmatrix) %*% dmatrix #beta update때 필요
first_term = solve(imatrix + rhodd)
sec_term = rho * t(dmatrix) #beta update때 필요



i = 1
while(i <= 100)
{
  
  tmp_beta = tf_func_beta(z = z_init,u = u_init,
                       rho = rho, y = y , first_term = first_term, sec_term = sec_term)
  
  tmp_z = tf_func_z(dmatrix = dmatrix, u = u_init, beta = tmp_beta, lam = lam,
                 rho = rho)
  
  tmp_u = tf_func_u(dmatrix = dmatrix ,beta = tmp_beta, z = tmp_z, u = u_init)
  
  u_init = tmp_u
  z_init = tmp_z
  
  i <- i+1
}

solution = list(tmp_beta , tmp_z , tmp_u)
solution[[1]]

plot(y)
points(1:100,solution[[1]],col = "red")

#############################################################
#genlasso
if(!require(devtools)){install.packages("devtools")} ; require(devtools)
assignInNamespace("version_info", 
                  c(devtools:::version_info, 
                    list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), 
                  "devtools")
install.packages("pkgbuild")
library(pkgbuild)
install.packages("devtools")
library(devtools)
devtools::install_github("r-lib/pkgbuild")
find_rtools(debug =  T)
devtools::install_github("glmgen/genlasso")
require(genlasso)


x = diag(1,n,n) ## x = I
out = genlasso(y,x,dmatrix) 
gen_beta=coef(out, lambda = 5)$beta
plot(solution[[1]])
points(1:100, gen_beta, col = "red")


