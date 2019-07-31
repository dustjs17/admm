rm(list = ls())

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("2_update_function.R")

#init set
n = 1000
p = 0.99
sigma = 20
b = 0.5
##noise
set.seed(723)
zt = rnorm(n,0,sigma)
##trend slope
init_vt = 0 
vt = c(); vt = cbind(init_vt,vt)
set.seed(723)
for(i in 2:n){
  vt[i] = sample(c(vt[i-1],runif(1,-b,b)),1,prob = c(0.99,0.01))
}
##true trend
init_xt = 0
xt = c(); xt = cbind(init_xt,xt)
for (i in 2:n) {
  xt[i] = xt[i-1] + vt[i-1]
}
##time series data
yt = xt + zt  
plot(xt,type = "l")
plot(yt,type = "l")
##ADMM
init_z = c(rep(0,n-2)); init_u = c(rep(1,n-2))
rho = 0.1;lam = 5000
dmatrix = mat_func(n)
rhodt = rho*t(dmatrix); rhodtd = rhodt %*% dmatrix #beta update
imatrix = diag(1,n,n)
solve_term = solve(imatrix + rhodtd)


for (iter in 1:50000000) {
  # cat("*****iter:::", iter, "\n")
  
  tmp_beta = update_beta(imatrix = imatrix,solve_term = solve_term,
                         y = yt,rhodt = rhodt,z = init_z,u = init_u)
  
  tmp_z = update_z(dmatrix = dmatrix, beta=tmp_beta, u = init_u,lam = lam,rho = rho)
  
  tmp_u = update_u(u = init_u, dmatrix = dmatrix, beta = tmp_beta, z = tmp_z)
  
  init_u = tmp_u; init_z = tmp_z
}

plot(tmp_beta, type = "l",main = "optim-lambda=35000")





