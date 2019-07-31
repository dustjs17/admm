rm(list = ls())

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("3_update_gev_function.R")

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

#init set
n = 100
p = 0.99
sigma = 2
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
yt = xt + zt + 20  
##gev set
true_theta=c(30,0.1)
true_mu = rep(0,n)
z = diag(n)

dmatrix = mat_func(n)
rho = 0.5;lam = 5
z_init = rep(0,n-2) ; u_init = rep(1,n-2)

start = list()
start$scale = sqrt(6 * var(yt))/pi
start$loc = mean(yt) - 0.58 * start$scale #이거......
tvec = c(true_mu,0,0)
tvec[1:n] = start$loc  #init mu value
tvec[n+1] = start$scale #init sigma value

for (iter in 1:10000) {
  # cat("*****iter:::", iter, "\n")
  
  #mu, sigma, k optim
  old_tvec = gevreg(x = yt,z = z)$par
  old_mu = old_tvec[1:n]
  
  # z update in ADMM
  tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init,lam = lam, rho = rho)
  
  # u update in ADMM
  tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
  
  u_init = tmp_u ; z_init = tmp_z; tvec = old_tvec
}

plot(old_mu, type = "l")




