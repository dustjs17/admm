rm(list = ls())
if(!require(Deriv)){install.packages("Deriv")};library(Deriv)
if(!require(evd)){install.packages("evd")};library(evd)

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("admmgev.R")

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

#set value
n = 100
true_theta=c(30,0.1)
true_beta = rep(0,n)

y = sin(seq(0,2*pi,length = n))*10; mt = y +100  #sin
x= rep(0,n); set.seed(1)
for (i in 1:n) {
  x[i] = rgev(1,loc =mt[i],scale=true_theta[1],shape=true_theta[2])
}
z = diag(n)

# init value
dmatrix = mat_func(n)
rho = 0.5;lam = 2
z_init = rep(0,n-2) ; u_init = rep(1,n-2)
stepsize = 0.01 

start <- list()
start$scale <- sqrt(6 * var(x))/pi
start$loc <- mean(x) - 0.58 * start$scale #이거......
tvec = c(true_beta,0,0)
tvec[1:n] = start$loc  #init mu value
tvec[n+1] = start$scale #init sigma value


#update
for (iter in 1:500) {
  # cat("*****iter:::", iter, "\n")
  
  #mu, sigma, k optim
  old_tvec = gevreg(x = x,z = z)
  old_mu = old_tvec[1:n]
  
  # z update in ADMM
  tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init,lam = lam, rho = rho)
  
  # u update in ADMM
  tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
  
  u_init = tmp_u ; z_init = tmp_z
  tvec = old_tvec
}

plot(old_mu, type = "l",main = "optim")
points(x,col="red")
plot(x)

optim_mu=old_mu

data.frame(symbolic = new_mu, optim = optim_mu )


