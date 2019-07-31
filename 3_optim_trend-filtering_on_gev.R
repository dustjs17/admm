rm(list = ls())

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("3_update_gev_function.R")

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

#set value
n = 100
true_theta=c(40,0.1)
true_beta = rep(0,n)
z = diag(n)

dmatrix = mat_func(n)

# resolution in ADMM
epri = 1e-4
edul = 1e-4
A = dmatrix
B = -diag(n-2)

# parameter used in optim()
ctr_list = list()
ctr_list$maxit = 20
ctr_list$reltol = 1e-6

# generate x in gev dist
## 1.stationay
mt = rep(100,n)
## 2. increasing
mt = seq(1,10,length.out = n)
## 3.sin
set.seed(1)
y = sin(seq(0,2*pi,length = n))*20; mt = y +100
## 4. sin + linear trend
set.seed(1)
y = sin(seq(0,2*pi,length = n))*20 + seq(1,n); mt = y +100

x= rep(0,n)
for (i in 1:n) {
  x[i] = rgev(1,loc =mt[i],scale=true_theta[1],shape=true_theta[2])
}

#model selection
init_rs = c()

lamb = seq(0.02,0.09,0.01)
# i = 3
par(mfrow=c(2,2))
#update
for (i in 1:length(lamb)) {
  lam = lamb[i]
  init_rho = 0.5
  z_init = rep(0,n-2) ; u_init = rep(1,n-2); init_mu = rep(Inf, n)
  
  start <- list()
  start$scale <- sqrt(6 * var(x))/pi
  start$loc <- mean(x) - 0.58 * start$scale 
  tvec = c(true_beta,0,0)
  tvec[1:n] = start$loc  #init mu value
  tvec[n+1] = start$scale #init sigma value
  tvec[n+2] = 0.1
  for (iter in 1:10000){
     cat("*****iter:::", iter, "\n")
    #mu, sigma, k optim
    old_tvec = gevreg(x = x,z = z,ctr_list)
    old_mu = old_tvec[1:n]
    
    # z update in ADMM
    tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init,lam = lam, rho = init_rho)
    
    AA = init_rho*t(A)%*%B
    sk = norm(drop(AA%*%(tmp_z -z_init)),"2")
    rk = norm(drop(A%*%old_mu + B%*%tmp_z),"2")
    
    if (rk<epri & sk <edul) break
    
    # u update in ADMM
    tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
    
    u_init = tmp_u ; z_init = tmp_z
    tvec = old_tvec 
    }
  
  plot(old_mu, type = "l",main = paste("lambda=",lam),ylim = c(30,170))
  points(x,col="red")
  
  #model selection
  dmu = dmatrix %*% old_mu
  stra = dmu - tmp_z
  nonzero = length(which(abs(dmu) > 1e-1))
  logl = - sum(dgev(x, loc = old_tvec[1:n], scale = old_tvec[n+1], shape = old_tvec[n+2], log = T) )
  
  aic = 2*logl + 2*nonzero
  bic = 2*logl + log(n)*nonzero
  tmp_rs = cbind(nonzero,aic,bic)
  init_rs = rbind(init_rs,tmp_rs)
}

cbind(lamb,init_rs)







