rm(list = ls())
setwd("C:/Users/UOS/Documents/GitHub/admm")
source("3_update_gev_function_.R")
if(!require(dplyr)){install.packages("dplyr")};require(dplyr)

load("C:/Users/UOS/서울시립대학교/전종준 - lab_work/Lab_process/HIJ/kma_data/Pr_46.RData")
str(Pr_46)
st = distinct(Pr_46,stnlds)

par(mfrow=c(2,2))

for (i in 1:nrow(st)) {
  # i = 1
  x = Pr_46 %>% subset(stnlds == st[i,])
  x = x$pr
  
  # optim control
  optim_controlList = list()
  optim_controlList$maxit = 1e+3
  
  #init_value
  n = length(x)
  z = diag(n)
  # plot(x)
  true_beta = rep(0,n)
  
  dmatrix = mat_func(n)
  rho = 0.5;lam = 0.02
  z_init = rep(0,n-2) ; u_init = rep(1,n-2)
  stepsize = 0.01 
  
  start <- list()
  start$scale <- sqrt(6 * var(x))/pi
  start$loc <- mean(x) - 0.58 * start$scale #이거......
  tvec = c(true_beta,0,0)
  tvec[1:n] = start$loc  #init mu value
  tvec[n+1] = start$scale
  
  #update
  s_time = Sys.time()
  for (iter in 1:300) {
    
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
  e_time = Sys.time()
  dif = e_time-s_time
  # cat("지역:::",st[i,] ,dif,"\n")
  
  plot(old_mu, type = "l",main = paste("지역:",st[i,],", lamda=",lam))
  
}



