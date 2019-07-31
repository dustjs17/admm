rm(list = ls())
setwd("C:/Users/UOS/Documents/GitHub/admm")
source("3_update_gev_function.R")
if(!require(dplyr)){install.packages("dplyr")};require(dplyr)
if(!require(evd)){install.packages("evd")};require(evd)

load("C:/Users/UOS/서울시립대학교/전종준 - lab_work/Lab_process/HIJ/kma_data/Pr_46.RData")
st = distinct(Pr_46,stnlds)

# parameter used in optim()
ctr_list = list()
ctr_list$maxit = 20
ctr_list$reltol = 1e-6

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

# resolution in ADMM
epri = 1e-4
edul = 1e-4

par(mfrow=c(2,2))

i = 1
# for (i in 1:nrow(st)) {
  x = Pr_46 %>% subset(stnlds == st[i,])
  x = x$pr

  #init_value
  n = length(x)
  z = diag(n)
  true_beta = rep(0,n)
  
  dmatrix = mat_func(n)
  init_rho = 0.5;lam = 0.1
  z_init = rep(0,n-2) ; u_init = rep(1,n-2)
  
  A = dmatrix
  B = -diag(n-2)
  
  start <- list()
  start$scale <- sqrt(6 * var(x))/pi
  start$loc <- mean(x) - 0.58 * start$scale 
  tvec = c(true_beta,0,0)
  tvec[1:n] = start$loc
  tvec[n+1] = start$scale
  tvec[n+2] = 0.1
  
  #update
  # s_time = Sys.time()
  for (iter in 1:10000) {
    
    #mu, sigma, k optim
    old_tvec = gevreg(x = x,z = z,ctr_list = ctr_list)
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
    
    cat("iter::",iter,"\n")
    
  }
  # e_time = Sys.time()
  # dif = e_time-s_time
  # cat("지역:::",st[i,] ,dif,"\n")
  
  plot(old_mu, type = "l",main = paste("지역:",st[i,],", lamda=",lam))
  points(x,col = "red")
  
# }



