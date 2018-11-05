rm(list = ls())
if(!require(Deriv)){install.packages("Deriv")};library(Deriv)
if(!require(evd)){install.packages("evd")};library(evd)

setwd("C:/Users/uos/Dropbox/admm")
source("library/admmgev.R")

# derivatives for mles
Grad = Deriv(logl,c("mu","sigma","k"),combine="cbind") ##partial derivative
Hmat = Deriv(logl,c("mu","sigma","k"),n=c(hessian=2),combine="cbind")


#constant- random variable
n = 50; true.vec = c(120,40,0.1) 
set.seed(1)
x= rgev(n,loc=true.vec[1] ,scale=true.vec[2],shape=true.vec[3]) 
fgev(x)

#increase - random variable
n = 50; true.vec = c(120,40,0.1)
y = seq(-10,10,length.out = n); mt = y +120
x= rep(0,n); set.seed(1)
for (i in 1:n) {
  x[i] = rgev(1,loc =mt[i],scale=true.vec[2],shape=true.vec[3])
}

fgev(x)

# init value
old_mu = rep(100,n)
old_sk = c(30,0.1) #(40,0.1)
stepsize = 0.01 
# stepsize = 0.01

dmatrix = mat_func(n)
rho = 0.5; lam = 5 
# lam = 100000
z_init = c(rep(0,n-2)) ; u_init = c(rep(1,n-2))
sec_term = rho * t(dmatrix)
rhodtd = sec_term %*% dmatrix
gvec_mu = rep(0,n)
init_hmat_mu = matrix(rep(0,n*n), ncol = n)


#update
for (iter in 1:10) {
  # cat("*****iter:::", iter, "\n")
  
  #mu, sigma, k iteration
  for (h in 1:100) {
    
    #mu iteration
    for(j in 1:10){
      mu = old_mu; sigma=old_sk[1]; k=old_sk[2]
      
      init_grad_mu = eval(Grad)[,1]  # gradient
      diag(init_hmat_mu) = eval(Hmat)[,1] 
      new_mu = old_mu -solve(init_hmat_mu + rhodtd) %*% (init_grad_mu + sec_term%*%(dmatrix%*%old_mu - z_init + u_init))
      
      #kkt condition check
      mu = new_mu
      grad3 = init_grad_mu + init_hmat_mu%*%(new_mu -old_mu) +sec_term%*%(dmatrix%*%new_mu - z_init + u_init) 
      
      old_mu = new_mu
      # plot(new_mu,main = "mu" )
    }
      
      # 1.mu check
    # if(h %%10 ==0){
      # cat("**sub_iter:::",h,"\n")
      # cat("mu:::", new_mu ,'\n')
      # cat("mu_gradient:::", grad3,"\n")
    # }
    
      #sigma, k iteration
      for (i in 1:100) {
        mu = old_mu; sigma=old_sk[1]; k=old_sk[2]
        gvec_sk = apply(eval(Grad),2,sum)[2:3] #gradient
        hmat_sk = matrix(apply(eval(Hmat),2,sum),3,3)[-1,-1] #hessian
        new_sk = old_sk - stepsize *solve(hmat_sk) %*% gvec_sk
        
        #kkt condition check
        sigma = new_sk[1] ; k = new_sk[2]; 
        grad2 = apply(eval(Grad),2,sum)[2:3]
      
        old_sk = new_sk
        
      }
      
      #2.theta check
    # if(h %% 10 ==0){
    #   cat("k:::", new_sk[2],"\n")
    #   cat("sigma:::", new_sk[1] ,'\n')
    #   cat("theta_gradient:::", grad2,'\n') 
    # }
    
    }
  
    #step1 kkt condition check
    grad = eval(Grad)[,1]
    sub_grad = sec_term%*%(dmatrix%*%new_mu - z_init + u_init)  ##### kkt 확인 조건식 다시 생각하기
    grad_mu = sub_grad + grad
    step1_grad = c(grad_mu ,grad2)
    
    # cat("step1_grad:::",step1_grad,'\n')
    
  # z update in ADMM
  tmp_z = func_z(dmatrix = dmatrix, u = u_init, mu = old_mu, lam = lam, rho = rho)
  
  #ste2 kkt condition check
  z_grad = rho*(tmp_z - dmatrix%*%old_mu -u_init)
  # cat("step2_grad:::",z_grad,'\n',"#######################################################", "\n")
  
  # u update in ADMM
  tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
  
  
  u_init = tmp_u ; z_init = tmp_z
  
}
plot(new_mu, type = "l")



