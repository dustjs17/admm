rm(list = ls())
if(!require(Deriv)){install.packages("Deriv")};library(Deriv)
if(!require(evd)){install.packages("evd")};library(evd)

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("admmgev.R")

# derivatives for mles
Gmu = Deriv(logl,"mu"); Hmu = Deriv(logl,"mu",n=c(hessian=2))
Gtheta = Deriv(logl,c("sigma","k"),combine = "cbind") 
Htheta = Deriv(logl,c("sigma","k"),n=c(hessian=2),combine = "cbind")

n= 100;true.vec = c(100,30,0.1)
#1. constant situation
set.seed(1)
x= rgev(n,loc=true.vec[1] ,scale=true.vec[2],shape=true.vec[3])

#2. linear
y = seq(-10,10,length.out = n); mt = y +100
# plot(y)
x= rep(0,n); set.seed(1)
for (i in 1:n) {
  x[i] = rgev(1,loc =mt[i],scale=true.vec[2],shape=true.vec[3])
}

#3. sin
y = sin(seq(0,2*pi,length = n))*10
plot(y)
mt = y +100

x= rep(0,n); set.seed(1)
for (i in 1:n) {
  x[i] = rgev(1,loc =mt[i],scale=true.vec[2],shape=true.vec[3])
}

#4. quadratic
y = (seq(1,100,length = n)-50)^2/250
plot(y)
mt = y + 100
x= rep(0,n); set.seed(1)
for (i in 1:n) {
  x[i] = rgev(1,loc =mt[i],scale=true.vec[2],shape=true.vec[3])
}

#5. piecewise
y = c(rep(1,20),seq(1,10,by = 0.5),rep(11,20),seq(12,22,by =0.5),rep(23,20))
plot(y)
mt = y + 100
x= rep(0,n); set.seed(1)
for (i in 1:n) {
  x[i] = rgev(1,loc =mt[i],scale=true.vec[2],shape=true.vec[3])
}


# fgev(x)

# init value
old_mu = mt 
old_sk = c(30,0.1) #(40,0.1)
stepsize = 0.01 

dmatrix = mat_func(n)
rho = 0.5; lam = 2 ; tol = 1e-08
z_init = rep(0,n-2) ; u_init = rep(1,n-2)
sec_term = rho * t(dmatrix); rhodtd = sec_term %*% dmatrix
init_hmat_mu = matrix(0,n,n)


#update
for (iter in 1:500) {
  # cat("*****iter:::", iter, "\n")
  
  #mu, sigma, k iteration
  for (h in 1:10) {
    
    #mu iteration
    for(j in 1:10){
      mu = old_mu; sigma=old_sk[1]; k=old_sk[2]
      init_grad_mu = eval(Gmu)  
      diag(init_hmat_mu) = eval(Hmu) 
      new_mu = old_mu - tryCatch(solve(init_hmat_mu+rhodtd),
                                 error = function(e){solve(init_hmat_mu+ rhodtd + diag(tol,nrow(init_hmat_mu)))})%*% 
        (init_grad_mu + sec_term %*%(dmatrix%*%old_mu - z_init + u_init))
     
       #kkt condition check
      # mu = new_mu
      # grad3 = init_grad_mu + init_hmat_mu%*%(new_mu -old_mu) +sec_term%*%(dmatrix%*%new_mu - z_init + u_init)
      
      old_mu = new_mu
      
    }
    # 1.mu check
    # if(h %%10 ==0){
    # cat("**sub_iter:::",h,"\n")
    # cat("mu:::", new_mu ,'\n')
    # cat("mu_gradient:::", grad3,"\n")
    # }
    
    #sigma, k iteration
    for (i in 1:10) {
      mu = old_mu; sigma=old_sk[1]; k=old_sk[2]
      gvec_sk = apply(eval(Gtheta),2,sum) 
      hmat_sk = matrix(apply(eval(Htheta),2,sum),2,2) 
      new_sk = old_sk - stepsize *tryCatch(solve(hmat_sk),
                                           error = function(e){solve(hmat_sk +diag(tol,nrow(hmat_sk)))}) %*% gvec_sk

      #kkt condition check
      # sigma = new_sk[1] ; k = new_sk[2]; 
      # grad2 = apply(eval(Grad),2,sum)[2:3]
      
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
  # grad = eval(Grad)[,1]
  # sub_grad = sec_term%*%(dmatrix%*%new_mu - z_init + u_init)  ##### kkt 확인 조건식 다시 생각하기
  # grad_mu = sub_grad + grad
  # step1_grad = c(grad_mu ,grad2)
  
  # cat("step1_grad:::",step1_grad,'\n')
  
  # z update in ADMM
  tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init,lam = lam, rho = rho)
  
  # cat("constraint::",min(dmatrix%*%old_mu - tmp_z),max(dmatrix%*%old_mu - tmp_z),'\n')
  # ste2 kkt condition check
  # z_grad = rho*(tmp_z - dmatrix%*%old_mu -u_init)
  # cat("step2_grad:::",z_grad,'\n',"#######################################################", "\n")
  
  # u update in ADMM
  tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
  
  u_init = tmp_u ; z_init = tmp_z
  
}

plot(new_mu, type = "l",main = "symbolic")
points(x,col="red")
plot(x)



