rm(list = ls())

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("3_update_gev_function.R")
if(!require(InspectChangepoint)){install.packages("InspectChangepoint")};library(InspectChangepoint)

n = 100
true_vec = c(100,30,0.1)

#sin function
set.seed(1)
y = sin(seq(0,2*pi,length.out = n))*10; mt = y +100
x = rep(0,n)
for (i in 1:n) {
  # i = 1
  x[i] = rgev(1,loc = mt[i],scale = true_vec[2],shape = true_vec[3]) 
}


# init value
old_mu = mt 
old_sk = c(30,0.1) #(40,0.1)
stepsize = 0.01 

dmatrix = mat_func(n)
rho = 0.5; lam = 0.1 ; tol = 1e-08
z_init = rep(0,n-2) ; u_init = rep(1,n-2)
init_hmat_mu = matrix(0,n,n)

sec_term = rho * t(dmatrix); rhodtd = sec_term %*% dmatrix

# derivatives for mles
Gmu = Deriv(logl,"mu"); Hmu = Deriv(logl,"mu",n=c(hessian=2))
Gtheta = Deriv(logl,c("sigma","k"),combine = "cbind") 
Htheta = Deriv(logl,c("sigma","k"),n=c(hessian=2),combine = "cbind")

#update
for (iter in 1:400) {
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
      
      old_mu = new_mu
      
    }
    
    #sigma, k iteration
    for (i in 1:10) {
      mu = old_mu; sigma=old_sk[1]; k=old_sk[2]
      gvec_sk = apply(eval(Gtheta),2,sum) 
      hmat_sk = matrix(apply(eval(Htheta),2,sum),2,2) 
      new_sk = old_sk - stepsize *tryCatch(solve(hmat_sk),
                                           error = function(e){solve(hmat_sk +diag(tol,nrow(hmat_sk)))}) %*% gvec_sk
      
      old_sk = new_sk
    }
    
  }
 
  # kkt condition,loss value check - mu,theta
  # if(iter %% 10 == 0){
  #   mu = old_mu; sigma = old_sk[1]; k = old_sk[2]
  # 
  #   # kkt - mu check
  #   kkt_mu = eval(Gmu) + sec_term %*% (dmatrix %*% mu -z_init + u_init)
  #   # kkt - theta check
  #   kkt_theta = apply(eval(Gtheta),2,sum)
  #   # loss value 
  #   first_loss = apply(func_gev(x = x, sigma = sigma,k = k, mu = mu),2,sum)
  #   third_loss = apply(func_thr(rho = rho,d = dmatrix,mu = mu,z = z_init,u = u_init),
  #                      2, sum)
  #   loss = first_loss + third_loss
  #   
  #   
  #   cat("**iter::",iter,"\n")
  #   cat("gev_loss::", loss, "\n")
  #   # cat("current_mu::", max(mu), min(mu),"\n")
  #   cat("mu_gradient::", max(kkt_mu), min(kkt_mu),"\n")
  #   cat("theta_gradient::", max(kkt_theta), min(kkt_theta),'\n')
  # }
  
  # z update in ADMM
  tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init,lam = lam, rho = rho)
  
  # kkt condition - z check
  # if(iter %% 10 ==0){
  #   
  #   kkt_z = rho * (tmp_z - (dmatrix %*% old_mu) - u_init)
  #   z_loss = lam*sum(abs(tmp_z))
  #   third_loss = apply(func_thr(rho = rho,d = dmatrix,mu = old_mu,z = tmp_z,u = u_init),
  #                      2,sum)
  #   # sec_loss = z_loss + third_loss  
  #   cat("number of nonzero::",length(which(!(tmp_z == 0))),'\n')
  #   cat("z_value::",max(tmp_z),min(tmp_z),"\n")
  #   cat("z_loss::",z_loss,'\n')
  #   cat("third_loss::",third_loss,"\n")
  #   cat("z_gradient::",max(kkt_z),min(kkt_z),'\n',"#######################################","\n")
  #   
  #   
  #   plot(new_mu,type = "l",main = paste0("iter:",iter))
  # }
  
  # u update in ADMM
  tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
  
  u_init = tmp_u ; z_init = tmp_z
  
}

# plot(new_mu,type = "l",main = paste0(" symbolic ,",  "lambda = ",lam))
symbolic = new_mu
symbolic_sk = old_sk

#condition check
(dmatrix%*% new_mu)- diag(1,(n-2),(n-2)) %*% z_init
