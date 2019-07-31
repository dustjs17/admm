#lasso update function
setwd("C:/Users/UOS/Documents/GitHub/admm")

#beta function
func_beta <- function(z,u,rho,xty, invxtxrho) {
  beta <- invxtxrho %*% (xty + rho * (z-u)) 
  return(beta)
}

#z function
func_z <- function(beta,u,lam,rho) {
  z <- ifelse(abs(u +beta) > (lam/rho) , 
              u + beta - sign(u + beta) * (lam /rho), 0)
  return(z)
}

#u function
func_u <- function(beta,z,u) {
  u <- u + (beta-z)
  return(u)
}

#######################trend-filtering##########################
#dmatrix
mat_func = function(n) {
  m = matrix(0,n-2,n)   ## dmatrix : (n-2)*n 
  for (i in 1:(n-2)) {
    m[i,i] = 1
    m[i,i+1] = -2
    m[i,i+2] = 1
  }
  return(m)
}

#update function
update_beta = function(imatrix,solve_term,y,rhodt,z,u){
  beta = solve_term  %*% (y + rhodt %*% (z-u))
  return(beta)
}

update_z = function(dmatrix,beta,u,rho,lam){
  z = ifelse(abs(dmatrix%*%beta + u)>lam/rho,
             dmatrix%*%beta + u - sign(dmatrix%*%beta + u)*lam/rho,0)
  return(z)
}


update_u = function(u,dmatrix,beta,z){
  u = u + dmatrix%*%beta - z
  return(u)
}
