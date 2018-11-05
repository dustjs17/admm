logl=expression(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k))
n =100
#calculate loss
loss = matrix(0,n,n)
loss_func = function(x, mu , sigma , k){
  for (i in 1:n ) {
    loss[,i] = (log(sigma)+(1+1/k)*(log(1+k*(x[i]-mu)/sigma))+((1+k*(x[i]-mu)/sigma)^(-1/k)))
  }
  gev_loss = apply(loss, 1, sum) ##row = 1
  return(gev_loss)
}

#z_function
func_z <- function(dmatrix,mu,u,lam,rho)
{
  z = ifelse(abs(dmatrix %*% mu + u) > (lam/rho) ,
             (dmatrix %*% mu) + u - sign(u + (dmatrix %*% mu)) * (lam /rho), 0)
  return(z)
}

#u_function
func_u <- function(dmatrix,mu,z,u) 
{
  u <- u + (dmatrix %*% mu) - z
  return(u)
}

#dmatrix
mat_func = function(n) {
  m = matrix(0,n-2,n)   ## dmatrix : (n-2)*n 
  for (i in 2:n-2) {
    m[i,i] = 1
    m[i,i+1] = -2
    m[i,i+2] = 1
  }
  return(m)
}

