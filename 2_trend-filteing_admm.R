rm(list = ls())

setwd("C:/Users/UOS/Documents/GitHub/admm")
source("2_update_function.R")


#init value
# sol_beta = seq(0,100,length.out = n) # beta : n*1
#set value
n = 100
sol_beta = sin(seq(0,2*pi,length = n))*10
set.seed(2019)
y = sol_beta + rnorm(n,0,1)
init_z = c(rep(0,n-2)); init_u = c(rep(1,n-2))

rho = 0.5;lam = 100
dmatrix = mat_func(n)
rhodt = rho*t(dmatrix); rhodtd = rhodt %*% dmatrix #beta update
imatrix = diag(1,n,n)
solve_term = solve(imatrix + rhodtd)


for(iter in 1:500000){
  
  tmp_beta = update_beta(imatrix = imatrix,solve_term = solve_term,
                         y = y,rhodt = rhodt,z = init_z,u = init_u)
  
  tmp_z = update_z(dmatrix = dmatrix, beta=tmp_beta, u = init_u,lam = lam,rho = rho)
  
  tmp_u = update_u(u = init_u, dmatrix = dmatrix, beta = tmp_beta, z = tmp_z)
  
  init_u = tmp_u; init_z = tmp_z
  
}

plot(tmp_beta,type = "l")



#############################################################
#genlasso
if(!require(devtools)){install.packages("devtools")} ; require(devtools)
if(!require(processx)){install.packages("processx")} ; require(processx)
assignInNamespace("version_info", 
                  c(devtools:::version_info, 
                    list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), 
                  "devtools")
install.packages("pkgbuild")
library(pkgbuild)
install.packages("devtools")
library(devtools)
devtools::install_github("r-lib/pkgbuild")
find_rtools(debug =  T)
devtools::install_github("glmgen/genlasso",force = TRUE)
require(genlasso)


x = diag(1,n,n) ## x = I
out = genlasso(y,x,dmatrix) 
gen_beta=coef(out, lambda = lam)$beta
points(1:100, gen_beta, col = "red",pch = 3)
plot(gen_beta)

gen_beta - tmp_beta
