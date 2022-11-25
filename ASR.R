m = c(3794901,15191619,19919840,20056779,19819518, 18257225,17722067,19511370,
      22179956,22479229,19805793,17224356,13307234,10654272,9409940,8725574,7414559,
      4900234,4259173)
f = 1/sum(m)*m
J= length(f)

mass_func = function(Nj, rho){
  sigma = rho*Nj
  lower = floor(Nj - 3 * sigma)
  upper = ceiling(Nj + 3 *sigma)
  p1 = pnorm(( 0.5 - 3* sigma) / sigma)
  int = (lower+ 1): (upper - 1)
  p2 = pnorm( (int + 0.5 - Nj)/ sigma) - 
    pnorm( (int - 0.5 - Nj)/ sigma)
  p3 = 1 - pnorm((3* sigma - 0.5) / sigma)
  p =c(p1,p2,p3)
  result = list(p, lower, upper)
  names(result) = c("p","lower","upper")
  return(result)
}

##Calculate \hat{N} and \hat{var}(\hat{N})
##Nj is the true population size in j-th group
ControlCV = function(Nj, rho){
  J=length(Nj)
  N =rep(0, J)
  V_N =rep(0, J)
  dense = lapply(Nj, mass_func, rho = rho)
  for (j in 1: J){
    N[j] = sample(dense[[j]]$lower: dense[[j]]$upper, 
                  1, replace = TRUE, prob = dense[[j]]$p)
    emp = mass_func(N[j], rho)
    k = emp$lower: emp$upper
    p = emp$p
    V_N[j]= sum(k^2 * p) - (sum( k * p))^2
  }
  result = list(N, V_N)
  names(result) = c("N", "V_N")
  return(result)
}

###Calculate %RB, var, CV of an estimator
performance = function(est, true){
  perRB_est = 100 * (mean(est)- true)/true
  var_est= var(est)
  CV_est=sqrt( var_est)/mean(est)
  result = c(perRB_est, var_est, CV_est)
  names(result) = c("RB","Var", "CV")
  return(result)
}


##########################Estimate rate
##rho is a scaler
rate = function(lambda, Nj, rho){
  Sample_est = ControlCV(Nj, rho)
  N = Sample_est$N
  X = rpois(J, lambda * Nj)
  Rj = X / N
  R = sum(f * Rj)###get \hat{R}
  
  V_N = Sample_est$V_N
  Rj_bc = Rj * (1- V_N/ N^2)
  R_bc = sum(f * Rj_bc)###get \hat{R}_{bc}
  
  bj= V_N / N^2
  var_Rbc = sum( f^2 *Rj_bc* (1 - bj) ^ 2 * (Rj_bc *bj + 1/N + 3 /N *bj - Rj_bc * bj^2))
  ###get variance of \hat{R}_{bc}
  
  result =c(R, R_bc,var_Rbc)
  names(result) = c("R", "R_bc","var_Rbc")
  return(result)
}



bias_performance<-function(Np, lambda, rho, K){
  Nj = floor( Np * f)
  Nj[J] = Np - sum(Nj[1:(J-1)])
  for (i in 1: length(rho)){
  rate= replicate(K, rate(lambda, Nj, rho[i]))
  R.per= performance(rate[1,], lambda)
  RB_R = R.per[1]
  
  Rbc.per= performance(rate[2,], lambda)
  RB_Rbc = Rbc.per[1]
  var_Rbc_t = Rbc.per[2]
 
  var.per= performance(rate[3,], var_Rbc_t)
  Evar = mean(rate[3,])
  RB_var = var.per[1]
  CV_var = var.per[3]
  per = c(Np, lambda, rho[i], RB_R, RB_Rbc, var_Rbc_t, Evar, RB_var, CV_var)
  write.table(t(per), file = "rate.csv", sep=",", na = "",
              row.names=FALSE,col.names=FALSE,append = TRUE)
  } 
}

setwd("~/ratio1")
rho = seq(0.05,0.3, by = 0.05)
K= 10000
bias_performance(100000, 0.001, rho, K)
bias_performance(50000, 0.001, rho, K)
bias_performance(10000, 0.001, rho, K)

bias_performance(100000, 1/5000, rho, K)
bias_performance(50000, 1/5000, rho, K)
bias_performance(10000, 1/5000, rho, K)

bias_performance(100000, 0.0001, rho, K)
bias_performance(100000, 0.00005, rho, K)
