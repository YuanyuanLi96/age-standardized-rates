m = c(3794901,15191619,19919840,20056779,19819518, 18257225,17722067,19511370,
      22179956,22479229,19805793,17224356,13307234,10654272,9409940,8725574,7414559,
      4900234,4259173)
f = 1/sum(m)*m
J= length(f)


##################Function for hypothesis testing
testing = function(alpha, esti, v_esti){
  crit = qnorm( 1 - alpha)
  crit_two = qnorm( 1 - alpha/2)
  t_stat = (esti- 1)/ sqrt(v_esti)
  result = list()
  result$t_stat = t_stat
  result$prob_reject = sapply(crit, function(x)sum(t_stat > x)/ length(t_stat))
  result$ci_prob= sapply(crit_two, function(x)sum(-x <= t_stat & t_stat <= x)/ 
                           length(t_stat))
  return(result)
} 


################Three methods for estimating N, V_N: SRS, STR, ControlCV
#Performance relate to sample size n
SRS = function(Pop, rho){
  ybar = sapply(1:J, function(x)
    sum(sample(Pop$a,Pop$n,replace=FALSE) == x)/Pop$n)
  N = Pop$Np * ybar
  if (all(N>0)==F){
    cat("Not all N is positive!","\n")
    cat("N: ", N,"\n")
    cat("ybar: ", ybar,"\n")
  }
  V_N = Pop$Np * N * (Pop$Np-N) / (Pop$Np-1) * ( 1/Pop$n - 1/Pop$Np)
  result = list(N, V_N, ybar)
  names(result) = c("N", "V_N","ybar")
  return(result)
}

#Performance relate to sample size n
STR = function(Pop, rho){
  ph_hat = matrix(0, J, H)
  for(h in 1:H){
    samp = sample(Pop$Ph[[h]], Pop$nh[h], replace = FALSE)
    ph_hat[,h] = sapply(1:J, function(x)sum(samp == x)/Pop$nh[h])
  }
  N = as.vector(ph_hat %*% Pop$Nst)
  if (all(N>0)==F){
    cat("Not all N is positive!", "\n")
    cat("N: ", N,"\n")
    print(ph_hat)
  }
  V_N =  as.vector((ph_hat * (1- ph_hat)) %*%
                     ( 1/ (Pop$nh-1) * (1 - Pop$nh/Pop$Nst) * Pop$Nst^2))
  result = list(N, V_N)
  names(result) = c("N", "V_N")
  return(result)
}


#Performance do not relate to sample size n
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



ControlCV = function(Pop, rho){
  J=length(Pop$Nj)
  N =rep(0, J)
  V_N =rep(0, J)
  dense = lapply(Pop$Nj, mass_func, rho = rho)
  for (j in 1: J){
    N[j] = sample(dense[[j]]$lower: dense[[j]]$upper, 
                  1, replace = TRUE, prob = dense[[j]]$p)
    emp = mass_func(N[j], rho)
    k = emp$lower: emp$upper
    p = emp$p
    V_N[j]= sum(k^2 * p) - (sum( k * p))^2
  }
  if (all(N>0)==F){
    cat("Not all N is positive!", "\n")
    cat("N: ", N,"\n")
    print(ph_hat)
  }
  result = list(N, V_N)
  names(result) = c("N", "V_N")
  return(result)
}




##########################Estimate rate
rate = function(lambda, method, Pop, rho){
  X = rpois(J, lambda * Pop$Nj)
  Sample_est = method(Pop, rho)
  N = Sample_est$N
  Rj = X / N
  R = sum(f * Rj)
  V_N = Sample_est$V_N
  Rj_bc = Rj * (1- V_N/ N^2)
  bj= V_N / N^2
  R_bc = sum(f * Rj_bc)
  uj =  Rj_bc*(1-bj)^2 *(Rj_bc* bj +1/ N + 3 /N *bj-  Rj_bc * bj^2)
  
  b = sum(f *bj * Rj)
  g = sum(f^2 *(Rj / N + Rj^2 / N^2 * V_N))
  var_Rbc = sum( f^2 * uj)
  result =c(R, b, g, R_bc,var_Rbc)
  names(result) = c("R", "b", "g", "R_bc","var_Rbc")
  return(result)
}


###################Estimate ratio
###Calculate %RB, var, CV of an estimator
performance = function(est, true){
  perRB_est = 100 * (mean(est)- true)/true
  var_est= var(est)
  CV_est=sqrt( var_est)/mean(est)
  result = c(perRB_est, var_est, CV_est)
  names(result) = c("%RB","Var", "CV")
  return(result)
}

##lambda1 is population rate for the numerator, lambda2 is population rate for the denominator 
ratio = function(lambda1,lambda2, K, alpha , method, P1, P2, rho1,rho2, table_name){
  rate1= replicate(K, rate(lambda1, method, P1, rho1))
  rate2= replicate(K, rate(lambda2, method, P2, rho2))
  R1 = rate1["R",]
  b = rate1["b",]
  hat_f= rate1["g",]
  R1_bc = rate1["R_bc",]
  
  R2 = rate2["R",]
  d = rate2["b",]
  g = rate2["g",]
  R2_bc = rate2["R_bc",]
  var1_bc = rate1["var_Rbc",]
  var2_bc = rate2["var_Rbc",]
  r = lambda1 / lambda2
  
  ###1. hat_r
  hat_r = R1 / R2
  hat_r_p= performance(hat_r,r)
  per1 = c(P1$Np, P2$Np, rho1,rho2, r, "hat_r", hat_r_p, rep(NA, 2))
  write.table(t(per1), file =table_name, sep=",", na = "",
              row.names=FALSE,col.names=FALSE,append = TRUE)
 
  ###2. r_bc
  r_bc = hat_r - b/ R2 + R1 *d/ R2^2 - R1* g/ R2 ^3
  r_bc_p= performance(r_bc,r)
  var_rbc =(hat_f+ hat_r^2 *g - 4 * hat_r *b *d)/R2^2
  var_rbc_p = performance(var_rbc,r_bc_p[2])
  
  per2 = c(P1$Np, P2$Np,rho1,rho2, r,"r_bc", r_bc_p, var_rbc_p[-2])
  write.table(t(per2), file = table_name, sep=",", na = "",
              row.names=FALSE,col.names=FALSE,append = TRUE)
 
  ###3. r*
  r_star = R1_bc / R2_bc
  r_star_p= performance(r_star,r)
  var_rstar = (var1_bc*R2^2+var2_bc*R1^2)/R2^4
  var_rstar_p = performance(var_rstar,r_star_p[2])
  
  per3 = c(P1$Np, P2$Np,rho1,rho2, r,"r*", r_star_p, var_rstar_p[-2])
  write.table(t(per3), file = table_name, sep=",", na = "",
              row.names=FALSE,col.names=FALSE,append = TRUE)
  
  
  ##4. r*_bc
  rs_bc= r_star-R1_bc / R2_bc^3 *var2_bc
  rs_bc_p = performance(rs_bc, r)
  var_rsbc = (var1_bc*R2^2+var2_bc*R1^2)/R2^4 #same as var_rstar 
  var_rsbc_p = performance(var_rsbc,rs_bc_p[2])
  
  per4 = c(P1$Np, P2$Np,rho1,rho2, r,"r*_bc", rs_bc_p, var_rsbc_p[-2])
  write.table(t(per4), file = table_name, sep=",", na = "",
              row.names=FALSE,col.names=FALSE,append = TRUE)
  
}


#----------------------------------------------------------
# Simulation
#----------------------------------------------------------
set.seed(100)
#Generate population: Np, Na(j), n
#####Np: sample size; H: # of strata
population= function(Np,H,samp_ratio){
  result=list()
  result$Np= Np
  Nj = floor( Np * f)
  Nj[J] = Np - sum(Nj[1:(J-1)])
  result$Nj= Nj
  a = rep(1:J, Nj)
  result$a = a
  result$n = samp_ratio * Np
  #Stratified sampling setting
  prob = (1:H) / 10
  STR_samp = sample(1: H, Np, replace = TRUE, prob = prob)
  result$Ph = lapply(1:H, function(x)a[which(STR_samp == x)])
  Nst = as.vector(table(STR_samp))
  result$Nst= Nst
  result$nh = floor(samp_ratio * Nst)
  return(result)
}



######Fix Population
H = 4
Np = c(50000, 100000)
Pop1=population(Np[1],H,0.02)
Pop2=population(Np[2],H,0.01)
alpha = c(0.01, 0.05, 0.1)

#setwd("~/paper/ratio/ratio2")

# Table 1
r = c(1.0, 1.2, 1.5, 2.0)
lambda2 = 0.001
lambda1 = lambda2 * r
K=10000
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, SRS, Pop2, Pop2,0,0, "table1.csv"))
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, STR, Pop2, Pop2, 0,0, "table1.csv"))

# Table 2
r = c(0.25, 4)
lambda2 = 0.001
lambda1 = lambda2 * r
rho= c(0.05, 0.1, 0.15, 0.2)
for (i in 1:length(rho)){
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, ControlCV, Pop2, Pop2,rho[i],rho[i], "table2.csv"))
}

# Table 3
for (i in 1:length(rho)){
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, SRS, Pop2, Pop2,rho[i],rho[i], "table3.csv"))
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, STR, Pop2, Pop2, rho[i],rho[i], "table3.csv"))
}

#test effect of rho
rho= c(0.05, 0.2)
cvfunc = function(Pop, rho){
 Vj = rep(0,J)
 Ej = rep(0,J)
 CVj = rep(0,J)
 for(j in 1:J){
 mass= mass_func(Pop$Nj[j], rho)
 k=mass$lower:mass$upper
 Vj[j] = sum(k^2*mass$p)-(sum(k*mass$p))^2
 Ej[j] = sum(k*mass$p)
 CVj[j] = sqrt(Vj[j])/Ej[j]
 }
 return(c(min(CVj), max(CVj)))
}

#CVs of \hat{N}
cvfunc(Pop1, rho[1])
cvfunc(Pop1, rho[2])
cvfunc(Pop2, rho[1])
cvfunc(Pop2, rho[2])

# Table 4
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, ControlCV,
                                 Pop1, Pop1,rho[1],rho[2],"table4n.csv"))
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, ControlCV, 
                                 Pop1, Pop1, rho[2],rho[1],"table4n.csv"))
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, ControlCV, 
                                 Pop2, Pop2, rho[1],rho[2],"table4n.csv"))
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, ControlCV, 
                                 Pop2, Pop2, rho[2],rho[1],"table4n.csv"))

#Table 5
##fix P1, increase P2
tryN2 = c(1,5)*100000
Pop1 = population(100000,H,0.01)
for (i in 1: length(tryN2)){
Pop2 = population(tryN2[i],H,0.01)
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, SRS, Pop1,Pop2, 0,0,"table5n.csv"))
sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, STR, Pop1,Pop2, 0,0,"table5n.csv"))
}

# Table 6
##fix P2, decrease P1
tryN1 = c(100000,50000)
Pop2 = population(100000,H,0.01)
sample_ratio=c(0.01,0.02)
for (i in 1: length(tryN1)){
  Pop1 = population(tryN1[i],H,sample_ratio[i])
  sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, SRS, Pop1,Pop2, 0,0,"table6n.csv"))
  sapply(lambda1, function(x)ratio(x, lambda2, K, alpha, STR, Pop1,Pop2, 0,0,"table6n.csv"))
}





