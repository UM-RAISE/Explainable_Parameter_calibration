
#############################################
# Title: Simulation Study
# Example 2
# Full_Calib
#############################################

rm(list=ls(all=TRUE))

library(rstudioapi)
library(foreach)
library(doParallel)
library(lhs)
library(MaxPro)

#############################################
# Problem Instance 
#############################################

# path
current_file = rstudioapi::getActiveDocumentContext()$path
work.path = dirname(current_file)
setwd(work.path)

source("parallel.R")
source("Test_functions_MaxPro.R")
# source("Test_functions_maximinLHS.R")

# example
name = "ackley" # schwef (example1), ackley (example2), hart3 (example3)
func = ackley

# method: A_Ex_Sliced_Calib, A_Rand_Sliced_Calib, A_Full_Calib, A_TR_Calib
method = "A_Full_Calib"

# algorithmic parameters
reps = 20 # number of experiments
Nx_train = 200
Nx_test = 200

trial_dim = 10
sigma0 = 0.5
trial_tau0 = 1/(sigma0^2)
trial_p = 0.05

if (name == "schwef") {
  
  d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
  load(file=paste0(work.path,"/Data/Initial_Data_schwef_D_30.RData"))
  
} else if (name == "ackley") {
  
  d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
  load(file=paste0(work.path,"/Data/Initial_Data_ackley_D_30.RData"))
  
} else if (name == "hart3") {
  
  d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
  load(file=paste0(work.path,"/Data/Initial_Data_hart3_D_30.RData"))
}

filename = paste0(method,".R")
source(filename)

#############################################
# Algorithm
#############################################

cores_temp = detectCores() - 1
if (cores_temp >= reps) {
  cores = reps
} else {
  cores = cores_temp 
}

cl = makeCluster(cores)
registerDoParallel(cl)

entire_tic <- proc.time()[3]

rslt = foreach(r=1:reps,.multicombine=TRUE) %dopar% {
  
  result = multiResultClass()
  
  #==================================
  # package
  #==================================
  library(lhs)
  library(DiceKriging)
  library(DiceOptim)
  library(sensitivity)
  library(MaxPro)
  library(mvtnorm)
  library(truncnorm)

  #==================================
  # main
  #==================================  
  set.seed(r)
  
  # data
  # temp = init_XY(name,Ninit=ninit,D=D,d=d,Nx=Nx_train,seed=r)
  # X = temp$Xlist
  # y = temp$ylist
  # obs.err = temp$errlist 
  
  X = initial_data$Xlist_init[[r]]
  y = initial_data$ylist_init[[r]]
  obs.err = initial_data$errlist_init[[r]]
  
  # main
  tic <- proc.time()[3]
  os <- optim.EI(f=func, ninit, end, obs.err, trial_dim, trial_tau0, trial_p, X, y, D, d, Nx=Nx_train, name, seed=r)
  toc <- proc.time()[3]
  time_record <- toc - tic
  
  #==================================
  # training
  #==================================
  yrslt = os$y
  yprog = bov(os$y)
  Xrslt <- os$X
  tmp = as.vector(yprog)
  idx = numeric(0)
  for (i in 1:end) {
    idx[i] = which(yrslt == tmp[i])
  }
  Xprog = Xrslt[idx,]
  m = which.min(yrslt)
  ybest = yrslt[m]
  Xbest = Xrslt[m,]
  
  #==================================
  # test
  #==================================
  set.seed(r+100)
  obs.err_test = rnorm(Nx_test,mean=0,sd=sigma)
  ybest_test = ff(obs.err_test, func(matrix(Xbest,nrow=1), d, D), Nx=Nx_test)
  
  #==================================
  # record
  #==================================
  result$X = X
  result$y = y
  result$obs.err = obs.err
  
  result$Xrslt = Xrslt
  result$yrslt = yrslt
  result$Xprog = Xprog
  result$yprog = yprog
  result$Xbest = Xbest
  result$ybest = ybest
  
  result$ybest_test = ybest_test
  result$obs.err_test = obs.err_test
  
  result$time_record = time_record
  result$major_rslt = os$major_matrix
  result$minor_rslt = os$minor_matrix
  result$Tk_rslt = os$Tk
  result$weight_rslt = os$weight
  result$nu_rslt = os$nu
  result$tau_rslt = os$tau
  result$mu_rslt = os$mu
  
  setwd(paste0(work.path,"/Results/Raw"))
  save.image(file=paste0(name,"_",method,"_D_",D,"_d_",trial_dim,"_sigma0_",gsub("\\.", "", sigma0),"_p_",gsub("\\.", "", trial_p),"_",r,".RData"))
  
  return(result)
}
stopCluster(cl)

# elapse time
entire_toc <- proc.time()[3]
entire_time = entire_toc - entire_tic
print(entire_time)

rslt$entire_time = entire_time

# save
setwd(paste0(work.path,"/Results"))
save(rslt, file=paste0(name,"_",method,"_D_",D,"_d_",trial_dim,"_sigma0_",gsub("\\.", "", sigma0),"_p_",gsub("\\.", "", trial_p),".RData"))

#############################################
# Result
#############################################

table = matrix(NA,1,14)

ybest_vec = numeric(0)
ybest_test_vec = numeric(0)
Xbest_vec = matrix(NA,reps,D)
time_record_vec = numeric(0)

for (r in c(1:reps)) {
  ybest_vec[r] = rslt[[r]]$ybest
  ybest_test_vec[r] = rslt[[r]]$ybest_test
  Xbest_vec[r,] = rslt[[r]]$Xbest
  time_record_vec[r] = rslt[[r]]$time_record
}

table[1,] = c(mean(ybest_vec),sd(ybest_vec),
              mean(ybest_test_vec),sd(ybest_test_vec),
              mean(Xbest_vec[,1]),sd(Xbest_vec[,1]),
              mean(Xbest_vec[,2]),sd(Xbest_vec[,2]),
              mean(Xbest_vec[,3]),sd(Xbest_vec[,3]),
              mean(time_record_vec),sd(time_record_vec),sum(time_record_vec),rslt$entire_time)  

colnames(table) = c("train_MSE (mean)","train_MSE (std)",
                    "test_MSE (mean)","test_MSE (std)",
                    "theta1 (mean)","theta1 (std)",
                    "theta2 (mean)","theta2 (std)",
                    "theta3 (mean)","theta3 (std)",
                    "time_onecore (mean)", "time_onecore (std)","total_time_onecore",
                    "total_time_multicore")
rownames(table) = c("result")

table[1,colnames(table) == c("theta1 (mean)")] = table[1,colnames(table) == c("theta1 (mean)")]*(10-(-5))-5
table[1,colnames(table) == c("theta2 (mean)")] = table[1,colnames(table) == c("theta2 (mean)")]*(10-(-5))-5
table[1,colnames(table) == c("theta3 (mean)")] = NA

table[1,colnames(table) == c("theta1 (std)")] = table[1,colnames(table) == c("theta1 (std)")]*(10-(-5))
table[1,colnames(table) == c("theta2 (std)")] = table[1,colnames(table) == c("theta2 (std)")]*(10-(-5))
table[1,colnames(table) == c("theta3 (std)")] = NA

print(table)

