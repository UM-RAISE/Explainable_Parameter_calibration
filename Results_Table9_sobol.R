
#############################################
# Title: Aggregate Results
#############################################

rm(list=ls(all=TRUE))

library(rstudioapi)
library(lhs)
library(sensitivity)
library(MaxPro)

#############################################
# Problem Instance (Table 9) - Proposed
#############################################

# path
current_file = rstudioapi::getActiveDocumentContext()$path
work.path = dirname(current_file)
setwd(paste0(work.path,"/Results"))

# name 
name_n = "hart3" # schwef,ackley,hart3
method_n = "A_Ex_Sliced_Calib"

# algorithmic parameters
reps_test = 20 # number of experiments
Nx_train = 200
Nx_test = 200

trial_dim_n = 10
sigma0_n = 0.5
trial_tau0_n = 1/(sigma0_n^2)
trial_p_n = 0.05

# example
if (name_n == "schwef") {
  d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
} else if (name_n == "ackley") {
  d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
} else if (name_n == "hart3") {
  d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
}
  
# load 
load(file=paste0(name_n,"_",method_n,"_D_",D,"_d_",trial_dim_n,"_sigma0_",gsub("\\.", "", sigma0_n),"_p_",gsub("\\.", "", trial_p_n),".RData"))  

table = matrix(NA,D,4)

mu_mat = matrix(NA,D,reps_test)

for (r in c(1:reps_test)) {
  temp = rslt[[r]]$mu
  mu_mat[,r] = temp[end,]
}

table[,1] = apply(mu_mat,1,mean,na.rm=T)
table[,2] = apply(mu_mat,1,sd,na.rm=T)

#############################################
# Problem Instance (Table 9) - Sobol
#############################################

setwd(work.path)
source("Test_functions_sobol.R")
func = hart3

errlist = list()
for (r in c(1:reps_test)) {
  errlist[[r]] = rslt[[r]]$obs.err
}

n <- 10000 # 10^5 is needed

test_function = function(X_arg,obs.err,Nx) {
  
  temp = ff(obs.err, func(X_arg,d=d,D=D),Nx)
  return(temp)
}

# sensitivity analysis
result_sobol_main = matrix(NA,D,reps_test)
result_sobol_total= matrix(NA,D,reps_test)

for (r in c(1:reps_test)) {
  
  X1 <- data.frame(matrix(runif(D*n), nrow = n))
  X2 <- data.frame(matrix(runif(D*n), nrow = n))
  
  # sobol indices
  sobol_result <- sobol2002(model=test_function, X1, X2, obs.err = errlist[[r]], Nx=Nx_train)
  sobol_main = sobol_result$S$original
  sobol_total = sobol_result$T$original
  
  result_sobol_main[,r] = sobol_main
  result_sobol_total[,r] = sobol_total
}

table[,3] = apply(result_sobol_total,1,mean,na.rm=T)
table[,4] = apply(result_sobol_total,1,sd,na.rm=T)

colnames(table) = c("Proposed (mean)","Proposed (std)","Sobol (mean)","Sobol (std)")
rownames(table) = paste("theta",1:D)

print(table)

setwd(paste0(work.path,"/Results"))
write.csv(table, file=paste0("sobol_",name_n,".txt"))

