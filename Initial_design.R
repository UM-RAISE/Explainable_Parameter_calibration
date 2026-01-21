
#############################################
# Title: Initial Design 
# Reps: 100 experiments 
#############################################

rm(list=ls(all=TRUE))

library(rstudioapi)
library(lhs)
library(MaxPro)

#############################################
# Problem Instance 
#############################################

# path
current_file = rstudioapi::getActiveDocumentContext()$path
work.path = dirname(current_file)
setwd(work.path)

source("Test_functions_MaxPro.R")

set.seed(2022)

# example
name = "hart3" # schwef (example1), ackley (example2), hart3 (example3)
func = hart3

reps = 100 # number of experiments
Nx_train = 200

if (name == "schwef") {
  d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
} else if (name == "ackley") {
  d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
} else if (name == "hart3") {
  d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
}

#############################################
# Algorithm
#############################################

initial_data = list()
Xlist_init = list() 
ylist_init = list()
errlist_init = list()

for (j in 1:reps){
  
  #==================================
  # observational error 
  obs.err = rnorm(Nx_train,mean=0,sd=sigma)
  #==================================
  
  InitialDesign<-MaxProLHD(n = ninit, p = D)$Design 
  DOX<-MaxPro(InitialDesign)
  X <- DOX$Design; colnames(X) = paste0("X",1:D)
  # X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
  y <- numeric(0)
  for (i in 1:ninit){
    temp = ff(obs.err, func(matrix(X[i,],nrow=1), d, D), Nx=Nx_train)
    y[i] = temp
  }
  Xlist_init[[j]] = X
  ylist_init[[j]] = y
  errlist_init[[j]] = obs.err
}

initial_data$Xlist_init = Xlist_init
initial_data$ylist_init = ylist_init
initial_data$errlist_init = errlist_init

# save
setwd(paste0(work.path,"/Data"))
save(initial_data, file=paste0("Initial_Data_",name,"_D_",D,".RData"))

