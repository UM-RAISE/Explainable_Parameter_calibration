
#############################################
# Title: Aggregate Results
#############################################

rm(list=ls(all=TRUE))

library(rstudioapi)

#############################################
# Problem Instance (Table 6)
#############################################

# path
current_file = rstudioapi::getActiveDocumentContext()$path
work.path = dirname(current_file)
setwd(paste0(work.path,"/Results"))

# name 
name.vec = c("schwef","ackley","hart3")  
sigma0.vec = c(0.4,0.5,0.6) 

# algorithmic parameters
reps_test = 20 # number of experiments
Nx_train = 200
Nx_test = 200

trial_dim_n = 10
# sigma0_n = 0.5
# trial_tau0_n = 1/(sigma0_n^2)
trial_p_n = 0.05

table = matrix(NA,length(name.vec)*length(sigma0.vec),14)

ii = 1
for (name_n in name.vec) {
  
  # example
  if (name_n == "schwef") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
  } else if (name_n == "ackley") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
  } else if (name_n == "hart3") {
    d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
  }
  
  for (sigma0_n in sigma0.vec) {
    
    # load 
    load(file=paste0(name_n,"_A_Ex_Sliced_Calib_D_",D,"_d_",trial_dim_n,"_sigma0_",gsub("\\.", "", sigma0_n),"_p_",gsub("\\.", "", trial_p_n),".RData"))  
    
    ybest_vec = numeric(0)
    ybest_test_vec = numeric(0)
    Xbest_vec = matrix(NA,reps_test,D)
    time_record_vec = numeric(0)
    
    for (r in c(1:reps_test)) {
      ybest_vec[r] = rslt[[r]]$ybest
      ybest_test_vec[r] = rslt[[r]]$ybest_test
      Xbest_vec[r,] = rslt[[r]]$Xbest
      time_record_vec[r] = rslt[[r]]$time_record
    }
    
    table[ii,] = c(mean(ybest_vec),sd(ybest_vec),
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
    rownames(table) = c("Ex1_sigma0_04","Ex1_sigma0_05","Ex1_sigma0_06",
                        "Ex2_sigma0_04","Ex2_sigma0_05","Ex2_sigma0_06",
                        "Ex3_sigma0_04","Ex3_sigma0_05","Ex3_sigma0_06")
    ii = ii + 1  
  }
}

print(table)

#############################################
# Scaling
#############################################

table[1:3,colnames(table) == c("theta1 (mean)")] = table[1:3,colnames(table) == c("theta1 (mean)")]*(500-(-500))-500
table[1:3,colnames(table) == c("theta2 (mean)")] = table[1:3,colnames(table) == c("theta2 (mean)")]*(500-(-500))-500
table[1:3,colnames(table) == c("theta3 (mean)")] = NA

table[1:3,colnames(table) == c("theta1 (std)")] = table[1:3,colnames(table) == c("theta1 (std)")]*(500-(-500))
table[1:3,colnames(table) == c("theta2 (std)")] = table[1:3,colnames(table) == c("theta2 (std)")]*(500-(-500))
table[1:3,colnames(table) == c("theta3 (std)")] = NA

table[4:6,colnames(table) == c("theta1 (mean)")] = table[4:6,colnames(table) == c("theta1 (mean)")]*(10-(-5))-5
table[4:6,colnames(table) == c("theta2 (mean)")] = table[4:6,colnames(table) == c("theta2 (mean)")]*(10-(-5))-5
table[4:6,colnames(table) == c("theta3 (mean)")] = NA

table[4:6,colnames(table) == c("theta1 (std)")] = table[4:6,colnames(table) == c("theta1 (std)")]*(10-(-5))
table[4:6,colnames(table) == c("theta2 (std)")] = table[4:6,colnames(table) == c("theta2 (std)")]*(10-(-5))
table[4:6,colnames(table) == c("theta3 (std)")] = NA

print(table)
write.table(table, file=paste0(work.path,"/Results/Table_sigma0.txt"), row.names=T, col.names=T, sep="|")

