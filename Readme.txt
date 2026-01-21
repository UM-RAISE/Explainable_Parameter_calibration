- Paper: Explainable Parameter Calibration via Importance-Driven Sequential Design with an Application to Building Energy Systems
- Journal: Annals of Applied Statistics
- Description: Reproduce the results in Table3 
- Instruction: 

(1) "Main" files
- File name: Table3_ExX_MethodName.R
- No. of files: 12 
- Setting path: each code will automatically set the working directory as the location where main R files are.
- Registering CPU cores: if no. of logical processors detected (i.e., detectCores()) is larger than 20 (= no. of experiments), each code will register 20 logical processors for experiments. Otherwise, each code will register detectCores()-1 logical processors for experiments.
- No. of experiments: 20 runs for each algorithm (Note: one experiment occupies one logical processor).
- Result: if you run each main file, the console in RStudio shows the result. Additionally, all the results will be saved in "Results" folder, so that "Results_Table3_method.R" produces the entire results in Table3.
- Author's computing environment: 64-bit Windows OS with the Intel Xeon CPU E5-2697 @ 2.60 GHz processor and 128GB RAM, R version: 4.2.2 (setting this R version is particularly important to reproduce the results).
- Specific file name and corresponding computing time:

File name | Mean time (std) with 1 exp/processor (unit:sec) | Total time if using only 1 processor sequentially (unit:sec) | Total time when 20 processors are used in parallel (unit:sec)
Table3_Ex1_Ex_Sliced_Calib.R   | 14151 (1917) | 283020 | 18369
Table3_Ex1_Rand_Sliced_Calib.R |  4775 (1046) |  95496 |  7147 
Table3_Ex1_Full_Calib.R        | 31906 (1645) | 638111 | 34509
Table3_Ex1_TR_Calib.R          | 27177 (1806) | 543540 | 30581 

Table3_Ex2_Ex_Sliced_Calib.R   | 16022 (1116) | 320444 | 18231
Table3_Ex2_Rand_Sliced_Calib.R |  7193 ( 690) | 143867 |  8363
Table3_Ex2_Full_Calib.R        | 31819 (1576) | 636372 | 34493
Table3_Ex2_TR_Calib.R          | 29605 (1163) | 592102 | 31491

Table3_Ex3_Ex_Sliced_Calib.R   | 25006 (1447) | 500130 | 28325
Table3_Ex3_Rand_Sliced_Calib.R |  8070 ( 524) | 161398 |  9080
Table3_Ex3_Full_Calib.R        | 30679 (1630) | 613571 | 32613
Table3_Ex3_TR_Calib.R          | 27755 (1806) | 555092 | 31275

- Specific outputs: 

train_MSE (mean) 
train_MSE (std)
test_MSE (mean)
test_MSE (std)
theta1 (mean)
theta1 (std)
theta2 (mean)
theta2 (std)
theta3 (mean)
theta3 (std)
time_onecore (mean)
time_onecore (std) 
total_time_onecore
total_time_multicore

(2) "Aggregate results" file (Optional)
- File name: Results_Table3_method.R
- Setting path: the code will automatically set the working directory as the location where the file is.
- This file will load all the results from 12 RData files from "Results" folder. 
- Specific outputs: 

train_MSE (mean) 
train_MSE (std)
test_MSE (mean)
test_MSE (std)
theta1 (mean)
theta1 (std)
theta2 (mean)
theta2 (std)
theta3 (mean)
theta3 (std)
time_onecore (mean)
time_onecore (std) 
total_time_onecore
total_time_multicore