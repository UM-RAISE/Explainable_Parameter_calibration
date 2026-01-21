
#############################################
# Title: Figures
#############################################

rm(list=ls(all=TRUE))

library(rstudioapi)
library(scales)
library(ggplot2)
library(rstudioapi)

## path
current_file = rstudioapi::getActiveDocumentContext()$path
work.path = dirname(current_file)
setwd(paste0(work.path,"/Results"))

name.vec = c("schwef","ackley","hart3")  
method.vec = c("A_Ex_Sliced_Calib","A_Rand_Sliced_Calib","A_Full_Calib","A_TR_Calib") 
dim.vec = c(10,20,25,30) 
p.vec = c(0,0.05,0.25,0.5) 
sigma0.vec = c(0.4,0.5,0.6) 

reps_test = 20 # number of experiments
trial_dim_n = 10
sigma0_n = 0.5
trial_tau0_n = 1/(sigma0_n^2)
trial_p_n = 0.05

################################################################################
# 2-1.Comparison with other BO alternatives
# --- MSE and calibrated values
# --- MSE over rounds
################################################################################

for (name_n in name.vec) {
  
  # example
  if (name_n == "schwef") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
  } else if (name_n == "ackley") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
  } else if (name_n == "hart3") {
    d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
  }
  
  dat = c()
  
  for (method_n in method.vec) { 
    
    load(file=paste0(name_n,"_",method_n,"_D_",D,"_d_",trial_dim_n,"_sigma0_",gsub("\\.", "", sigma0_n),"_p_",gsub("\\.", "", trial_p_n),".RData"))  
    xx = ninit:end
    
    yprog_mat = matrix(NA,reps_test,end-ninit+1)
    for (r in c(1:reps_test)) {
      yprog_mat[r,] = rslt[[r]]$yprog[xx]
    }
    
    yprog_mat_mean = apply(yprog_mat,2,mean,na.rm=TRUE)
    yprog_mat_sd = apply(yprog_mat,2,sd,na.rm=TRUE)
    
    dat_temp = rbind(yprog_mat_mean, yprog_mat_sd)
    dat = rbind(dat, dat_temp)
  }
  
  dat = cbind(ninit:end,t(dat))
  colnames(dat) = c("iter","sliced_mean","sliced_sd",
                    "rand_sliced_mean","rand_sliced_sd",
                    "full_mean","full_sd",
                    "tr_mean","tr_sd") 
  dat = data.frame(dat)
  assign( paste0("data_",name_n), dat)
}

head(data_schwef)
head(data_ackley)
head(data_hart3)

################################################################################
# 2-2.Quantifying Parameter Importance
# --- Weights
# --- Probability Mass
################################################################################

method_n = "A_Ex_Sliced_Calib"

for (name_n in name.vec) {
  
  # example
  if (name_n == "schwef") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
  } else if (name_n == "ackley") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
  } else if (name_n == "hart3") {
    d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
  }
  
  load(file=paste0(name_n,"_",method_n,"_D_",D,"_d_",trial_dim_n,"_sigma0_",gsub("\\.", "", sigma0_n),"_p_",gsub("\\.", "", trial_p_n),".RData")) 
  xx = (ninit+1):end
  
  weight_temp = matrix(NA,nrow=reps_test,ncol=end-ninit)
  
  for (ii in c(1:D)) {
    
    ridx = 1
    for (jj in c(1:reps_test)) {
      weight_temp[ridx,] = rslt[[jj]]$mu[,ii][(ninit+1):end]
      ridx = ridx + 1
    }
    assign(paste0("weight",ii), weight_temp)
  }
  
  dat.weight = c()

  for (ii in c(1:D)) {
    
    weight.temp = eval(parse(text=paste0("weight",ii)))
    weight.temp.mean = apply(weight.temp,2,mean,na.rm=TRUE)
    weight.temp.sd = apply(weight.temp,2,sd,na.rm=TRUE)

    dat.weight = cbind(dat.weight, weight.temp.mean)
  }
  
  dat.weight = cbind((ninit+1):end,dat.weight)

  colnames(dat.weight) = c("iter",paste0("theta",1:D))
  dat.weight = as.data.frame(dat.weight)

  assign( paste0("data_weight_",name_n), dat.weight)
}

head(data_weight_schwef)
head(data_weight_ackley)
head(data_weight_hart3)

################################################################################
# 2-3.Sensitivity Analysis in Algorithmic Parameters
# (1) Effect of varying d
################################################################################

trial_dim_n = 10
sigma0_n = 0.5
trial_p_n = 0.05

for (name_n in name.vec) {
  
  # example
  if (name_n == "schwef") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
  } else if (name_n == "ackley") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
  } else if (name_n == "hart3") {
    d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
  }
  
  dat = c()
  
  for (trial_dim_n in dim.vec) { 
    
    load(file=paste0(name_n,"_",method_n,"_D_",D,"_d_",trial_dim_n,"_sigma0_",gsub("\\.", "", sigma0_n),"_p_",gsub("\\.", "", trial_p_n),".RData")) 
    xx = ninit:end
    
    yprog_mat = matrix(NA,reps_test,end-ninit+1)
    for (r in c(1:reps_test)) {
      yprog_mat[r,] = rslt[[r]]$yprog[xx]
    }
    
    yprog_mat_mean = apply(yprog_mat,2,mean,na.rm=TRUE)
    yprog_mat_sd = apply(yprog_mat,2,sd,na.rm=TRUE)

    dat_temp = rbind(yprog_mat_mean, yprog_mat_sd)
    dat = rbind(dat, dat_temp)
  }
  
  dat = cbind(ninit:end,t(dat))
  colnames(dat) = c("iter","d10_mean","d10_sd",
                    "d20_mean","d20_sd",
                    "d25_mean","d25_sd",
                    "d30_mean","d30_sd") 
  dat = data.frame(dat)
  assign( paste0("data_d_",name_n), dat)
}

head(data_d_schwef)
head(data_d_ackley)
head(data_d_hart3)

################################################################################
# 2-3.Sensitivity Analysis in Algorithmic Parameters
# (2) Effect of varying p
################################################################################

trial_dim_n = 10
sigma0_n = 0.5
trial_p_n = 0.05

for (name_n in name.vec) {
  
  # example
  if (name_n == "schwef") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
  } else if (name_n == "ackley") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
  } else if (name_n == "hart3") {
    d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
  }
  
  dat = c()
  
  for (trial_p_n in p.vec) { 
    
    load(file=paste0(name_n,"_",method_n,"_D_",D,"_d_",trial_dim_n,"_sigma0_",gsub("\\.", "", sigma0_n),"_p_",gsub("\\.", "", trial_p_n),".RData")) 
    xx = ninit:end
    
    yprog_mat = matrix(NA,reps_test,end-ninit+1)
    for (r in c(1:reps_test)) {
      yprog_mat[r,] = rslt[[r]]$yprog[xx]
    }
    
    yprog_mat_mean = apply(yprog_mat,2,mean,na.rm=TRUE)
    yprog_mat_sd = apply(yprog_mat,2,sd,na.rm=TRUE)
    
    dat_temp = rbind(yprog_mat_mean, yprog_mat_sd)
    dat = rbind(dat, dat_temp)
  }
  
  dat = cbind(ninit:end,t(dat))
  colnames(dat) = c("iter","p0_mean","p0_sd",
                    "p005_mean","p005_sd",
                    "p025_mean","p025_sd",
                    "p05_mean","p05_sd") 
  dat = data.frame(dat)
  assign( paste0("data_p_",name_n), dat)
}

head(data_p_schwef)
head(data_p_ackley)
head(data_p_hart3)

################################################################################
# 2-3.Sensitivity Analysis in Algorithmic Parameters
# (3) Effect of varying sigma0
################################################################################

trial_dim_n = 10
sigma0_n = 0.5
trial_p_n = 0.05

for (name_n in name.vec) {
  
  # example
  if (name_n == "schwef") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,300); sigma=0.1
  } else if (name_n == "ackley") {
    d=2; D=30; ninit=5*D; end=500; ylimit=c(0,5); sigma=0.1
  } else if (name_n == "hart3") {
    d=3; D=30; ninit=5*D; end=500; ylimit=c(0,0.4); sigma=0.01
  }
  
  dat = c()
  
  for (sigma0_n in sigma0.vec) { 
    
    load(file=paste0(name_n,"_",method_n,"_D_",D,"_d_",trial_dim_n,"_sigma0_",gsub("\\.", "", sigma0_n),"_p_",gsub("\\.", "", trial_p_n),".RData")) 
    xx = ninit:end
    
    yprog_mat = matrix(NA,reps_test,end-ninit+1)
    for (r in c(1:reps_test)) {
      yprog_mat[r,] = rslt[[r]]$yprog[xx]
    }
    
    yprog_mat_mean = apply(yprog_mat,2,mean,na.rm=TRUE)
    yprog_mat_sd = apply(yprog_mat,2,sd,na.rm=TRUE)
    
    dat_temp = rbind(yprog_mat_mean, yprog_mat_sd)
    dat = rbind(dat, dat_temp)
  }
  
  dat = cbind(ninit:end,t(dat))
  colnames(dat) = c("iter","sigma04_mean","sigma04_sd",
                    "sigma05_mean","sigma05_sd",
                    "sigma06_mean","sigma06_sd") 
  dat = data.frame(dat)
  assign( paste0("data_sigma_",name_n), dat)
}

head(data_sigma_schwef)
head(data_sigma_ackley)
head(data_sigma_hart3)

################################################################################
# ggplot - 2-1.Comparison with other BO alternatives
################################################################################

sqrtn = sqrt(reps_test)
rib1 <- aes(x = iter, y = sliced_mean,
            ymax = sliced_mean + sliced_sd/sqrtn, 
            ymin = sliced_mean - sliced_sd/sqrtn)
rib2 <- aes(x = iter, y = rand_sliced_mean,
            ymax = rand_sliced_mean + rand_sliced_sd/sqrtn, 
            ymin = rand_sliced_mean - rand_sliced_sd/sqrtn)
rib3 <- aes(x = iter, y = full_mean,
            ymax = full_mean + full_sd/sqrtn, 
            ymin = full_mean - full_sd/sqrtn)
rib4 <- aes(x = iter, y = tr_mean,
            ymax = tr_mean + tr_sd/sqrtn, 
            ymin = tr_mean - tr_sd/sqrtn)

# data_schwef, data_ackley, data_hart3
# legend.position = c(.77, .4), c(.78, .86)
legend.list = c("Ex-Sliced-Calib","Rand-Sliced-Calib","Full-Calib","TR-Calib")
legend.list2 <- c("Ex-Sliced-Calib"="#F8766D","Rand-Sliced-Calib"="#7CAE00","Full-Calib"="#00BFC4","TR-Calib"="#C77CFF")

ggplot(data = data_hart3) + 
  geom_line(size = 4, aes(x = iter, y = tr_mean, color="TR-Calib")) +
  geom_ribbon(rib4, alpha = 0.3, fill="#C77CFF") +
  geom_line(size = 4, aes(x = iter, y = full_mean, color="Full-Calib")) + 
  geom_ribbon(rib3, alpha = 0.3, fill="#00BFC4") +
  geom_line(size = 4, aes(x = iter, y = rand_sliced_mean, color="Rand-Sliced-Calib")) + 
  geom_ribbon(rib2, alpha = 0.3, fill="#7CAE00") + 
  geom_line(size = 4, aes(x = iter, y = sliced_mean, color="Ex-Sliced-Calib")) + 
  geom_ribbon(rib1, alpha = 0.3, fill="#F8766D") +
  theme_minimal() +
  labs(x = "Number of evaluations", y = "Best MSE") + # color="legend.list2"
  theme(axis.text = element_text(size = 55),
        axis.title.x = element_text(size = 75),
        axis.title.y = element_text(size = 75),
        axis.line = element_line(color="gray20", size=1),
        legend.position = c(.7, .85), # c(.705, .85) c(.72, .85) c(.7, .85)
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.text = element_text(size=50),
        legend.background = element_rect(color = NA, fill = "transparent")) +
  geom_vline(xintercept = 150, linetype="dashed", size=2, color="gray40") +
  scale_x_continuous(breaks=round(seq(150,500,by=50),1)) +
  scale_color_manual(breaks=legend.list, values=legend.list2) 

################################################################################
# ggplot - 2-2.Quantifying Parameter Importance
# - Probability Mass
################################################################################

# data_weight_schwef, data_weight_ackley
# legend.position = c(.22, .85), c(.8, .85)
legend.list = c("theta1","theta2","theta3-30") # c("theta 1","theta 2","theta 3","theta 4-30")
legend.list2 <- c("theta1"="#F8766D","theta2"="#7CAE00","theta3-30"="gray60")

ggplot(data = data_weight_ackley) + 
  geom_line(size = 2, aes(x = iter, y = theta3, color="theta3-30")) + 
  geom_line(size = 2, aes(x = iter, y = theta4, color="theta3-30")) +   
  geom_line(size = 2, aes(x = iter, y = theta5, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta6, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta7, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta8, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta9, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta10, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta11, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta12, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta13, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta14, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta15, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta16, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta17, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta18, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta19, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta20, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta21, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta23, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta24, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta25, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta26, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta27, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta28, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta29, color="theta3-30")) +
  geom_line(size = 2, aes(x = iter, y = theta30, color="theta3-30")) +
  geom_line(size = 4, aes(x = iter, y = theta1, color="theta1")) + 
  geom_line(size = 4, aes(x = iter, y = theta2, color="theta2")) + 
  theme_minimal() +
  # ylim(c(0.032,0.05)) +
  labs(x = "Number of evaluations", y = "Posterior mean") + 
  theme(axis.text = element_text(size = 55),
        axis.title.x = element_text(size = 75),
        axis.title.y = element_text(size = 75),
        axis.line = element_line(color="gray20", size=1),
        legend.position = c(.748, .25), # c(.748, .3) c(.748, .25)
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.text = element_text(size=60),
        legend.background = element_rect(color = NA, fill = "transparent")) +
  geom_vline(xintercept = 151, linetype="dashed", size=2, color="gray60") +
  scale_x_continuous(breaks=round(seq(150,500,by=50),1)) +
  scale_color_manual(breaks=legend.list, values=legend.list2) 

# data_weight_hart3
# legend.position = c(.8, .85)
legend.list = c("theta1","theta2","theta3","theta4-30") # c("theta 1","theta 2","theta 3","theta 4-30")
legend.list2 <- c("theta1"="#F8766D","theta2"="#7CAE00","theta3"="#00BFC4","theta4-30"="gray60")

ggplot(data = data_weight_hart3) + 
  geom_line(size = 2, aes(x = iter, y = theta4, color="theta4-30")) +   
  geom_line(size = 2, aes(x = iter, y = theta5, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta6, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta7, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta8, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta9, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta10, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta11, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta12, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta13, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta14, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta15, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta16, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta17, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta18, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta19, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta20, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta21, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta23, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta24, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta25, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta26, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta27, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta28, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta29, color="theta4-30")) +
  geom_line(size = 2, aes(x = iter, y = theta30, color="theta4-30")) +
  geom_line(size = 4, aes(x = iter, y = theta1, color="theta1")) + 
  geom_line(size = 4, aes(x = iter, y = theta2, color="theta2")) + 
  geom_line(size = 4, aes(x = iter, y = theta3, color="theta3")) + 
  theme_minimal() +
  # ylim(c(0.03,0.07)) +
  labs(x = "Number of evaluations", y = "Posterior mean") + 
  theme(axis.text = element_text(size = 55),
        axis.title.x = element_text(size = 75),
        axis.title.y = element_text(size = 75),
        axis.line = element_line(color="gray20", size=1),
        legend.position = c(.748, .36), # c(.8, .85) c(.748, .25) c(.748, .436)
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.text = element_text(size=60),
        legend.background = element_rect(color = NA, fill = "transparent")) +
  geom_vline(xintercept = 151, linetype="dashed", size=2, color="gray60") +
  scale_x_continuous(breaks=round(seq(150,500,by=50),1)) +
  scale_color_manual(breaks=legend.list, values=legend.list2) 

################################################################################
# ggplot - 2-3.Sensitivity Analysis in Algorithmic Parameters
# (1) Effect of varying d
################################################################################

sqrtn = sqrt(reps_test)
rib1 <- aes(x = iter, y = d10_mean,
            ymax = d10_mean + d10_sd/sqrtn, 
            ymin = d10_mean - d10_sd/sqrtn)
rib2 <- aes(x = iter, y = d20_mean,
            ymax = d20_mean + d20_sd/sqrtn, 
            ymin = d20_mean - d20_sd/sqrtn)
rib3 <- aes(x = iter, y = d25_mean,
            ymax = d25_mean + d25_sd/sqrtn, 
            ymin = d25_mean - d25_sd/sqrtn)
rib4 <- aes(x = iter, y = d30_mean,
            ymax = d30_mean + d30_sd/sqrtn, 
            ymin = d30_mean - d30_sd/sqrtn)

# data_d_schwef, data_d_ackley, data_d_hart3
# legend.position = c(.77, .35), c(.77, .8), c(.77, .8) 
legend.list = c("d=10","d=20","d=25","d=30")
legend.list2 <- c("d=10"="#F8766D","d=20"="#7CAE00","d=25"="#00BFC4","d=30"="#C77CFF")

ggplot(data = data_d_hart3) + 
  geom_line(size = 4, aes(x = iter, y = d30_mean, color="d=30")) + 
  geom_ribbon(rib4, alpha = 0.3, fill="#C77CFF") +
  geom_line(size = 4, aes(x = iter, y = d25_mean, color="d=25")) + 
  geom_ribbon(rib3, alpha = 0.3, fill="#00BFC4") +
  geom_line(size = 4, aes(x = iter, y = d20_mean, color="d=20")) + 
  geom_ribbon(rib2, alpha = 0.3, fill="#7CAE00") + 
  geom_line(size = 4, aes(x = iter, y = d10_mean, color="d=10")) + 
  geom_ribbon(rib1, alpha = 0.3, fill="#F8766D") +
  theme_minimal() +
  labs(x = "Number of evaluations", y = "Best MSE") + # color="legend.list2"
  theme(axis.text = element_text(size = 55),
        axis.title.x = element_text(size = 75),
        axis.title.y = element_text(size = 75),
        axis.line = element_line(color="gray20", size=1),
        legend.position = c(.85, .85),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.text = element_text(size=50),
        legend.background = element_rect(color = NA, fill = "transparent")) +
  geom_vline(xintercept = 150, linetype="dashed", size=2, color="gray40") +
  scale_x_continuous(breaks=round(seq(150,500,by=50),1)) +
  scale_color_manual(breaks=legend.list, values=legend.list2) 

################################################################################
# ggplot - 2-3.Sensitivity Analysis in Algorithmic Parameters
# (3) Effect of varying p
################################################################################

sqrtn = sqrt(reps_test)
rib1 <- aes(x = iter, y = p005_mean,
            ymax = p005_mean + p005_sd/sqrtn, 
            ymin = p005_mean - p005_sd/sqrtn)
rib2 <- aes(x = iter, y = p05_mean,
            ymax = p05_mean + p05_sd/sqrtn, 
            ymin = p05_mean - p05_sd/sqrtn)
rib3 <- aes(x = iter, y = p025_mean,
            ymax = p025_mean + p025_sd/sqrtn, 
            ymin = p025_mean - p025_sd/sqrtn)
rib4 <- aes(x = iter, y = p0_mean,
            ymax = p0_mean + p0_sd/sqrtn, 
            ymin = p0_mean - p0_sd/sqrtn)

# data_p_schwef, data_p_ackley, data_p_hart3
# legend.position = c(.77, .35), c(.77, .8), c(.77, .8) 
legend.list = c("p=0.5","p=0.25","p=0.05","p=0")
legend.list2 <- c("p=0.05"="#F8766D","p=0"="#7CAE00","p=0.25"="#00BFC4","p=0.5"="#C77CFF")

ggplot(data = data_p_hart3) + 
  geom_line(size = 4, aes(x = iter, y = p05_mean, color="p=0.5")) + 
  geom_ribbon(rib2, alpha = 0.3, fill="#C77CFF") + 
  geom_line(size = 4, aes(x = iter, y = p025_mean, color="p=0.25")) + 
  geom_ribbon(rib3, alpha = 0.3, fill="#00BFC4") +
  geom_line(size = 4, aes(x = iter, y = p0_mean, color="p=0")) + 
  geom_ribbon(rib4, alpha = 0.3, fill="#7CAE00") +
  geom_line(size = 4, aes(x = iter, y = p005_mean, color="p=0.05")) + 
  geom_ribbon(rib1, alpha = 0.3, fill="#F8766D") +
  theme_minimal() +
  labs(x = "Number of evaluations", y = "Best MSE") + # color="legend.list2"
  theme(axis.text = element_text(size = 55),
        axis.title.x = element_text(size = 75),
        axis.title.y = element_text(size = 75),
        axis.line = element_line(color="gray20", size=1),
        legend.position = c(.83, .85),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.text = element_text(size=50),
        legend.background = element_rect(color = NA, fill = "transparent")) +
  geom_vline(xintercept = 150, linetype="dashed", size=2, color="gray40") +
  scale_x_continuous(breaks=round(seq(150,500,by=50),1)) +
  scale_color_manual(breaks=legend.list, values=legend.list2) 

################################################################################
# ggplot - 2-3.Sensitivity Analysis in Algorithmic Parameters
# (3) Effect of varying sigma0
################################################################################

sqrtn = sqrt(reps_test)
rib1 <- aes(x = iter, y = sigma04_mean,
            ymax = sigma04_mean + sigma04_sd/sqrtn, 
            ymin = sigma04_mean - sigma04_sd/sqrtn)
rib2 <- aes(x = iter, y = sigma05_mean,
            ymax = sigma05_mean + sigma05_sd/sqrtn, 
            ymin = sigma05_mean - sigma05_sd/sqrtn)
rib3 <- aes(x = iter, y = sigma06_mean,
            ymax = sigma06_mean + sigma06_sd/sqrtn, 
            ymin = sigma06_mean - sigma06_sd/sqrtn)

# data_sigma_schwef, data_sigma_ackley, data_sigma_hart3
# legend.position = c(.77, .35), c(.77, .8), c(.77, .8) 
legend.list = c("sigma0=0.6","sigma0=0.5","sigma0=0.4")
legend.list2 <- c("sigma0=0.5"="#F8766D","sigma0=0.4"="#7CAE00","sigma0=0.6"="#00BFC4")

ggplot(data = data_sigma_hart3) + 
  # geom_line(size = 4, aes(x = iter, y = gamma1_mean, color="gamma=1")) + 
  # geom_ribbon(rib4, alpha = 0.3, fill="#C77CFF") +
  geom_line(size = 4, aes(x = iter, y = sigma06_mean, color="sigma0=0.6")) + 
  geom_ribbon(rib3, alpha = 0.3, fill="#00BFC4") +
  geom_line(size = 4, aes(x = iter, y = sigma04_mean, color="sigma0=0.4")) + 
  geom_ribbon(rib1, alpha = 0.3, fill="#7CAE00") + 
  geom_line(size = 4, aes(x = iter, y = sigma05_mean, color="sigma0=0.5")) + 
  geom_ribbon(rib2, alpha = 0.3, fill="#F8766D") +
  theme_minimal() +
  labs(x = "Number of evaluations", y = "Best MSE") + # color="legend.list2"
  theme(axis.text = element_text(size = 55),
        axis.title.x = element_text(size = 75),
        axis.title.y = element_text(size = 75),
        axis.line = element_line(color="gray20", size=1),
        legend.position = c(.782, .85), # c(.782, .85) c(.792, .85)
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.text = element_text(size=50),
        legend.background = element_rect(color = NA, fill = "transparent")) +
  geom_vline(xintercept = 150, linetype="dashed", size=2, color="gray40") +
  scale_x_continuous(breaks=round(seq(150,500,by=50),1)) +
  scale_color_manual(breaks=legend.list, values=legend.list2) 

################################################################################
# ggplot - Estimated vs True Importance
################################################################################

setwd(paste0(work.path,"/Results"))

weight_df = read.table(file="importance_boxplot.txt",header=FALSE,sep="|")
colnames(weight_df) = c("example","method","parameter","mean","sd")

p = ggplot(weight_df[weight_df$example=="Example III",], aes(x=parameter, y=mean, fill=method)) 
p +
  geom_bar(stat="identity", color="black", position=position_dodge(0.5), width = 0.45) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3, size=1,
                position=position_dodge(0.5)) +
  labs(x="Parameter", y = "Importance measure") +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('#F8766D','#00BFC4'), labels= c("Ex-Sliced-Calib", "Sobol")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 55),
        axis.title.x = element_text(size = 75),
        axis.title.y = element_text(size = 75),
        axis.line = element_line(color="gray20", size=1),
        legend.position = c(0.3, 0.88), # c(.782, .85) c(.792, .85)
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.text = element_text(size=55),
        legend.background = element_rect(color = NA))

