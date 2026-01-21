
#################################################################
# test functions
#################################################################

ackley <- function(xx,d=d,D=D)
{
  ##########################################################################
  # ACKLEY FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # a = constant (optional), with default value 20
  # b = constant (optional), with default value 0.2
  # c = constant (optional), with default value 2*pi
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(0,...,0)
  # function type: Many local minima
  ##########################################################################
  
  a = 20 # default
  b = 0.2 # default
  c = 2*pi # default
  
  xx <- xx*15 - 5 
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  }  
  if (is.matrix(xx) == 1) {
    xxxx <- xx[,(d-1):D] # non-important parameters    
  } else {
    xxxx <- xx[(d-1):D] # non-important parameters  
  } 
  
  sum1 <- sum(xxx^2)
  sum2 <- sum(cos(c*xxx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  min_val = -5
  max_val = 10
  aa = min_val + 2*(max_val-min_val)/50
  bb = 4*1000
  cc = max_val - 2*(max_val-min_val)/50
  dd = 4*1000
  
  y <- term1 + term2 + a + exp(1)
  y <- y + ifelse(xxxx[1]>=min_val && xxxx[1]<=aa, -(bb/(aa-min_val))*(xxxx[1]-min_val) + bb, 0) + ifelse(xxxx[1]>=cc && xxxx[1]<=max_val, (cc/(max_val-cc))*(xxxx[1]-max_val) + cc, 0) + 
    ifelse(xxxx[2]>=min_val && xxxx[2]<=aa, -(bb/(aa-min_val))*(xxxx[2]-min_val) + bb, 0) + ifelse(xxxx[2]>=cc && xxxx[2]<=max_val, (cc/(max_val-cc))*(xxxx[2]-max_val) + cc, 0) +
    ifelse(xxxx[3]>=min_val && xxxx[3]<=aa, -(bb/(aa-min_val))*(xxxx[3]-min_val) + bb, 0) + ifelse(xxxx[3]>=cc && xxxx[3]<=max_val, (cc/(max_val-cc))*(xxxx[3]-max_val) + cc, 0)
  
  return(y)
}

levy <- function(xx, d, D) 
{
  ##########################################################################
  # LEVY FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(1,...,1)
  # function type: many local minima
  ##########################################################################
  
  xx <- (xx-0.5)*2*10
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  }
  if (is.matrix(xx) == 1) {
    xxxx <- xx[,(d-1):D] # non-important parameters    
  } else {
    xxxx <- xx[(d-1):D] # non-important parameters  
  } 
  
  w <- 1 + (xxx - 1)/4
  
  term1 <- (sin(pi*w[1]))^2 
  term3 <- (w[d]-1)^2 * (1+1*(sin(2*pi*w[d]))^2)
  
  wi <- w[1:(d-1)]
  sum <- sum((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
  
  min_val = -10
  max_val = 10
  aa = min_val + 2*(max_val-min_val)/50
  bb = 0.6*1000
  cc = max_val - 2*(max_val-min_val)/50
  dd = 0.6*1000
  
  y <- term1 + sum + term3 
  y <- y + ifelse(xxxx[1]>=min_val && xxxx[1]<=aa, -(bb/(aa-min_val))*(xxxx[1]-min_val) + bb, 0) + ifelse(xxxx[1]>=cc && xxxx[1]<=max_val, (cc/(max_val-cc))*(xxxx[1]-max_val) + cc, 0) + 
    ifelse(xxxx[2]>=min_val && xxxx[2]<=aa, -(bb/(aa-min_val))*(xxxx[2]-min_val) + bb, 0) + ifelse(xxxx[2]>=cc && xxxx[2]<=max_val, (cc/(max_val-cc))*(xxxx[2]-max_val) + cc, 0) +
    ifelse(xxxx[3]>=min_val && xxxx[3]<=aa, -(bb/(aa-min_val))*(xxxx[3]-min_val) + bb, 0) + ifelse(xxxx[3]>=cc && xxxx[3]<=max_val, (cc/(max_val-cc))*(xxxx[3]-max_val) + cc, 0) 
  return(y)
}

rastr <- function(xx, d, D)
{
  ##########################################################################
  # RASTRIGIN FUNCTION
  # input: 
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(0,...,0)
  # function type: many local minima
  ##########################################################################
  
  xx <- xx*7 - 3  
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  }
  
  sum <- sum(xxx^2 - 10*cos(2*pi*xxx))
  
  y <- (10*d + sum)
  return(y)
}

schwef <- function(xx, d, D)
{
  ##########################################################################
  # SCHWEFEL FUNCTION (log scale)
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(420.9687,...,420.9687) * log is recommended
  # function type: many local minima
  ##########################################################################
  
  xx <- (xx-0.5)*2*500
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  } 
  if (is.matrix(xx) == 1) {
    xxxx <- xx[,(d-1):D] # non-important parameters    
  } else {
    xxxx <- xx[(d-1):D] # non-important parameters  
  } 
  
  sum <- sum(xxx*sin(sqrt(abs(xxx)))) 
  
  min_val = -500
  max_val = 500
  aa = min_val + 2*(max_val-min_val)/50
  bb = 250*1000
  cc = max_val - 2*(max_val-min_val)/50
  dd = 250*1000
  
  y <- 418.9829*d - sum
  y <- y + ifelse(xxxx[1]>=min_val && xxxx[1]<=aa, -(bb/(aa-min_val))*(xxxx[1]-min_val) + bb, 0) + ifelse(xxxx[1]>=cc && xxxx[1]<=max_val, (cc/(max_val-cc))*(xxxx[1]-max_val) + cc, 0) + 
    ifelse(xxxx[2]>=min_val && xxxx[2]<=aa, -(bb/(aa-min_val))*(xxxx[2]-min_val) + bb, 0) + ifelse(xxxx[2]>=cc && xxxx[2]<=max_val, (cc/(max_val-cc))*(xxxx[2]-max_val) + cc, 0) +
    ifelse(xxxx[3]>=min_val && xxxx[3]<=aa, -(bb/(aa-min_val))*(xxxx[3]-min_val) + bb, 0) + ifelse(xxxx[3]>=cc && xxxx[3]<=max_val, (cc/(max_val-cc))*(xxxx[3]-max_val) + cc, 0)
  return(y)
}

perm0db <- function(xx, d, D)
{
  ##########################################################################
  # PERM FUNCTION 0, d, beta
  # input:
  # xx = c(x1, x2, ..., xD)
  # b = constant (optional), with default value 10
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(1,1/2,...,1/d)
  # function type: bowl-shaped
  ##########################################################################
  
  b = 10 # default
  
  xx <- (xx-0.5)*2*d
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  } 
  
  ii <- c(1:d)
  jj <- matrix(rep(ii,times=d), d, d, byrow=TRUE)
  
  xxmat <- matrix(rep(xxx,times=d), d, d, byrow=TRUE)
  inner <- rowSums((jj+b)*(xxmat^ii-(1/jj)^ii))	
  outer <- sum(inner^2)
  
  y <- outer
  return(y)
}

spheref <- function(xx, d, D)
{
  ##########################################################################
  # SPHERE FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # b = constant (optional), with default value 10
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(0,...,0)
  # function type: bowl-shaped
  ##########################################################################
  
  xx <- (xx-0.5)*2*5.12
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  } 
  
  sum <- sum(xxx^2)
  
  y <- (1/10)*(sum)
  return(y)
}

dixonpr <- function(xx, d, D)
{
  ##########################################################################
  # DIXON-PRICE FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt_i=2^{-(2^i-2)/2^i}, i=1,...,d
  # function type: valley-Shaped
  ##########################################################################
  
  xx <- (xx-0.5)*2*10
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  }
  
  x1 <- xxx[1]
  term1 <- (x1-1)^2
  
  ii <- c(2:d)
  xi <- xxx[2:d]
  xold <- xxx[1:(d-1)]
  sum <- sum(ii * (2*xi^2 - xold)^2)
  
  y <- term1 + sum
  return(y)
}

rosen <- function(xx, d, D)
{
  ##########################################################################
  # ROSENBROCK FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(1,...,1)
  # function type: valley-Shaped
  ##########################################################################
  
  xx <- xx*4 - 2  # (xx-0.5)*2*2.048
  
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  }
  
  xi <- xxx[1:(d-1)]
  xnext <- xxx[2:d]
  
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  
  y <- (sum)
  return(y)
}

michal <- function(xx, d, D)
{
  ##########################################################################
  # MICHALEWICZ FUNCTION
  # xx = c(x1, x2, ..., xD)
  # m = constant (optional), with default value 10
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=-1.8013 at d=2, xopt=(2.20,1.57)
  # global minimum: f(xopt)=-4.687658 at d=5, f(xopt)=-9.66015 at d=10
  # function type: steep ridges/drops
  ##########################################################################
  
  m = 10 # default
  
  xx <- xx*pi
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  } 
  
  ii <- c(1:length(xxx))
  sum <- sum(sin(xxx) * (sin(ii*xxx^2/pi))^(2*m))
  
  y <- -sum
  return(y)
}

branin <- function(xx, d, D)
{
  ##########################################################################
  # BRANIN FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # a = constant (optional), with default value 1
  # b = constant (optional), with default value 5.1/(4*pi^2)
  # c = constant (optional), with default value 5/pi
  # r = constant (optional), with default value 6
  # s = constant (optional), with default value 10
  # t = constant (optional), with default value 1/(8*pi)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0.397887 at xopt=(-pi,12.275),(pi,2.275),(9.42478,2.475)
  # function type: others
  ##########################################################################
  
  a=1; b=5.1/(4*pi^2); c=5/pi; r=6; s=10; t=1/(8*pi) # default
  
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:2] # important parameters    
  } else {
    xxx <- xx[1:2] # important parameters  
  } 
  
  xxx[1] = 15*xxx[1]-5
  xxx[2] = 15*xxx[2]
  
  x1 <- xxx[1]
  x2 <- xxx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  return(y)
}

goldpr <- function(xx, d, D)
{
  ##########################################################################
  # GOLDSTEIN-PRICE FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=3 at xopt=(0,-1)
  # function type: others
  ##########################################################################
  
  xx <- 4*xx-2
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:2] # important parameters    
  } else {
    xxx <- xx[1:2] # important parameters  
  } 
  
  x1 <- xxx[1]
  x2 <- xxx[2]
  
  fact1a <- (x1 + x2 + 1)^2
  fact1b <- 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2
  fact1 <- 1 + fact1a*fact1b
  
  fact2a <- (2*x1 - 3*x2)^2
  fact2b <- 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2
  fact2 <- 30 + fact2a*fact2b
  
  y <- fact1*fact2
  return(y)
}

hart3 <- function(xx, d, D) {
  
  ##########################################################################
  # HARTMANN 3-DIMENSIONAL FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=-3.86278 at xopt=(0.114614,0.555649,0.852547)
  # function type: others
  ##########################################################################  
  
  xx <- xx
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:3] # important parameters    
  } else {
    xxx <- xx[1:3] # important parameters  
  } 
  if (is.matrix(xx) == 1) {
    xxxx <- xx[,(d-1):D] # non-important parameters    
  } else {
    xxxx <- xx[(d-1):D] # non-important parameters  
  } 
  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  
  A <- c(3.0, 10, 30,
         0.1, 10, 35,
         3.0, 10, 30,
         0.1, 10, 35)
  
  A <- matrix(A, 4, 3, byrow=TRUE)
  P <- 10^(-4) * c(3689, 1170, 2673,
                   4699, 4387, 7470,
                   1091, 8732, 5547,
                   381, 5743, 8828)
  P <- matrix(P, 4, 3, byrow=TRUE)
  
  xxmat <- matrix(rep(xxx,times=4), 4, 3, byrow=TRUE)
  inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
  outer <- sum(alpha * exp(-inner))
  
  min_val = 0
  max_val = 1
  aa = min_val + 2*(max_val-min_val)/50
  bb = 0.4*1000
  cc = max_val - 2*(max_val-min_val)/50
  dd = 0.4*1000
  
  y <- -outer + 3.86278
  y <- y + ifelse(xxxx[1]>=min_val && xxxx[1]<=aa, -(bb/(aa-min_val))*(xxxx[1]-min_val) + bb, 0) + ifelse(xxxx[1]>=cc && xxxx[1]<=max_val, (cc/(max_val-cc))*(xxxx[1]-max_val) + cc, 0) + 
    ifelse(xxxx[2]>=min_val && xxxx[2]<=aa, -(bb/(aa-min_val))*(xxxx[2]-min_val) + bb, 0) + ifelse(xxxx[2]>=cc && xxxx[2]<=max_val, (cc/(max_val-cc))*(xxxx[2]-max_val) + cc, 0) +
    ifelse(xxxx[3]>=min_val && xxxx[3]<=aa, -(bb/(aa-min_val))*(xxxx[3]-min_val) + bb, 0) + ifelse(xxxx[3]>=cc && xxxx[3]<=max_val, (cc/(max_val-cc))*(xxxx[3]-max_val) + cc, 0)
  return(y)
}

hart6 <- function(xx, d, D)
{
  ##########################################################################
  # HARTMANN 6-DIMENSIONAL FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=-3.32237 at xopt=(0.20169,0.150011,0.476874,0.275332,0.311652,0.6573)
  # function type: others
  ##########################################################################
  
  xx <- xx
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:6] # important parameters    
  } else {
    xxx <- xx[1:6] # important parameters  
  } 
  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(10, 3, 17, 3.5, 1.7, 8,
         0.05, 10, 17, 0.1, 8, 14,
         3, 3.5, 1.7, 10, 17, 8,
         17, 8, 0.05, 10, 0.1, 14)
  A <- matrix(A, 4, 6, byrow=TRUE)
  P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
                   2329, 4135, 8307, 3736, 1004, 9991,
                   2348, 1451, 3522, 2883, 3047, 6650,
                   4047, 8828, 8732, 5743, 1091, 381)
  P <- matrix(P, 4, 6, byrow=TRUE)
  
  xxmat <- matrix(rep(xxx,times=4), 4, 6, byrow=TRUE)
  inner <- rowSums(A[,1:6]*(xxmat-P[,1:6])^2)
  outer <- sum(alpha * exp(-inner))
  
  y <- - outer + 3.32237 + ifelse(xxx[1]>=0 && xxx[1]<=0.01, 3000, 0) + ifelse(xxx[1]>=0.98 && xxx[1]<=1, 3000, 0)
  # y <- -(2.58 + outer) / 1.94
  
  return(y)
}

hart4 <- function(xx, d, D)
{
  ##########################################################################
  # HARTMANN 4-DIMENSIONAL FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=-3.32237 at xopt=(0.20169,0.150011,0.476874,0.275332)
  # function type: others
  ##########################################################################
  
  xx <- xx
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:4] # important parameters    
  } else {
    xxx <- xx[1:4] # important parameters  
  } 
  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(10, 3, 17, 3.5, 1.7, 8,
         0.05, 10, 17, 0.1, 8, 14,
         3, 3.5, 1.7, 10, 17, 8,
         17, 8, 0.05, 10, 0.1, 14)
  A <- matrix(A, 4, 6, byrow=TRUE)
  P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
                   2329, 4135, 8307, 3736, 1004, 9991,
                   2348, 1451, 3522, 2883, 3047, 6650,
                   4047, 8828, 8732, 5743, 1091, 381)
  P <- matrix(P, 4, 6, byrow=TRUE)
  
  xxmat <- matrix(rep(xxx,times=4), 4, 4, byrow=TRUE)
  inner <- rowSums(A[,1:4]*(xxmat-P[,1:4])^2)
  outer <- sum(alpha * exp(-inner))
  
  y <- (1.1 - outer) / 0.8387 + 3.135474
  
  return(y)
}

permdb <- function(xx, d, D)
{
  ##########################################################################
  # PERM FUNCTION d, beta
  # input:
  # xx = c(x1, x2, ..., xD)
  # b  = constant (optional), with default value 0.5
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(1,2,...,d)
  # function type: others
  ##########################################################################
  
  b=0.5
  
  xx <- xx*d
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  }
  
  ii <- c(1:d)
  jj <- matrix(rep(ii,times=d), d, d, byrow=TRUE)
  
  xxmat <- matrix(rep(xx,times=d), d, d, byrow=TRUE)
  inner <- rowSums((jj^ii+b)*((xxmat/jj)^ii-1))	
  outer <- sum(inner^2)
  
  y <- outer
  return(y)
}

powell <- function(xx, d, D)
{
  ##########################################################################
  # POWELL FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=0 at xopt=(0,...,0)
  # function type: others
  ##########################################################################
  
  xx <- xx*9-4
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  } 
  
  xxa <- xx[seq(1, d-3, 4)]
  xxb <- xx[seq(2, d-2, 4)]
  xxc <- xx[seq(3, d-1, 4)]
  xxd <- xx[seq(4, d, 4)]
  
  sumterm1 <- (xxa + 10*xxb)^2
  sumterm2 <- 5 * (xxc - xxd)^2
  sumterm3 <- (xxb - 2*xxc)^4
  sumterm4 <- 10 * (xxa - xxd)^4
  sum <- sum(sumterm1 + sumterm2 + sumterm3 + sumterm4)
  
  y <- sum
  return(y)
}

shekel <- function(xx, d, D)
{
  ##########################################################################
  # SHEKEL FUNCTION
  # input:
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=-10.1532 at xopt=(4,4,4,4) at m=5
  # global minimum: f(xopt)=-10.4029 at xopt=(4,4,4,4) at m=7
  # global minimum: f(xopt)=-10.5364 at xopt=(4,4,4,4) at m=10
  # function type: others
  ##########################################################################
  
  xx <- xx*10
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:4] # important parameters    
  } else {
    xxx <- xx[1:4] # important parameters  
  }  
  
  m <- 10 # default
  b <- 0.1 * c(1, 2, 2, 4, 4, 6, 3, 7, 5, 5)
  C <- c(4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0,
         4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6,
         4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0,
         4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6)
  C <- matrix(C, 4, 10, byrow=TRUE)
  Ct <- t(C)
  
  xxmat <- matrix(rep(xx,times=m), m, 4, byrow=TRUE)
  inner <- rowSums((xxmat-Ct[,1:4])^2)
  
  outer <- sum(1/(inner+b))
  
  y <- -outer
  return(y)
}

stybtang <- function(xx, d, D)
{
  ##########################################################################
  # STYBLINSKI-TANG FUNCTION
  # input: 
  # xx = c(x1, x2, ..., xD)
  # d = dimension of important parameters 
  # D = total dimension of parameters 
  # global minimum: f(xopt)=-39.16599d at xopt=(-2.903534,...,-2.903534)
  # function type: others
  ##########################################################################
  
  xx <- xx*10-5
  if (is.matrix(xx) == 1) {
    xxx <- xx[,1:d] # important parameters    
  } else {
    xxx <- xx[1:d] # important parameters  
  }
  
  sum <- sum(xxx^4 - 16*xxx^2 + 5*xxx)
  
  y <- sum/2
  return(y)
}

##########################################
ff <- function(x., y.,Nx) {
  return( (1/Nx)*sum( (x. + sqrt(y.))^2 ) )
}
##########################################

init_XY <- function(name, Ninit, D, d, Nx, seed) {
  
  Xlist = list() 
  ylist = list()
  errlist = list()
  
  if(name == 'ackley'){ 
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, ackley(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  } 
  
  if(name == 'levy'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.01)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, levy(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }
  
  if(name == 'rastr'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, rastr(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }
  
  if(name == 'schwef'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, schwef(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }
  
  if(name == 'perm0db'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, perm0db(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }  
  
  if(name == 'spheref'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, spheref(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }  
  
  if(name == 'dixonpr'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, dixonpr(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }   
  
  if(name == 'rosen'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, rosen(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }    
  
  if(name == 'michal'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, michal(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }    
  
  if(name == 'branin'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, branin(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }  
  
  if(name == 'goldpr'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, goldpr(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }  
  
  if(name == 'hart6'){ 
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, hart6(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }    
  
  if(name == 'hart4'){ 
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, hart4(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }    
  
  if(name == 'hart3'){ 
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.01)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, hart3(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }    
  
  if(name == 'permdb'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, permdb(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }   
  
  if(name == 'powell'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, powell(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }   
  
  if(name == 'shekel'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, shekel(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  } 
  
  if(name == 'stybtang'){  
    
    ################################
    # observational error 
    obs.err = rnorm(Nx,mean=0,sd=0.1)
    ################################
    # InitialDesign<-MaxProLHD(n = Ninit, p = D)$Design 
    # DOX<-MaxPro(InitialDesign)
    # X <- DOX$Design; colnames(X) = paste0("X",1:D)
    X <- maximinLHS(Ninit, D); colnames(X) = paste0("X",1:D)
    y <- numeric(0)
    for (i in 1:Ninit){
      temp = ff(obs.err, stybtang(matrix(X[i,],nrow=1), d, D), Nx)
      y[i] = temp
    }
    Xlist = X
    ylist = y
    errlist = obs.err
    
  }
  
  return(list(Xlist=Xlist, ylist=ylist, errlist=errlist)) 
}



