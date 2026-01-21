
ff <- function(x., y., Nx) {
  
  rslt = numeric(0)
  for (jj in c(1:length(y.))) {
    rslt[jj] = (1/Nx)*sum( (x. + sqrt(y.[jj]))^2 )
  }
  return( rslt )
}

schwef <- function(xx_temp, d, D)
{
  xx_temp = as.matrix(xx_temp)
  
  rslt = numeric(0)
  for (ii in c(1:nrow(xx_temp))) {
    
    xx = xx_temp[ii,]
    
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
    # y <- y + ifelse(xxxx[1]>=min_val && xxxx[1]<=aa, -(bb/(aa-min_val))*(xxxx[1]-min_val) + bb, 0) + ifelse(xxxx[1]>=cc && xxxx[1]<=max_val, (cc/(max_val-cc))*(xxxx[1]-max_val) + cc, 0) + 
    #   ifelse(xxxx[2]>=min_val && xxxx[2]<=aa, -(bb/(aa-min_val))*(xxxx[2]-min_val) + bb, 0) + ifelse(xxxx[2]>=cc && xxxx[2]<=max_val, (cc/(max_val-cc))*(xxxx[2]-max_val) + cc, 0) +
    #   ifelse(xxxx[3]>=min_val && xxxx[3]<=aa, -(bb/(aa-min_val))*(xxxx[3]-min_val) + bb, 0) + ifelse(xxxx[3]>=cc && xxxx[3]<=max_val, (cc/(max_val-cc))*(xxxx[3]-max_val) + cc, 0)
    
    rslt[ii] = y
  }
  return(rslt)
}

ackley <- function(xx_temp,d=d,D=D)
{
  xx_temp = as.matrix(xx_temp)
  
  rslt = numeric(0)
  for (ii in c(1:nrow(xx_temp))) {
    
    xx = xx_temp[ii,]
    
    a = 20 # default
    b = 0.2 # default
    c = 2*pi # default
    
    xx <- xx*15 - 5 #xx*14 - 7 #     xx*15 - 5 (xx-0.5)*2*32.768
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
    # y <- y + ifelse(xxxx[1]>=min_val && xxxx[1]<=aa, -(bb/(aa-min_val))*(xxxx[1]-min_val) + bb, 0) + ifelse(xxxx[1]>=cc && xxxx[1]<=max_val, (cc/(max_val-cc))*(xxxx[1]-max_val) + cc, 0) + 
    #   ifelse(xxxx[2]>=min_val && xxxx[2]<=aa, -(bb/(aa-min_val))*(xxxx[2]-min_val) + bb, 0) + ifelse(xxxx[2]>=cc && xxxx[2]<=max_val, (cc/(max_val-cc))*(xxxx[2]-max_val) + cc, 0) +
    #   ifelse(xxxx[3]>=min_val && xxxx[3]<=aa, -(bb/(aa-min_val))*(xxxx[3]-min_val) + bb, 0) + ifelse(xxxx[3]>=cc && xxxx[3]<=max_val, (cc/(max_val-cc))*(xxxx[3]-max_val) + cc, 0)
    
    rslt[ii] = y
  }
  return(rslt)
}

hart3 <- function(xx_temp, d, D) {
  
  xx_temp = as.matrix(xx_temp)
  
  rslt = numeric(0)
  for (ii in c(1:nrow(xx_temp))) {
    
    xx = xx_temp[ii,]
    
    # xx <- xx
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
    # y <- y + ifelse(xxxx[1]>=min_val && xxxx[1]<=aa, -(bb/(aa-min_val))*(xxxx[1]-min_val) + bb, 0) + ifelse(xxxx[1]>=cc && xxxx[1]<=max_val, (cc/(max_val-cc))*(xxxx[1]-max_val) + cc, 0) + 
    #   ifelse(xxxx[2]>=min_val && xxxx[2]<=aa, -(bb/(aa-min_val))*(xxxx[2]-min_val) + bb, 0) + ifelse(xxxx[2]>=cc && xxxx[2]<=max_val, (cc/(max_val-cc))*(xxxx[2]-max_val) + cc, 0) +
    #   ifelse(xxxx[3]>=min_val && xxxx[3]<=aa, -(bb/(aa-min_val))*(xxxx[3]-min_val) + bb, 0) + ifelse(xxxx[3]>=cc && xxxx[3]<=max_val, (cc/(max_val-cc))*(xxxx[3]-max_val) + cc, 0)
    
    rslt[ii] = y
  }
  return(rslt)
}