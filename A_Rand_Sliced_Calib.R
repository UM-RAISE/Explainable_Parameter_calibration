
bov <- function(y, end=length(y))
{
  prog <- rep(min(y), end)
  prog[1:min(end, length(y))] <- y[1:min(end, length(y))]
  for(i in 2:end) 
    if(is.na(prog[i]) || prog[i] > prog[i-1]) prog[i] <- prog[i-1]
  return(prog)
}

optim.EI <- function(f, ninit, end, obs.err, trial_dim, trial_tau0, trial_p, X, y, D, d, Nx, name, seed)
{

  set.seed(seed)
  
  major_matrix = matrix(NA,end,trial_dim)
  minor_matrix = matrix(NA,end,D-trial_dim)
  
  Tk = matrix(NA,end,D)
  weight = matrix(NA,end,D)
  nu = matrix(NA,end,D)
  tau = matrix(NA,end,D)
  mu = matrix(NA,end,D)
  
  for(ii in (ninit+1):end) {     
    print(ii)
    
    # target dimension
    target_dim = trial_dim 

    major = sort(sample(1:D,target_dim,replace=FALSE,prob=rep(1,D)/D))
    minor = (1:D)[-major]
    
    # slicing two subspaces
    XM <- X[,major] 
    Xm <- X[,minor] 
    
    major_matrix[ii,] = major
    minor_matrix[ii,] = minor
    
    # major space design
    tryCatch({
      kmodelM <- km(design=data.frame(XM),response=y,control=list(trace=FALSE), nugget=1e-6) 
      xnewM <- max_EI(kmodelM,lower=rep(0,length(major)),upper=rep(1,length(major)),control=list(pop.size=200,max.generations=30,wait.generations=10,BFGSburnin=20,trace=FALSE))$par # 
    }, error=function(e){cat("ERROR_GP_Major :",conditionMessage(e), "\n")})
    
    # minor space design
    xnewm <- Xm[which.min(y),]
    
    xnew = rep(0,D)
    xnew[major] = xnewM
    xnew[minor] = xnewm

    xnew = matrix(xnew,1,D)
    ynew <- ff(obs.err, f(xnew,d=d,D=D), Nx)  
    print(xnew)
    print(ynew)    
    
    X <- rbind(X, xnew)
    y <- c(y, ynew)
  }
  
  return(list(X=X, y=y, major_matrix=major_matrix, minor_matrix=minor_matrix, Tk=Tk, weight=weight, nu=nu, tau=tau, mu=mu))
}
