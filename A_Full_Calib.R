
bov <- function(y, end=length(y))
{
  prog <- rep(min(y), end)
  prog[1:min(end, length(y))] <- y[1:min(end, length(y))]
  for(i in 2:end) 
    if(is.na(prog[i]) || prog[i] > prog[i-1]) prog[i] <- prog[i-1]
  return(prog)
}

eps <- sqrt(.Machine$double.eps)

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
  
  tryCatch({
    kmodel <- km(design=data.frame(X),response=y,control=list(trace=FALSE),nugget=1e-6) 
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})   
  
  ## optimization loop of sequential acquisitions
  maxei <- c()
  for(ii in (ninit+1):end) {
    print(ii)
    tryCatch({
      xnew <- max_EI(kmodel,lower=rep(0,ncol(X)),upper=rep(1,ncol(X)),control=list(pop.size=200,max.generations=30,wait.generations=10,BFGSburnin=20,trace=FALSE))$par # 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})   
    
    tryCatch({
      ynew <- ff(obs.err, f(xnew,d=d,D=D), Nx)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
    
    print(xnew)
    print(ynew)
    xnew = matrix(xnew,1,D)

    tryCatch({
    kmodel <- update(object = kmodel, newX = xnew, newy = ynew, cov.reestim = TRUE)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
    
    X <- rbind(X, xnew)
    y <- c(y, ynew)
    
  }
  
  return(list(X=X, y=y, major_matrix=major_matrix, minor_matrix=minor_matrix, Tk=Tk, weight=weight, nu=nu, tau=tau, mu=mu))
}
