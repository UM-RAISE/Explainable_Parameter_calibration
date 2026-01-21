
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
  
  if (name == "schwef") {
    nugget_val = 10^-4
  } else {
    nugget_val = 10^-6
  }
  
  kmodelM = km(design=data.frame(X),response=y,control=list(trace=FALSE),nugget=nugget_val)
  rslt = TREGO.nsteps(model=kmodelM,fun=f,nsteps=end-ninit,lower=rep(0,D),upper=rep(1,D),d=d,D=D,trace=3,control=list(pop.size=200,max.generations=30,wait.generations=10,BFGSburnin=20,trace=FALSE)) 
    
  X = rbind(X, rslt$par)
  y = c(y, rslt$value)
  
  return(list(X=X, y=y, major_matrix=major_matrix, minor_matrix=minor_matrix, Tk=Tk, weight=weight, nu=nu, tau=tau, mu=mu))
}
