
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
  
  trial_dimension = trial_dim 
  tau0 = trial_tau0 
  p = trial_p
  
  major_matrix = matrix(NA,end,trial_dimension)
  minor_matrix = matrix(NA,end,D-trial_dimension)
  
  ## initialization
  Tk = matrix(NA,end,D)
  weight = matrix(NA,end,D)
  
  ## prior-posterior hyperparameters
  nu = matrix(NA,end,D)
  tau = matrix(NA,end,D)
  
  ## sampling distribution (likelihood) parameters 
  mu = matrix(NA,end,D)
  tau0 = tau0
  
  nu[ninit,] = rep(0.5,D)
  tau[ninit,] = rep(0.0001,D)
  mu[ninit,] = rep(0.5,D)
  
  for(ii in (ninit+1):end) {     
    print(ii)
    
    ## outlier eliminations 
    tryCatch({
      
      Q1 = summary(y)[2]
      Q3 = summary(y)[5]
      IQR = Q3-Q1
      no_outlier = as.logical( (y >= Q1-1.5*IQR)*(y <= Q3+1.5*IQR) )
      kmodel <- km(design=data.frame(X[no_outlier,]),response=y[no_outlier],control=list(trace=FALSE), nugget=1e-6) 

    }, error=function(e){cat("ERROR_Outlier:",conditionMessage(e), "\n")})
    
    ## sobol' sensitivity analysis 
    tryCatch({
    
      kriging.mean <- function(Xnew, m) predict.km(m, Xnew, "UK", se.compute=FALSE, checkNames=FALSE)$mean
      SA.metamodel <- fast99(model=kriging.mean, factors=D, n=500, q="qunif", q.arg=list(min=0, max=1), m=kmodel)
      main_effect = SA.metamodel$D1/SA.metamodel$V # main effect
      total_effect = 1 - SA.metamodel$Dt/SA.metamodel$V # total effect

      Tk[ii,] = total_effect
      
    }, error=function(e){cat("ERROR_Sobol :",conditionMessage(e), "\n")})
    
    ## update weight distributions using Thompson sampling
    tryCatch({
    
      for (k in 1:D) {
        
        nu[ii,k] = (tau[ii-1,k]*nu[ii-1,k] + tau0*Tk[ii,k]) / (tau[ii-1,k] + tau0)
        tau[ii,k] = tau[ii-1,k] + tau0
        mu[ii,k] = rnorm(1,mean=nu[ii,k],sd=1/sqrt(tau[ii,k]))
        weight[ii,k] = rtruncnorm(1,a=0.2,b=0.8,mean=mu[ii,k],sd=1/sqrt(tau0))
          
      }
      
    }, error=function(e){cat("ERROR_Thompson :",conditionMessage(e), "\n")})
    
    ## choose d parameters probabilistic  
    tryCatch({
      
      idx = rank(-weight[ii,], ties.method="random")
      major = c(1:D)[idx %in% c(1:trial_dimension)]
      minor = c(1:D)[-major]
        
      # slicing two subspaces
      XM <- X[,major] 
      Xm <- X[,minor] 
    
    }, error=function(e){cat("ERROR_Knapsack :",conditionMessage(e), "\n")})
    
    ## build a GP in the major subspace
    tryCatch({
      
      kmodelM = km(design=data.frame(XM), response=y, control=list(trace=FALSE), nugget=1e-6) 
      xnewM = max_EI(kmodelM, lower=rep(0,length(major)), upper=rep(1,length(major)),control=list(pop.size=200,max.generations=30,wait.generations=10,BFGSburnin=20,trace=FALSE))$par
    
    }, error=function(e){cat("ERROR_GP_Major :",conditionMessage(e), "\n")})
    
    ## sliced sequential design 
    tryCatch({
    
      if (runif(1)<p) {
        
        S = CandPoints(N=10000, p_cont=D-trial_dimension)
        res = MaxProAugment(ExistDesign=Xm, CandDesign=S, nNew=1)
        xnewm  = res$Design[nrow(res$Design),]
        
        xnew = rep(0,D)
        
        xnew[major] = xnewM
        xnew[minor] = xnewm   
        
        ynew = ff(obs.err, f(xnew,d=d,D=D), Nx)
        
      } else {
        
        xnew = rep(0,D)
        xnewm <- Xm[which.min(y),]
        
        xnew[major] = xnewM
        xnew[minor] = xnewm
        
        ynew <- ff(obs.err, f(xnew,d=d,D=D), Nx)
        
      }

    }, error=function(e){cat("ERROR_Sliced_Design :",conditionMessage(e), "\n")})
      
    ## store results 
    major_matrix[ii,] = major
    minor_matrix[ii,] = minor

    print(xnew)
    print(ynew)
    xnew = matrix(xnew,1,D)
    
    X <- rbind(X, xnew)
    y <- c(y, ynew)
  }
  
  return(list(X=X, y=y, major_matrix=major_matrix, minor_matrix=minor_matrix, Tk=Tk, weight=weight, nu=nu, tau=tau, mu=mu))
}
