
multiResultClass <- function(X=NULL,y=NULL,obs.err=NULL,Xrslt=NULL,yrslt=NULL,Xprog=NULL,yprog=NULL,Xbest=NULL,ybest=NULL,ybest_test=NULL,obs.err_test=NULL,time_record=NULL,major_rslt=NULL,minor_rslt=NULL,Tk_rslt=NULL,weight_rslt=NULL,nu_rslt=NULL,tau_rslt=NULL,mu_rslt=NULL)
{
  me <- list(
    X = X,
    y = y,
    obs.err = obs.err,
    Xrslt = Xrslt,
    yrslt = yrslt,
    Xprog = Xprog,
    yprog = yprog,
    Xbest = Xbest,
    ybest = ybest,
    ybest_test = ybest_test,
    obs.err_test = obs.err_test,
    time_record = time_record,
    major_rslt = major_rslt,
    minor_rslt = minor_rslt,
    Tk_rslt = Tk_rslt,
    weight_rslt = weight_rslt,
    nu_rslt = nu_rslt,
    tau_rslt = tau_rslt,
    mu_rslt = mu_rslt    
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

