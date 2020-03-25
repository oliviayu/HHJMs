
## derive H matrix in profile h-likelihood functions
## density of random effects are excluded here, will be included while estimating parameters
## Note: this function returns -H 

getHmat <- function(RespLog, pars){
  
      p <- length(pars)
      
      negH_long <-  getHessian(RespLog[[1]], pars)
      negH_surv <-  getHessian(RespLog[[2]], pars)
      
  return(list(negH_long=negH_long, negH_surv=negH_surv))
}


