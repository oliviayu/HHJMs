
# This function returns the conditional log-likelihood of 
# a longitudinal variable, y|u, modelled by a regular GLME model.

# condition on: fmReverse()

glme_loglike <- function(glmeObject){
  
  vars <- fmReverse(glmeObject$fm)
  resp <- vars$resp  # response variable
  rvX <- vars$rvX  # covariance for fixed pars
  rvZ <- vars$rvZ # covariance for random effects
  p <- length(rvX)  # dimension of fixed pars
  q <- length(rvZ) # dimension of random effects

  fixed <- paste(glmeObject$par, 0:(p-1), sep="") # name fixed pars
  raneff <- paste(glmeObject$ran.par, 1:q, sep="")  # name random effects
  linear_pred <- paste(c(fixed, raneff), rep("*", p+q), 
                       c(rvX, rvZ), sep="", collapse="+")
  
  if(glmeObject$family=="binomial"){          
    
    loglike <- paste(resp, "*(", linear_pred, ")-", 
                     "log(1+exp(", linear_pred, "))")
    sigma=c()
    
  } else if (glmeObject$family=="normal"){ 
    
    sigma <- glmeObject$sigma
    loglike <- paste("- 0.5*(", resp, "- (", linear_pred, 
                     "))^2/",sigma, "^2-log(",sigma,")-0.5*log(2*pi)",
                     sep="")
    
  } else if (glmeObject$family=="poisson"){
    
    loglike <- paste(resp, "*(", linear_pred, ")-exp(",
                     linear_pred, ")-log(factorial(",resp,"))")
    sigma=c()
    
  }

  return(list(loglike=loglike, fixed=fixed, raneff=raneff, 
              linear_pred=linear_pred,
              rvX=rvX, rvZ=rvZ,
              resp=resp, sigma=sigma))
}

