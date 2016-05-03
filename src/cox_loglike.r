
# This function returns the log likelihood function of 
# the Cox model that shares some parameters with 
# the mixed effect models.

cox_loglike <- function(survObject, Jllike){
   model <- survObject$fm
   status <- survObject$event
   resp <- fmReverse(model)$resp  # response variable (event time)
   rvX <- fmReverse(model)$rvX  # baseline covariates
   p <- length(rvX)   # dimension of fixed parameters
   par <- paste(survObject$par, 1:p, sep="")  # name fixed parameters
   
  
#  if(survObject$sharedPar=='random_effect'){  # if random effects are used as predictors in the Cox model 
    raneff <- unlist(Jllike$raneff)  # get random effects from mixed effect models
    q <- length(raneff)  # dimension of random effects
#   } else {  # other cases will be implemented in the future
#     error_message <- 'Error Message: Only shared random effects models are currently available.'
#     return(error_message)
#   } 
   
  Asso <- paste("Asso",1:q, sep="")  # name the coefficients of the random effects
  # linear predictor in Cox model
  linear_pred <- paste(paste(par,"*", rvX, collapse="+", sep=""), 
                       paste(Asso, "*", raneff, collapse="+", sep=""), 
                       sep="+")
  # log likelihood                          
  loglike <- paste(status,"*( log(h0) +", linear_pred, ") - 
                   exp(", linear_pred, ")*Th0", sep="")
  
  # initial values for the parameters to be estimated 
  str_val <- survObject$str_val
  names(str_val) <- c(par, Asso)
  
  return(list(resp=resp, loglike=loglike, par=c(par,Asso), linear_pred=linear_pred,
              str_val=str_val))
}

