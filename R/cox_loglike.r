
# This function returns the log likelihood function of 
# the Cox model that shares some parameters with 
# the mixed effect models.

cox_loglike <- function(survObject, weibPar=NULL){
   model <- survObject$fm
   status <- survObject$event
   resp <- fmReverse(model)$resp  # response variable (event time)
   rvX <- fmReverse(model)$rvX  # baseline covariates
   if(sum(toupper(rvX)=='NULL')>0 | sum(rvX=='1')>0){ 
    linear1 <- "0"
    par=c()
   }else{
     p <- length(rvX)   # dimension of fixed parameters
     par <- paste(survObject$par, 1:p-1, sep="")  # name fixed parameters
     linear1 <- paste(par,"*", rvX, collapse="+", sep="")
   }
   linear_pred <- linear1
  # log likelihood                
  if(is.null(survObject$distribution)){
    loglike <- paste(status,"*( log(h0) +", linear_pred, ") - 
                   exp(", linear_pred, ")*Th0", sep="")
    str_val <- survObject$str_val
    names(str_val) <- par
  } else if(survObject$distribution=="weibull"){

    # log-likelihood for weibull ph model
    loglike <- paste(status,"*( Wlogscale+log(Wshape)+log(",resp,")*(Wshape-1) +", linear_pred, ") - 
                   exp(Wlogscale+", linear_pred, ")*",resp,"^Wshape", sep="")
#     loglike <- paste(status,"*( Wlogscale+Wlogshape+log(",resp,")*(exp(Wlogshape)-1) +", linear_pred, ") - 
#                    exp(Wlogscale+", linear_pred, ")*",resp,"^exp(Wlogshape)", sep="")

    
    # log-likelihood for weibull AFT model
#     loglike <- paste(status,"*( Wlogscale+log(Wshape)+log(",resp,")*(Wshape-1) -(", linear_pred, 
#                      ")*Wshape) - exp(Wlogscale-(", linear_pred, ")*Wshape)*",
#                      resp,"^Wshape", sep="")
    
    if(is.null(weibPar)){ weibPar=c(0,1) }
    str_val <- c(survObject$str_val, weibPar)
    par=c(par,"Wlogscale", "Wshape")
    # par=c(par,"Wlogscale", "Wlogshape")
    names(str_val) <- par
  }
  
  return(list(resp=resp, loglike=loglike, par=par, linear_pred=linear_pred,
              str_val=str_val))
}

