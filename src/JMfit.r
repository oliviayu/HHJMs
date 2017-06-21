library(Deriv)
library(matrixcalc)
library(lbfgs)

JMfit <- function(
  #### arguments, must specified by users  ####
  glmeObject, 
  survObject,
  long.data, surv.data,  
  idVar, eventTime, survFit,
  
  # arguments by defualt
  method = "h-likelihood",
  itertol=1e-3, Ptol=1e-2, epsilon=10^{-6},
  iterMax=10,  ghsize=4, srcpath=NULL, parallel,
  Silent=T
){
  
  if(method=="h-likelihood"){
    result = JMfit_HL(glmeObject, survObject, long.data, surv.data,  
                     idVar, eventTime, survFit, itertol, Ptol,iterMax, Silent)
  } else if (method =="aGH"){
    result = JMfit_aGH(glmeObject, survObject, 
                       long.data, surv.data,  
                      idVar, eventTime, survFit, 
                      itertol=itertol, Ptol=Ptol, epsilon=epsilon,
                      iterMax, ghsize, srcpath=srcpath, 
                      parallel=parallel, Silent)

  } else{
    message("Must specify the method. If method=\"h-likelihood\" (by default), the h-likelihood method is used for joint inference. If method=\"aGH\", the adaptive GH method is used.")
    break
  }
  
  return(result)
}