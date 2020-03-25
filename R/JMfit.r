library(Deriv)
library(matrixcalc)
library(lbfgs)

#' Fit a joint model
#' 
#' @param glmeObject A list, specifying the GLME models to be fitted.
#' @param survObject A list, specifying the survival model (either Cox PH or Weibull model) to be fitted.
#' @param long.data longitudinal data
#' @param surv.data survival data
#' @param idVar subject id
#' @param eventTime observed event time
#' @param survFit an object returned by coxph() or survreg() function to represent a fitted survival model from the two-step method
#' @param method a vector indicating which method to apply. If method='h-likelihood' (by default), the h-likelihood method is used; if method='aGH', the adaptive Gauss-Hermite method is used.
#' @param itertol Convergence tolerance on the relative absolute change in log-likelihood function between successive iterations. Convergence is declared when the change is less than itertol. Default is itertol = 0.001.
#' @param Ptol Convergence tolerance on the average relative absolute change in parameter estimates between successive iterations. Convergence is declared when the change is less than Ptol. Default is Ptol = 0.01.
#' @param epsilon A small numerical value, used to calculate the numerical value of the derivative of score function. The default value is 1e-6.
#' @param iterMax The maximum number of iterations. The default value is 10.
#' @param ghsize The number of quadrature points used in the adaptive GH method. The default value is 4.
#' @param Silent logical: indicating if messages about convergence success or failure should be suppressed.
#' 
#' @export
JMfit <- function(
  #### arguments, must specified by users  ####
  glmeObject, survObject, long.data, surv.data, idVar, eventTime, survFit,
  # arguments with defualt values
  method = c("h-likelihood", "aGH"), itertol = 1e-3, Ptol = 1e-2, epsilon = 1e-6,
  iterMax = 10,  ghsize = 4, srcpath = NULL, parallel = F, Silent = T){
  
  method <- match.arg(method)
  
  if(method == "h-likelihood"){
    result <- JMfit_HL(glmeObject, survObject, long.data, surv.data,  
                       idVar, eventTime, survFit, itertol, Ptol,iterMax, Silent)
  } else if (method == "aGH"){
    result <- JMfit_aGH(glmeObject, survObject, 
                        long.data, surv.data,  
                        idVar, eventTime, survFit, 
                        itertol=itertol, Ptol=Ptol, epsilon=epsilon,
                        iterMax, ghsize, srcpath=srcpath, 
                        parallel=parallel, Silent)
    
  } else{
    message("Must specify the method. If method=\"h-likelihood\" (by default), the h-likelihood method is used for joint inference. If method=\"aGH\", the adaptive GH method is used.")
    break
  }
  
  result
}