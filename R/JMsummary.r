#' Summarize modelling results in table
#' 
#' @export
JMsummary <- function(testjm, newSD = NULL, digits = 3){
  est <- testjm$fixedest
  if(is.null(newSD)){
    sd <- testjm$fixedsd
  } else {
    sd <- newSD
  }
  Zvalue <- est/sd
  Pvalue <- (1-pnorm(abs(Zvalue), 0, 1))*2
  Coeff <- data.frame(Estimate = est, Std.Error = sd, Zvalue, Pvalue)
  print(round(Coeff,digits = digits))
}