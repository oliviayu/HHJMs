

HHJMsummary <- function(testjm, digits=4){
  est <- round(testjm$fixedest, digits)
  sd <- round(testjm$fixedsd, digits)
  zval <- round(testjm$fixedest/testjm$fixedsd, digits)
  pval <- round((1-pnorm(abs(zval), 0, 1))*2, digits)
  
  Coeff <- data.frame(Estimate=est, Std.Error=sd, z.value=zval, Pvalue=pval)
  print(Coeff)
}