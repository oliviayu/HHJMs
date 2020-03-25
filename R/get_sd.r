#' Estimate the standard errors of fixed parameters
#' 
get_sd <- function(long.data, surv.data, Bi, B,
                   fixedest0, sigma0, invSIGMA0,
                   RespLog, Jfixed, Jraneff, p, q, n){
  
  finalHmat <- getHmat(RespLog, pars=c( Jfixed, Jraneff))
  mat1 <- evalMat(finalHmat$negH_long, q=p+q, data = long.data, 
                  par.val =c(sigma0, fixedest0), raneff.val = B) 
  mat2 <- evalMat(finalHmat$negH_surv, q=p+q, data = surv.data, 
                  par.val = c(sigma0, fixedest0), 
                  raneff.val = Bi)
  mat3 <- bdiag(diag(0, p), -invSIGMA0*(n))
  Hval <-  as.matrix(-(mat1+mat2+mat3))
  covMat <- ginv(Hval)
  sd2 <- diag(covMat)[1:p]
  sd <- sqrt(sd2)
  
  return(sd)
}
