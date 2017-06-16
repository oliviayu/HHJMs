
# This function estimates the fixed parameters in the 
# joint models, by maximizing the profile h-likelihood function.
### Gradient OK.  (Jan 24, 2017)

cholHi_adaptGH <- function(RespLog=list(Jlik1, Jlik2), 
                      Jraneff,   
                      long.data, surv.data, # dataset
                      Bi,     # random effect values
                      invSIGMA, sigma,    # dispersion pars
                      ParVal=NULL, # other par values, given
                      idVar,
                      uniqueID,
                      Silent=T){ 
  
  # derive -H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars=Jraneff)
  q <- length(Jraneff)
  n <- nrow(Bi)
  par.val <- data.frame(as.list(c(ParVal, sigma)))
  
  nD3 <- -invSIGMA
  
  Chol.H = as.list(rep(NA,n))
  
  for(i in 1:n){
    subdat1 = subset(long.data, long.data[,idVar]==uniqueID[i])
    subdat2 = subset(surv.data, surv.data[,idVar]==uniqueID[i])
  
    # evaluate the H matrix
    nD1 <- evalMat(Hmats[[1]], q, subdat1, par.val, raneff.val=Bi[i,])   # longitudinal part
    nD2 <- evalMat(Hmats[[2]], q, subdat2, par.val, raneff.val=Bi[i,])  # survival part
    H <-  -(nD1+nD2+nD3)
    
    # Chol.H[[i]] = chol(H)
    Chol.H[[i]] = solve(H)
  }

  return(Chol.H)

}



