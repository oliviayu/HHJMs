#' Estimate random effects
#' 
get_raneff <- function(RespLog, 
                       long.data,surv.data, idVar, uniqueID,
                       Jraneff, invSIGMA0, sigma0, n, ni,q,N,
                       fixedest0,
                       Silent, scale=T){
  
    nB <- nBi <- c()      
    
    for(i in 1:n){
      subdat <- subset(long.data, long.data[, idVar]==uniqueID[i])
      subsurv.dat <- subset(surv.data, surv.data[, idVar]==uniqueID[i])
      
      bi <- estRaneff(RespLog=RespLog, 
                      raneff=Jraneff, 
                      long.data=subdat, surv.dat=subsurv.dat, 
                      invSIGMA=invSIGMA0, sigma=sigma0, 
                      ParVal=fixedest0,
                      Silent=Silent)    
      nBi <- rbind(nBi, bi)
      nB <- rbind(nB, matrix(rep(bi, ni[i]), ncol=q, byrow=T))           
      if(Silent==F)  cat("i=",i, "\n")    
    }

  if(scale==T){
    Bi <- as.matrix(nBi)%*%diag(1/apply(nBi,2,sd),q,q)
    B <- as.matrix(nB)%*%diag(1/apply(nBi,2,sd),q,q)
    cenBi <- apply(Bi, 2, mean)
    Bi <- Bi - matrix(rep(1, n),ncol=1)%*%matrix(cenBi, nrow=1)
    B <- B - matrix(rep(1, N),ncol=1)%*%matrix(cenBi, nrow=1)
  } else {
    Bi <- nBi
    B <- nB
  }
    
  Bi <- as.data.frame(Bi)
  B <- as.data.frame(B)
  names(B) <- names(Bi) <- Jraneff
  rownames(Bi) <- rownames(B) <- c()
  
  return(output=list(Bi=Bi, B=B))
}