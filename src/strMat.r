
# This function returns a matrix with elements in string.
# Moreover, it returns L(l)'L(l)=inv.SIGMA, which can be considered
# as Chloskey decomposition of inverse covariance matrix of 
# random effects.

strMat <- function(q2){
  
  L <- c()
  for(i in 1:q2){
    Li <- c(rep(0, i-1), paste("L", i, i:q2, sep=""))
    L <- rbind(L, Li)
  }
  
  M <- matrix(NA,q2,q2)  # M=L'L
  for(i in 1:q2){
    for(j in 1:q2){
      M[i,j] <- Simplify(paste(L[,i], L[,j],sep="*",collapse = "+"))
    }
  }
  
  Mpar <- t(L)[lower.tri(L,diag=T)]
  Mp <- sum(1:q2)

  dM <- array(NA, dim=c(q2,q2,Mp))
  
  for(i in 1:Mp){
    ttz <- unlist(lapply(M, function(x){Vderiv(x, Mpar[i])}))
    dM[,,i] <- matrix(as.character(ttz),q2,q2)
  }

  
  return(list(M=M, Mpar=Mpar, dM=dM))
  
}