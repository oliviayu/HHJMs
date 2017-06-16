
# This function returns a matrix with elements in string.
# Moreover, it returns L(l)'L(l)=SIGMA, which is a spherical
# parameterization with diagonal elements 1

strMat <- function(q2){
    L <- c()
    L <- rbind(L, c(1, rep(0,q2-1)))
    Mpar <- c()
    
    for( i in 2:q2){
      # l0 <- paste('L',i,1, sep="") 
      l0 <- 1
      Li <- c()
      for(j in 1:i){
        m0 <- paste('L',i, 2:min((j+1),i),sep="")
        Mpar <- c(Mpar, m0)
        
        if(j < i){
          sin0 <- rep(c("sin(", "cos("), c(length(m0)-1, 1))
        } else {
          sin0 <- rep("sin(",  length(m0))
        }
        l1 <- paste(sin0, m0, rep(")", length(m0)), sep="")
        Li <- c(Li, paste(c(l0, l1), collapse = "*"))
      }
      L <- rbind(L, c(Li,rep(0, q2-i)))
    }
    L <- t(L)

  
  M <- matrix(NA,q2,q2)  # M=L'L
  for(i in 1:q2){
    for(j in 1:q2){
      M[i,j] <- Simplify(paste(L[,i], L[,j],sep="*",collapse = "+"))
    }
  }
  diag(M) <- "1"
  
  Mpar <- unique(Mpar)
  Mp <- length(Mpar)
  
  dM <- array(NA, dim=c(q2,q2,Mp))
  
  for(i in 1:Mp){
    ttz <- unlist(lapply(M, function(x){Vderiv(x, Mpar[i])}))
    dM[,,i] <- matrix(as.character(ttz),q2,q2)
  }
  
  
  return(list(M=M, Mpar=Mpar, dM=dM))
  
}