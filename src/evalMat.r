
# evaluate a matrix 

evalMat <- function(mat, q, data=NULL, par.val=NULL, raneff.val=NULL){
  D <- matrix(0, q, q)
  
  if(!is.data.frame(par.val)){
    par.val <- data.frame(as.list(par.val))
  }
  
  for(i in 1:q){
    for(j in 1:q){
      kk <- (i-1)*q+j
      Di <- with(data, with(par.val, with(raneff.val,  
                eval(parse(text=mat[kk])))))

      D[i,j] <- sum(Di)
     }
  }
  
  return(D)
}