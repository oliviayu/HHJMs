

# calculate hessian matrix 

getHessian <- function(lik1, pars){
  lik <- parse(text=lik1)
  q <- length(pars)
  result <- as.list(rep(NA, q^2))
  

    for(i in 1:q){
      for(j in 1:q){
        k <- (i-1)*q+j
        result[[k]] <- D(D(lik, pars[i]), pars[j]) 
        names(result)[k] <- paste(pars[i], pars[j],sep=",")
      }
    }

  return(result)
}  