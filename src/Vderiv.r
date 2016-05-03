
# compute derivatives of an expression w.r.t 
# a vector of parameters respectively

Vderiv <- function(lik1, pars){
  
    lik <- parse(text=lik1)
    q <- length(pars)
    result <- as.list(rep(NA, q))
       
    for(i in 1:q){
      result[[i]] <- Simplify(D(lik, pars[i])) 
      names(result)[i] <- pars[i]
    }
        

  return(result)
}  