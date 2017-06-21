
# This function estimates the fixed parameters in the 
# joint models, by maximizing the profile h-likelihood function.
### Gradient OK.  (Jan 24, 2017)

estFixeff <- function(RespLog=list(Jlik1, Jlik2), 
                          fixed,     # pars to be estimated    
                          str_val,   # starting values of fixed pars
                          Dpars=Jraneff,   
                          long.data, surv.data, # dataset
                          B, Bi,     # random effect values
                          invSIGMA, sigma,    # dispersion pars
                          ParVal=NULL, # other par values, given
                          lower=NULL,
                          Silent=T){ 
  
  # derive -H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars=Dpars)
  q <- length(Dpars)
  n <- nrow(Bi)
  p <- length(fixed)
  
  
  # ff() returns the negative value of the adjusted 
  # profile h-likelihood to be optimized.
  
  ff <- function(xx){
    fy <- numeric(1)
    
    # assign values to parameters
    par.val <- Vassign(fixed, xx)    
    par.val <- data.frame(c(par.val, ParVal, sigma))
    
    # evaluate the H matrix
    nD1 <- evalMat(Hmats[[1]], q, long.data, par.val, raneff.val=B)   # longitudinal part
    nD2 <- evalMat(Hmats[[2]], q, surv.data, par.val, raneff.val=Bi)  # survival part
    nD3 <- -invSIGMA*n
    H <-  -(nD1+nD2+nD3)
    
    # evaluate h-likelihood
    lhlike1 <- with(long.data, with(par.val, with(B, eval(parse(text=RespLog[[1]])))))
    lhlike2 <- with(surv.data, with(par.val, with(Bi, eval(parse(text=RespLog[[2]])))))
    lhlike3 <-   - diag(0.5*as.matrix(Bi)%*%invSIGMA%*%t(Bi))
    
    # evaluate profile h-likelihood
    fy <- sum(lhlike1)+sum(lhlike2)+sum(lhlike3) -0.5*log(det(H/2/pi)) 
    
    return(-fy)
  }
  
  
  # gr() returns the gradients of the negative adjusted 
  # profile h-likelihood.
  
  gr <- function(xx){
    fy <- numeric(p)
    
    # assign values to parameters
    par.val <- Vassign(fixed, xx)    
    par.val <- data.frame(c(par.val, ParVal, sigma))
#     if(grepl("leftint", RespLog[[1]])){
#       long.data$leftint <- with(par.val, with(B, 
#                                               with(long.data, pnorm(delim_val, mean=eval(parse(text=Cmu)),
#                                                                     sd=eval(parse(text=Csigma))))))
#     }
    
    
    # evaluate inverse H matrix
    nD1 <- evalMat(Hmats[[1]], q, long.data, par.val, raneff.val=B)   # longitudinal part
    nD2 <- evalMat(Hmats[[2]], q, surv.data, par.val, raneff.val=Bi)  # survival part
    nD3 <- -invSIGMA*n
    H <-  as.matrix(-(nD1+nD2+nD3))
    invH <- solve(H)  
    
    # calculate and evaluate dH/dbeta, where beta are fixed par
    for(i in 1:p){
      dH <- dH2 <- rep(NA, q^2)
      
      for(j in 1:(q^2)){
        # from longitudinal models
        Dij <- Hmats[[1]][j]
        dH[j] <- Simplify(Vderiv(Dij, fixed[i]))
        # from survival models
        dH2[j] <- Vderiv(Hmats[[2]][j], fixed[i])
      }
      
      (dH_val1 <- evalMat(dH, q, long.data, par.val, raneff.val=B))
      (dH_val2 <- evalMat(dH2, q, surv.data, par.val, raneff.val=Bi))
      dH_val <- -(dH_val1 +dH_val2)  # note: Hmats returns -H matrix
      
      
      # calculate dh/dbeta: derivative of log h-likelihood  
      df1 <- Vderiv(RespLog[[1]], fixed[i])
      df1_val <- with(long.data, with(par.val, with(B, eval(parse(text=df1)))))
      
      df2 <- Vderiv(RespLog[[2]], fixed[i])
      df2_val <- with(surv.data, with(par.val, with(Bi, eval(parse(text=df2)))))
      
      fy[i] <- sum(df1_val)+sum(df2_val)-0.5*matrix.trace(as.matrix(invH%*%dH_val))
    }
    
    return(-fy)
  }
  
  
  ## start iteration
  str_val0 <- str_val
  message <- -1
  M <- 1

  if(Silent==F) check=0  else check=1
  
  if(is.null(lower)){
    lower= -Inf; method="BFGS"
  } else {
    method="L-BFGS-B"
  }
  
  # start iteration
  while(message != 0 & M<10){
    str_val0 <- sapply(str_val0, function(x)x+rnorm(1,0, min(1, abs(x/5))))
    
  
    result <- try(optim(par=str_val0, fn=ff, gr=gr,
                          method=method,
                          lower=lower,
                          control = list(
                            trace=1-check,
                            maxit=1000
                          )),
                    silent=T)
    
    
    error_mess <- attr(result, "class")
    
    if(length(error_mess)!=0 ){ 
      message = -1
    } else {
      message <- result$convergence 
      str_val0 <- result$par
    }
    
    if(Silent==F){ print(message); print(M); print(result)   }
    M <- M +1
  }
  
  
  if(message==0){
    gamma <- result$par
    names(gamma) <- names(str_val)
    fval <- result$value
    
#     cat("my gr:", gr(gamma), '\n')
#     eps=10^{-10}
#     numGr <- c()
#     for(i in 1:length(gamma)){
#       dist <- rep(0, length(gamma))
#       dist[i] <- eps
#       numGr <- c(numGr, (ff(gamma+dist)-ff(gamma))/eps)    
#     }
#     cat("numerical gr:", numGr, '\n')
    
#     plot(gr(gamma), ylim=range(c(gr(gamma), numGR)))
#     points(numGR, color="red")
    
  } else {
    stop("Iteration stops because fixed parameters can not be successfully estimated.")  
  }
  
  return(list(gamma=gamma, fval=fval))
}



