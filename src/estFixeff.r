
# This function estimates the fixed parameters in the 
# joint models, by maximizing the profile h-likelihood function.

estFixeff <- function(RespLog=list(Jlik1, Jlik2), 
                          fixed,     # pars to be estimated    
                          str_val,   # starting values of fixed pars
                          Dpars=Jraneff,   
                          long.data, surv.data, # dataset
                          B, Bi,     # random effect values
                          invSIGMA, sigma,    # dispersion pars
                          ParVal=NULL, # other par values, given
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
    fy <- sum(lhlike1)+sum(lhlike2)+sum(lhlike3)-0.5*log(det(H/2/pi))
    
    return(-fy)
  }
  
  
  # gr() returns the gradients of the negative adjusted 
  # profile h-likelihood.
  
  gr <- function(xx){
    fy <- numeric(p)
    
    # assign values to parameters
    par.val <- Vassign(fixed, xx)    
    par.val <- data.frame(c(par.val, ParVal, sigma))
    
    # evaluate inverse H matrix
    nD1 <- evalMat(Hmats[[1]], q, long.data, par.val, raneff.val=B)   # longitudinal part
    nD2 <- evalMat(Hmats[[2]], q, surv.data, par.val, raneff.val=Bi)  # survival part
    nD3 <- -invSIGMA*n
    H <-  as.matrix(-(nD1+nD2+nD3))
    invH <- ginv(H)  
    
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
  
  while(message < 0 & M<50){
    error_mess="try-error"
    k <- 0
    
    while(length(error_mess)!=0 & k<50){
      str_val0 <- sapply(str_val0, function(x)x+rnorm(1,0, min(1, abs(x/5))))
      
      result <- try(lbfgs::lbfgs(call_eval=ff, call_grad=gr,
                          vars=str_val0, #epsilon=1e-4, 
                         # delta=1e-4,
                          max_iterations=1500,
                          invisible = check), 
                    silent=T)
      

#       result <- try(BBoptim(str_val0, ff, gr, control=list(checkGrad=T)),
#                     silent=T)
      
      error_mess <- attr(result, "class")
      
      k <- k+1
    }
    
    if(k>=50){ 
      stop('Error occurs when estimating fixed parameters.')
    }
    
    message <- result$convergence
    str_val0 <- result$par
    if(Silent==F) print(message); print(M); print(result)
    M <- M+1
  }
  
  if(M<50){
    gamma <- result$par
    names(gamma) <- names(str_val)
    fval <- result$value
  } else {
    message("Iteration limit reached without convergence -- for fixed parameters.")  
    stop()
  }
  
  return(list(gamma=gamma, fval=fval))
}



