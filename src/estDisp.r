
# This function estimates the dispersion parameters in the 
# joint models, by maximizing the adjusted profile 
# h-likelihood function.

# condition on: getHmat(), Vassign(), evalMat()

estDisp <- function(RespLog=list(Jlik1, Jlik2), 
                       Jdisp,      # pars to be estimated    
                       sigma0,   # starting values 
                       invSIGMA0,  # starting values
                       Dpars=c(Jraneff,Jfixed,Spar),   
                       long.data, surv.data, # dataset
                       B, Bi,     # values of random effects, given
                       ParVal, # values of fixed par, given
                       Silent=T){ 
  
  # derive -H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars=Dpars)  
  q <- length(Dpars)
  q1 <- length(Jdisp)
  q2 <- ncol(Bi)
  n <- nrow(Bi)
  L2  <- strMat(q2)   # Return matrices of strings for cov(bi) 
  

  # ff() returns the negative value of the adjusted 
  # profile h-likelihood to be optimized.
  
  ff <- function(xx){  # xx are input values of dispersion parameters
    fy <- numeric(1)
    
    # assign values to the parameters
    par.val <- Vassign(Jdisp, xx[1:q1]^2)  # sigma must be >0 
    Lval <- Vassign(L2$Mpar, xx[-(1:q1)]) 
    mat <- evalMat(as.list(L2$M), q2, par.val=Lval)
    par.val <- data.frame(c(par.val, ParVal))
    
    # evaluate the H matrix
    nD1 <- evalMat(Hmats[[1]], q, long.data, par.val, raneff.val=B)   # longitudinal part
    nD2 <- evalMat(Hmats[[2]], q, surv.data, par.val, raneff.val=Bi)  # survival part
    nD3 <- bdiag(-mat*n, diag(0, q-q2,q-q2))  # random effect part
    H <-  as.matrix(-(nD1+nD2+nD3))

    # evaluate the h-likelihood
    lhlike1 <- with(long.data, with(par.val, with(B, eval(parse(text=RespLog[[1]])))))
    lhlike2 <- with(surv.data, with(par.val, with(Bi, eval(parse(text=RespLog[[2]])))))
    lhlike3 <-  0.5*log(det(mat)) - diag(0.5*as.matrix(Bi)%*%mat%*%t(as.matrix(Bi)))
    
    # evaluate the adjusted profile h-likelihood
    fy <- sum(lhlike1)+sum(lhlike2)+sum(lhlike3)-0.5*log(det(H/2/pi))
    
    return(-fy)
  }
  
  
  # gr() returns the gradients of the negative adjusted 
  # profile h-likelihood.
  
  gr <- function(xx){
    fy <- numeric(q1+sum(1:q2))
    
    # assign values to the parameters
    par.val <- Vassign(Jdisp, xx[1:q1]^2) 
    Lval <- Vassign(L2$Mpar, xx[-(1:q1)])
    mat <- evalMat(as.list(L2$M), q2, par.val=Lval)
    par.val <- data.frame(c(par.val, ParVal))
    
    # evaluate the inverse of H 
    nD1 <- evalMat(Hmats[[1]], q, long.data, par.val, raneff.val=B)   # longitudinal part
    nD2 <- evalMat(Hmats[[2]], q, surv.data, par.val, raneff.val=Bi)  # survival part
    nD3 <- bdiag(-mat*n, diag(0, q-q2,q-q2))
    H <-  as.matrix(-(nD1+nD2+nD3))
    invH <- ginv(H)
    

    ## part1: dispersion pars in mixed-effect models
    # evaluate d(H)/d(sigma), where sigma are dispersion parameters to be estimated
    for(i in 1:q1){
      dH  <- rep(NA, q^2) 
      
      for(j in 1:q^2){
        Dij <- Hmats[[1]][j]
        dH[j] <- Simplify(Vderiv(Dij, Jdisp[i]))
      }
      
      dH_val1 <- evalMat(dH, q, long.data, par.val, raneff.val=B)
      dH_val <- -(dH_val1)  
      
      ## dh /d sigma, where sigma is from LME models and 
      ## h is log h-likelihood function
      df1 <- Vderiv(RespLog[[1]], Jdisp[i])
      df1_val <- with(long.data, with(par.val, with(B, eval(parse(text=df1)))))
   
      fy0 <- sum(df1_val)-0.5*matrix.trace(as.matrix(invH%*%dH_val))
      fy[i] <- fy0*2*xx[i]
    }
    
    ## part2: dispersion pars of random effects
    # get dh/dL[i,j] -0.5 tr(H^(-1) %*% dH/d L[i,j])
    for(i in 1:length(L2$Mpar)){
      dMi <- L2$dM[,,i]
      dMival <- evalMat(mat=as.list(dMi), q=q2, par.val=Lval)
      
      # dh/dL[i,j]
      df_Mval <- sum(matrix.trace(ginv(mat)%*%dMival)/2 - 
                       diag(as.matrix(Bi)%*%dMival%*%t(as.matrix(Bi)))/2)
      # H^(-1) %*% dH/d L[i,j]
      dH_Mval <- invH%*%(bdiag(dMival,diag(0, q-q2,q-q2))*n )
      fy[q1+i] <- df_Mval-0.5*matrix.trace(as.matrix(dH_Mval))
    }
    
    return(-fy) 
  }
  
  # initial values 
  L <- chol(invSIGMA0)
  str_val <- unlist(c(sigma0, L[lower.tri(L,diag=T)]))
  
  message <- -1
  M <- 1
  
  
  if(Silent==F) check=0  else check=1
  
  # start iteration
  while(message < 0 & M<50){
    error_mess="try-error"
    k <- 0
    
    if(M<10){ str_val0 <- str_val  
    }else{ 
      L0 <- diag(1, q2, q2)
      str_val0 <- c(rep(1, q1), L0[lower.tri(L0,diag=T)])}
    

    while(length(error_mess)!=0 & k<50){

      result <- try(lbfgs::lbfgs(call_eval=ff, call_grad=gr,
                                 vars=str_val0, epsilon=1e-3, 
                                 delta=1e-4,
                                 max_iterations=1000,
                                 invisible = check), 
                    silent=T)
      
      # result <- BBoptim(str_val0, ff, gr, control=list(checkGrad=T))

      str_val0 <- sapply(str_val0, function(x){x+rnorm(1, 0, max(1, abs(x)))})
      
      error_mess <- attr(result, "class")
      
      k <- k+1
      if(Silent==F){
        cat('k=',k , error_mess, '\n')
        # print(result)
      } 
        
    }
    
    if(k>=50){ 
      stop('Error occurs when estimating dispersion parameters.')
    }
  
    message <- result$convergence
    str_val0 <- result$par
      
    if(Silent==F) print(message); print(M); print(result)

    M <- M+1
  }
  
    
  if(M<50){
    output <- (result$par[1:q1])^2
    names(output) <- Jdisp
      
    Lval <- Vassign(L2$Mpar, result$par[-(1:q1)])
    mat <- evalMat(as.list(L2$M), q2, par.val=Lval)
  } else {
    message("Iteration limit reached without convergence -- for dispersion parameters.")  
    stop()
  }
  
  return(list(sigma=output, invSIGMA=mat, Lval=Lval))
  
}


