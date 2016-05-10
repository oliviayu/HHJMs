
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
    # evaluate d(H)/d(xi), sigma=xi^2, where sigma are dispersion parameters to be estimated
    for(i in 1:q1){
      dH  <- rep(NA, q^2) 
      
      for(j in 1:q^2){
        Dij <- Hmats[[1]][j]
        dH[j] <- Simplify(Vderiv(Dij, Jdisp[i]))
      }
      
      dH_val1 <- evalMat(dH, q, long.data, par.val, raneff.val=B)
      dH_val <- -(dH_val1)*2*xx[i]  

      ## dh /d xi, with sigma=xi^2 and sigma is from LME models and 
      ## h is log h-likelihood function
      df1 <- Vderiv(RespLog[[1]], Jdisp[i])
      df1_val <- with(long.data, with(par.val, with(B, eval(parse(text=df1)))))
   
      fy[i] <- sum(df1_val)*2*xx[i]-0.5*matrix.trace(as.matrix(invH%*%dH_val))
      
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
  L <- t(chol(solve(cov(Bi))))
  str_val0 <- unlist(c(sqrt(sigma0), L[lower.tri(L,diag=T)]))
  message <- (-1)
  M <- 1
  
  if(Silent==F) check=0  else check=1
  
  # start iteration
  while(message < 0 & M<20){
  
    if(M<10){
        result <- try(lbfgs::lbfgs(call_eval=ff, call_grad=gr,
                                   vars=str_val0, epsilon=1e-4, 
                                   delta=1e-4,
                                   max_iterations=1500,
                                   invisible = check), 
                      silent=T)
    } else {
      result <- try(BB::BBoptim(str_val0, fn=ff, gr=NULL, 
                                control=list(checkGrad=F,
                                             ftol=1e-10,
                                             gtol=1e-2,
                                             trace=T)),
                    silent=T)
    }
        
    error_mess <- attr(result, "class")
        
    if(length(error_mess)!=0 ){ 
      message = -1
    } else {
      if(M<10) message <- result$convergence 
      if(M>=10) message <- -abs(result$convergence)
      str_val0 <- result$par
    }
     
    str_val0 <- sapply(str_val0, function(x)x+rnorm(1,0, min(1, abs(x/5))))
    if(Silent==F) print(message); print(M); print(result)   
    M <- M +1
  }

    
  if(message==0){
    output <- (result$par[1:q1])^2
    names(output) <- Jdisp
    Lval <- Vassign(L2$Mpar, result$par[-(1:q1)])
    mat <- evalMat(as.list(L2$M), q2, par.val=Lval)
  } else {
    stop("Iteration limit reached without convergence -- for dispersion parameters.")  
  }
  
  return(list(sigma=output, invSIGMA=mat, Lval=Lval))
  
}


