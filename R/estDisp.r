#' Estimate the dispersion parameters in the joint models, by maximizing the 
#' adjusted profile h-likelihood function.
estDisp <- function(RespLog = list(Jlik1, Jlik2), 
                    Jdisp,      # pars to be estimated    
                    sigma0,   # starting values 
                    invSIGMA0,  # starting values
                    Dpars = c(Jraneff,Jfixed,Spar),   
                    long.data, surv.data, # dataset
                    B, Bi,     # values of random effects, given
                    ParVal, # values of fixed par, given
                    Silent = T,
                    Lval0 = NULL, lower, upper){ 
  
  # derive -H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars=Dpars)  
  q <- length(Dpars)
  q1 <- length(Jdisp)
  q2 <- ncol(Bi)
  n <- nrow(Bi)
  L2  <- strMat(q2)   # Return matrices of strings for cov(bi) 
  
  # ff() returns the negative value of the adjusted 
  # profile h-likelihood to be optimized.
  # xx <- c(unlist(sigma0),0.5)
  ff <- function(xx){  # xx are input values of dispersion parameters
    fy <- numeric(1)
    
    # assign values to the parameters
    if(!is.null(Jdisp)){
      par.val <- Vassign(Jdisp, xx[1:q1])  
      Lval <- Vassign(L2$Mpar, xx[-(1:q1)]) 
    } else { 
      par.val <- c()
      Lval <- Vassign(L2$Mpar, xx) 
    }
    
    # invmat=SIGMA 
    invmat <- evalMat(as.list(L2$M), q2, par.val=Lval)
    mat <- solve(invmat)
    par.val <- as.list(c(par.val, ParVal))
    
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
    # cat("fy value=", fy,'\n')
    return(-fy)
  }
  
  
  ### gradient function
  q0 <- length(L2$Mpar)
  k <- q1+q0
  dhlike1 <- deriv(formula(paste("~",RespLog[[1]])), Jdisp)
  dhlike2 <- deriv(formula(paste("~",RespLog[[2]])), Jdisp)
  dH1 <- dH2 <- as.list(rep(NA, q^2))
  for(i in 1:q^2){
    dH1[[i]] <- deriv(formula(paste("~",Hmats[[1]][i])), Jdisp)
    dH2[[i]] <- deriv(formula(paste("~",Hmats[[2]][i])), Jdisp)
  }
  
  ############ gradient function
  gr <- function(xx){  # xx are input values of dispersion parameters
    fy <- numeric(k)
    
    # assign values to the parameters
    if(!is.null(Jdisp)){
      par.val <- Vassign(Jdisp, xx[1:q1])  # sigma must be >0 
      Lval <- Vassign(L2$Mpar, xx[-(1:q1)]) 
    } else { 
      par.val <- c()
      Lval <- Vassign(L2$Mpar, xx) 
    }
    # invmat=SIGMA 
    invmat <- evalMat(as.list(L2$M), q2, par.val=Lval)
    mat <- solve(invmat)
    par.val <- as.list(c(par.val, ParVal))
    
    # evaluate the H matrix
    nD1 <- evalMat(Hmats[[1]], q, long.data, par.val, raneff.val=B)   # longitudinal part
    nD2 <- evalMat(Hmats[[2]], q, surv.data, par.val, raneff.val=Bi)  # survival part
    nD3 <- bdiag(-mat*n, diag(0, q-q2,q-q2))  # random effect part
    H <-  as.matrix(-(nD1+nD2+nD3))
    invH <- solve(H)
    # derivative of sd parameters 
    val <- with(par.val, with(long.data, with(B, attr(eval(dhlike1),"gradient"))))
    dh_val1 <- apply(val, 2, sum)
    val2 <- with(par.val, with(surv.data, with(Bi, attr(eval(dhlike2), "gradient"))))
    dh_val2 <- apply(val2, 2, sum)
    myDmat <- myDmat2 <- matrix(NA, nrow=q^2, ncol=q1)
    for(i in 1:q^2){
      myDmat[i,] <- apply(with(par.val, with(long.data, with(B, attr(eval(dH1[[i]]),"gradient")))),2,sum)
      myDmat2[i,] <- apply(with(par.val, with(surv.data, with(Bi, attr(eval(dH2[[i]]),"gradient")))),2,sum)
    }
    traces <- rep(NA, q1)
    for(i in 1:q1){
      dH1_val <- matrix(myDmat[,i], q,q)
      dH2_val <- matrix(myDmat2[,i], q,q)
      traces[i] <-  matrix.trace(invH%*%(-dH1_val-dH2_val)) 
    }
    
    fy0 <- dh_val1+dh_val2 -0.5*traces
    fy1 <- - fy0
    
    # derivative of correction term
    fy2 <- rep(NA, q0)
    for(i in 1:q0){
      dM <- L2$dM[,,i]
      dM_val <- evalMat(dM, q2, par.val=Lval)
      gh1 <- sum(-0.5* matrix.trace(mat%*%dM_val)+
                   0.5*diag(as.matrix(Bi)%*%mat%*%dM_val%*%mat%*%t(as.matrix(Bi))))
      DH <- mat%*%dM_val%*%mat*n
      DH2 <- -bdiag(DH, diag(rep(0,q-q2)))
      gh2 <- -0.5*matrix.trace(as.matrix(invH%*%DH2))
      fy2[i] <- -(gh1+gh2)
    }
    
    fy <- c(fy1, fy2)
    # cat("gr value=", fy,'\n')
    return(fy)
  }
  
  
  message <- (-1)
  M <- 0
  
  if(Silent==F) check=0  else check=1
  if((q1+q2)==1){ method1=method2="Brent"} else{
    method="L-BFGS-B"
  }
  
  # print(L2$M)
  
  # start iteration
  while(message != 0 & M<20){
    if(is.null(Lval0)){
      Lval00= runif( q0, 0, pi)
    } else{
      Lval00 =  Lval0+rnorm(length(Lval0),0,0.01)
    }
    str_val0 <-  unlist(c(sigma0+rnorm(length(sigma0),0,0.1), 
                          Lval00))
    
    
    result <- try(optim(par=str_val0, fn=ff, gr=gr,
                        method=method,
                        lower=c(lower, rep(0,q0)), 
                        upper=c(upper,rep(pi,q0)),
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
      str_val <- result$par
    }
    
    if(Silent==F){ print(message); print(M); print(result)   }
    M <- M +1
  }
  
  
  if(message==0){
    if(q1>0){
      output <- result$par[1:q1]
      names(output) <- Jdisp
      Lval <- Vassign(L2$Mpar, result$par[-(1:q1)])
      invmat <- evalMat(as.list(L2$M), q2, par.val=Lval)
      mat <- solve(invmat)
    } else{ 
      output=NULL
      Lval <- Vassign(L2$Mpar, result$par)
      invmat <- evalMat(as.list(L2$M), q2, par.val=Lval)
      mat <- solve(invmat)
    }
    
  } else {
    stop("Iteration stops because dispersion parameters can not be successfully estimated. 
         The Hessian matrix of random effects and fixed parameters may not be invertible. 
         It's probably because the random effects are highly correlated with each other.")  
  }
  
  return(list(sigma=output, invSIGMA=mat, Lval=Lval))
}


