
# This function estimates the fixed parameters in the 
# joint models, by maximizing the profile h-likelihood function.
### Gradient OK.  (Jan 24, 2017)

estFixeff_adptGH <- function(RespLog=list(Jlik1, Jlik2), 
                             fixed=names(fixedest0),     # pars to be estimated    
                             str_val=fixedest0,   # starting values of fixed pars
                             Dpars=Jraneff,   
                             idVar,
                             long.data, surv.data, # dataset
                             GHsample,
                             GHzsamp,
                             invSIGMA0,
                             uniqueID, 
                             ghsize,
                             ParVal=NULL, # other par values, given
                             lower, # lower/upper bound for dispersion parameters
                             upper,
                             Silent=T){ 
  
  # derive -H matrix, where H is defined in the adjusted profile h-likelihood
  q <- length(Dpars)
  p <- length(fixed)
  n = nrow(surv.data)
  L2  <- strMat(q)
  q0 <- length(L2$Mpar)
  weights = GHzsamp$weights
  
  # xx=str_val
  ff <- function(xx){
    fy <- numeric(1)
    
    # assign values to parameters
    par.val <- Vassign(fixed, xx[1:p])    
    par.val <- data.frame(c(par.val, ParVal))
    
    # invmat=SIGMA;
    Lval <- Vassign(L2$Mpar, xx[-(1:p)]) 
    invmat <- evalMat(as.list(L2$M), q, par.val=Lval)
    invSIGMA <- solve(invmat)
    
    
    # evaluate approximate log-likelihood
    fn = rep(NA, n)
    
    for(i in 1:n){
      # i=1
      subdat1 = subset(long.data, long.data[,idVar]==uniqueID[i])
      subdat2 = subset(surv.data, surv.data[,idVar]==uniqueID[i])
      samples = as.data.frame(GHsample[[i]]$points)
      names(samples)=Dpars
      likefn = rep(NA, ghsize^q)
      for(j in 1:ghsize^q){
        # j=2
        lhlike1 = with(subdat1, with(par.val, with(samples[j,], eval(parse(text=RespLog[[1]])))))
        lhlike2 = with(subdat2, with(par.val, with(samples[j,], eval(parse(text=RespLog[[2]])))))
        lhlike3 = exp(-as.matrix(samples[j,])%*%invSIGMA%*%t(samples[j,])/2+sum(GHzsamp$points[j,]^2))
        likefn[j] = exp(sum(lhlike1)+lhlike2)*lhlike3
      }
      
      fn[i] = log(sum(likefn*weights)) 
    }

    fy = sum(fn) + n/2*log(det(invSIGMA))
#     print(n/2*log(det(invSIGMA)))
#     print(sum(fn))
    
    return(-fy)
  }
  
  
  gr.long <- deriv(formula(paste("~", RespLog[[1]])), fixed)
  gr.surv <- deriv(formula(paste("~", RespLog[[2]])), fixed)
  
  gr <- function(xx){
    fy <- numeric(p+q0)
    
    # assign values to parameters
    par.val <- Vassign(fixed, xx[1:p])    
    par.val <- data.frame(c(par.val, ParVal))
    
    # invmat=SIGMA;
    Lval <- Vassign(L2$Mpar, xx[-(1:p)]) 
    invmat <- evalMat(as.list(L2$M), q, par.val=Lval)
    invSIGMA <- solve(invmat)
    

    gn = matrix(NA, nrow=n, ncol=p+q0)
    
    for(i in 1:n){
      # i=1
      subdat1 = subset(long.data, long.data[,idVar]==uniqueID[i])
      subdat2 = subset(surv.data, surv.data[,idVar]==uniqueID[i])
      samples = as.data.frame(GHsample[[i]]$points)
      names(samples)=Dpars
      
      likefn =  rep(NA, ghsize^q)
      gri = matrix(NA, ghsize^q, p+q0)
      
      norm_term = 0
      
      for(j in 1:ghsize^q){
        # j=1
        lhlike1 = with(subdat1, with(par.val, with(samples[j,], eval(parse(text=RespLog[[1]])))))
        lhlike2 = with(subdat2, with(par.val, with(samples[j,], eval(parse(text=RespLog[[2]])))))
        lhlike3 =exp(-as.matrix(samples[j,])%*%invSIGMA%*%t(samples[j,])/2)
        lhlike_all =exp(sum(lhlike1)+lhlike2)*lhlike3
        likefn[j] = lhlike_all*exp(sum(GHzsamp$points[j,]^2))
        norm_term = norm_term+lhlike_all
        
        val1 = with(par.val, with(subdat1, with(samples[j,], attr(eval(gr.long),"gradient"))))
        val2 = with(par.val, with(subdat2, with(samples[j,], attr(eval(gr.surv), "gradient"))))
        
        ## derivative w.r.t parameters in cov(b)
        ## check, correct!
        fy2 <- rep(NA, q0)
        for(s in 1:q0){
          dM <- L2$dM[,,s]
          dM_val <- evalMat(dM, q, par.val=Lval)
          fy2[s] <- -0.5* matrix.trace(invSIGMA%*%dM_val)+
                       0.5*diag(as.matrix(samples[j,])%*%invSIGMA%*%dM_val%*%invSIGMA%*%t(as.matrix(samples[j,])))
        }
        
        gri[j,] = c(apply(val1,2,sum)+val2, fy2)*likefn[j]
      }
      
      gn[i,] = as.vector(t(gri)%*%weights)/norm_term
    }
    

    fy = apply(gn,2,sum)*det(invSIGMA0)^{-1/2}*2^{q/2}
    
    return(-fy)
  }
  
  
  ## start iteration
  str_val0 = str_val
  Alower = c(lower, rep(0,q0)) 
  Aupper = c(upper, rep(pi, q0))
  # print(ff(c(str_val0, runif( q0, 0, pi))))
  # print(gr(c(str_val0, runif( q0, 0, pi))))
  
  message <- -1
  M <- 1

  if(Silent==F) check=0  else check=1
  
  # start iteration
  while(message != 0 & M<50){
    
    result <- try(optim(par=c(str_val0, runif( q0, 0, pi)), fn=ff, gr=gr,
                        method="L-BFGS-B",
                        lower=Alower,
                        upper=Aupper,
                        control = list(
                          trace=1-check,
                          maxit=500
                        ), hessian=F),
                  silent=T)
    # result
    error_mess <- attr(result, "class")
    
    if(length(error_mess)!=0 ){ 
      message = -1
    } else {
      message <- result$convergence 
    }
    str_val0 <- sapply(str_val, function(x)x+rnorm(1,0, min(1, abs(x/2))))
    
    if(Silent==F){ print(message); print(M); print(result)   }
    M <- M +1
  }
  
  
  if(message==0){
    gamma <- result$par[1:p]
    names(gamma) <- fixed
    fval <- result$value
    
    Lval <- Vassign(L2$Mpar, result$par[-(1:p)])
    invmat <- evalMat(as.list(L2$M), q, par.val=Lval)
    mat <- solve(invmat) # invSIGMA
    
#         str_val0 <- sapply(str_val0, function(x)x+rnorm(1,0, min(1, abs(x/5))))
#         gamma=str_val0
#         gamma=xx
#        myGR = gr(gamma)
#         cat("my gr:", myGR, '\n')
#         eps=10^{-10}
#         numGr <- c()
#         for(i in 1:length(gamma)){
#           dist <- rep(0, length(gamma))
#           dist[i] <- eps
#           numGr <- c(numGr, (ff(gamma+dist)-ff(gamma))/eps)
#           cat("i=",i,'\n')
#         }
#         cat("numerical gr:", numGr, '\n')
#     #     
#         plot(myGR, ylim=range(c(myGR, numGr)))
#         points(numGr, pch=6, col="red")
#         myGR-numGr
    
  } else {
    stop("Iteration stops because fixed parameters can not be successfully estimated.")  
  }
  
  return(list(gamma=gamma, fval=fval, invSIGMA=mat,
              Lval=Lval))
}



