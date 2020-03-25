library(foreach)
library(parallel)
library(doParallel)

get_aGH_sd2 <- function(long.data, surv.data,
                       fixedest0, 
                       RespLog, p, q, GHzsamp,GHsample,
                       Dpars=Jraneff,ghsize,
                       uniqueID, idVar, invSIGMA,
                       epsilon=10^{-6}, otherval=NULL, srcpath=NULL, parallel=F){
  
  #############
  fixed = names(fixedest0)
  p = length(fixedest0)
  gr.long <- deriv(formula(paste("~", RespLog[[1]])), fixed)
  gr.surv <- deriv(formula(paste("~", RespLog[[2]])), fixed)
  
  weights = GHzsamp$weights
  n = nrow(surv.data)
  
  ## function for calculating S(theta, bi)
  gr <- function(xx){
    fy <- numeric(p)
    # assign values to parameters
    par.val <- Vassign(fixed, xx[1:p])
    par.val = c(par.val, otherval)
    # invmat=SIGMA;
    
    gn = matrix(NA, nrow=n, ncol=p)
    
    for(i in 1:n){
      # i=1
      subdat1 = subset(long.data, long.data[,idVar]==uniqueID[i])
      subdat2 = subset(surv.data, surv.data[,idVar]==uniqueID[i])
      samples = as.data.frame(GHsample[[i]]$points)
      names(samples)=Dpars
      
      likefn =  rep(NA, ghsize^q)
      gri = matrix(NA, ghsize^q, p)
      
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
        

        gri[j,] = c(apply(val1,2,sum)+val2)*likefn[j]
      }
      
      gn[i,] = as.vector(t(gri)%*%weights)/c(norm_term)
    }
    
    fy = apply(gn,2,sum)*det(invSIGMA)^{-1/2}*2^{q/2}
    return(-fy)
  }
  
  ## esimate Hessian matrix by numerical derivative 
  est = unlist(c(fixedest0))
  
  if(parallel==T){
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    
    Hmat = foreach(exponent = 1:(p), 
                   .combine = rbind
    )  %dopar%  
    {
      setwd(srcpath)
      file.sources = list.files(pattern="*.r$")
      sapply(file.sources,source,.GlobalEnv)
      
      Delta = rep(0, p)
      Delta[exponent] = epsilon
      gr1= gr(est+Delta)
      gr2= gr(est-Delta)
      (gr1-gr2)/(epsilon*2)
    }
    
    stopCluster(cl)
  } else {
    Hmat = matrix(NA, nrow=p, ncol=p)
    for(i in 1:(p)){
      # i=1
      Delta = rep(0, p)
      Delta[i] = epsilon
      Hmat[i,]=(gr(est+Delta)-gr(est))/epsilon
      # cat("i=",i,'\n')
    }
  }
  
  sd = sqrt(diag(solve(Hmat)))
  
  return(sd)
}