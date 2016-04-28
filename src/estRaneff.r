
# This function estimates the random effects in the 
# joint models, by maximizing the h-likelihood function.

estRaneff <- function(RespLog, 
                    raneff=Jraneff, 
                    long.data=subdat, surv.dat=subsurv.dat,
                    invSIGMA=invSIGMA0, sigma=sigma0, 
                    ParVal, # other par values, given 
                    Silent=T){ 
  
   # ff() returns the value of h-likelihood
  
      ff <- function(xx){  
        fy <- numeric(1)
        
        # assign values to parameters
        B.val <- Vassign(raneff, xx)              
        par_val <- c(B.val, ParVal, sigma)
      
        # values from longitudinal data
        val <- with(par_val, 
                    with(long.data, eval(parse(text=RespLog[[1]]))))
        # values from survival data
        val2 <- with(par_val,
                    with(surv.dat, eval(parse(text=RespLog[[2]]))))
        # values from the distribution of random effect
        Sig <-  -0.5*xx%*%invSIGMA%*%xx +0.5*log(det(invSIGMA))
        
        fy <- sum(val)+sum(val2)+Sig
        
        return(-fy)
        
      }
      
     # start iteration
      message <- 1
      M <- 0
      bval <- rep(0,length(raneff))
      
      while(message!=0 & M<50){
        error_mess="try-error"
        
        while(length(error_mess)!=0){
          result <- try( optim(bval,  ff, control = list(maxit = 500)), silent=F)
          error_mess <- attr(result, "class")
          bval <- rep(0, length(raneff))
        }
        
        message <- result$convergence
        M <- M+1
        if(Silent==F) cat(paste("M=",M,", message=", message,".",sep=""))
      }
    
     if(M<50){
       bi <- result$par
     } else {
       message("Iteration limit reached without convergence -- for random effects.")  
       break
     }
       
  return(bi)
}


