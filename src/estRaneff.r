
# This function estimates the random effects in the 
# joint models, by maximizing the h-likelihood function.
## gradient OK! (Jan 24, 2017)

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
        # print(val)
        # values from survival data
        val2 <- with(par_val,
                    with(surv.dat, eval(parse(text=RespLog[[2]]))))
        # print(val2)
        # values from the distribution of random effect
        Sig <-  -0.5*xx%*%invSIGMA%*%xx +0.5*log(det(invSIGMA))
        
        fy <- sum(val)+sum(val2)+Sig
        # print(fy)
        return(-fy)
      }
      
      
      k= length(raneff)
#       gr.long <- Vderiv(RespLog[[1]], raneff)
#       gr.surv <- Vderiv(RespLog[[2]], raneff)
      gr.long <- deriv(formula(paste("~", RespLog[[1]])), raneff)
      gr.surv <- deriv(formula(paste("~", RespLog[[2]])), raneff)
      
      gr <- function(xx){
        fy <- numeric(k)
        # assign values to parameters
        B.val <- Vassign(raneff, xx)              
        par_val <- c(B.val, ParVal, sigma)
        gr.val1 <- gr.val2 <- rep(NA, k)
        
        ##
        val <-  with(par_val, with(long.data, attr(eval(gr.long),"gradient")))
        # print(val)
        gr.val1 <- as.vector(apply(val, 2, sum))
        #cat("long gr:", gr.val1, '\n')
        # gradient values from survival data
        val2 <-  with(par_val, with(surv.dat, attr(eval(gr.surv), "gradient")))
        #cat("surv gr:", val2,'\n')
        gr.val2 <- as.vector(val2)
        
#         for(i in 1:k){
#           # gradients from longitudinal data
#           val <-  with(par_val, with(long.data, eval(parse(text=gr.long[[i]]))))
#           # print(val)
#           gr.val1[i] <- sum(val)
#           # gradient values from survival data
#           val2 <-  with(par_val, with(surv.dat, eval(parse(text=gr.surv[[i]]))))
#           gr.val2[i] <- val2
#         }
        
        # print("check ok")
        # gradient values from the distribution of random effect
        gr.val3 <-  as.vector(-invSIGMA%*%xx)
        
#         print(gr.val1)
#         print(gr.val2)
#         print(gr.val3)
        fy <- gr.val1+gr.val2+gr.val3
        # print(fy)

        return(-fy)
      }
      
      
     # start iteration
      message <- 1
      M <- 0
      bval <- rep(0,length(raneff))
      
      while(message!=0 & M<50){
        error_mess="try-error"
        
        result <- try( optim(par=bval,  fn=ff, gr=gr,
                             method="L-BFGS-B",
                            control = list(maxit = 2000,
                                           trace=0)), silent=T)
#           result <- try(BB::BBoptim(bval, fn=ff, gr=gr, 
#                                     control=list(checkGrad=T,
#                                    trace=T)),
#                             silent=F)
#         
        error_mess <- attr(result, "class")

          
        M <- M+1
#         bval <- rnorm(length(raneff),0,1)
#         ff(bval)
#         gr(bval)
#         bval
#         result
        if(length(error_mess)==1) message=-1 else message <- result$convergence

        if(Silent==F) cat(paste("M=",M,", message=", message,".",sep=""))
      }
    
     if(message==0){
       bi <- result$par
#        cat("my gr:", gr(bi), '\n')
#        eps=10^{-10}
#        cat("numerical gr:", c((ff(bi+c(eps,0,0,0))-ff(bi))/eps,
#                (ff(bi+c(0,eps,0,0))-ff(bi))/eps,
#                (ff(bi+c(0,0,eps,0))-ff(bi))/eps,
#                (ff(bi+c(0,0,0,eps))-ff(bi))/eps
#              ),'\n')
     } else {
       stop("Iteration stops because random effects can not be successfully estimated.")  
     }
       
      
  return(bi)
}


