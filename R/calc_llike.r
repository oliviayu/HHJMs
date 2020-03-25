
calc_llike <- function(long.data, surv.data, RespLog,
                       Bi, B, fixedest, new_sigma, 
                       new_invSIGMA,Jraneff,q,n){
  
    log_val1 <- sum(with(long.data, with(B, with(as.list(c(fixedest,new_sigma)),
                eval(parse(text=RespLog[[1]]))))))+sum(with(surv.data, with(Bi, with(
                as.list(c(fixedest, new_sigma)), eval(parse(text=RespLog[[2]]))))))
    log_val2 <- sum(-diag(as.matrix(Bi)%*%new_invSIGMA%*%t(as.matrix(Bi)))/2+
                0.5*log(det(new_invSIGMA))-q/2*log(2*pi))  
    hloglike_value <- log_val1+log_val2
    
    Hmats <- getHmat(RespLog, pars=Jraneff)
    nH1 <- evalMat(Hmats[[1]], q, long.data, 
                   par.val=c(fixedest, new_sigma), raneff.val=B)   # longitudinal part
    nH2 <- evalMat(Hmats[[2]], q, surv.data, 
                   par.val=c(fixedest, new_sigma), raneff.val=Bi)  # survival part
    nH3 <- -new_invSIGMA*n
    Hval <-  -(nH1+nH2+nH3)
    loglike_value <- hloglike_value -0.5*log(det(Hval/2/pi))
    
    return(loglike_value)
    
}