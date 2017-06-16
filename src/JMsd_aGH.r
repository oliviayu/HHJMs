
JMsd_aGH <- function(HHJMfit_obj, ghsize=4, Silent=T, epsilon=10^{-6}){
    q = ncol(HHJMfit_obj$Bi)  
    RespLog=HHJMfit_obj$RespLog
    Jraneff=HHJMfit_obj$Jraneff 
    long.data=HHJMfit_obj$long.data
    surv.data=HHJMfit_obj$surv.data
    Bi=HHJMfit_obj$Bi 
    invSIGMA0=solve(HHJMfit_obj$covBi)
    sigma=HHJMfit_obj$sigma  
    ParVal=HHJMfit_obj$fixedest
    idVar=HHJMfit_obj$idVar
    uniqueID=HHJMfit_obj$uniqueID
    n=nrow(surv.data)
    
    GHzsamp0 = mgauss.hermite(n=ghsize, mu=rep(0,q), sigma=diag(1/2, q))
    idsigma = cholHi_adaptGH(RespLog, Jraneff, 
                             long.data, surv.data,
                             Bi, invSIGMA0, 
                             sigma,  ParVal,
                             idVar, uniqueID, Silent)
    # generate GH samples by subject
    GHsample0 <- as.list(rep(NA,n))
    for(i in 1:n){
      GHsample0[[i]] = mgauss.hermite(n=ghsize, mu=as.numeric(Bi[i,]), sigma=idsigma[[i]])
    }
    
    GHsd2 = try(get_aGH_sd2(long.data, surv.data,
                            fixedest0=ParVal,
                            RespLog, p=length(ParVal), q, 
                            GHzsamp0, GHsample0,
                            Dpars=Jraneff, ghsize,
                            uniqueID, idVar, 
                            invSIGMA=invSIGMA0, epsilon,
                            otherval=sigma), silent = T)
   return(GHsd2)
}