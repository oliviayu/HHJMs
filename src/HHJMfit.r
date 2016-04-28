## newly updated on Apr 22nd, 2016

# This function is for joint modelling of multiple 
# (generalized) linear mixed-effect models and a Cox model, 
# via h-likelihood method. It also addresses left-censoring 
# problem in a continuous longitudinal data, e.g. NAb.

library(Deriv)
library(matrixcalc)
library(lbfgs)

HHJMfit <- function(
      #### arguments must be specified by users  ####
      glmeObject, 
      survObject,
      long.data, surv.data,  # data 
      idVar, 

      # arguments by defualt
      SIGMA=NULL,   # starting values
      iterVal=0.05, 
      iterMax=10, 
      nblock=100, # for estimating baseline hazard
      Silent=T  # check outputs, for debugging
      ){
  
  ########################## settings for Longitudinal modeles
  Jllike <- Jglme_loglike(glmeObject)  
  Jraneff <- unlist(Jllike$raneff) # random effects in glme's
  Jfixed <- unlist(Jllike$fixed)   # fixed parameters in glme's
  Jdisp <- na.omit(unlist(Jllike$sigma))  # dispersion parameters
  # joint log-likelihood of longitudinal variables, simplified
  (Jlik1 <- Simplify( paste("(", unlist(Jllike$loglike), ")", collapse="+")) )  
  p <- length(Jfixed)  #  dimension of fixed parameters in glme's
  q <- length(Jraneff) # diemnsion of random effects
  GlmeParValue <- as.list(unlist(Jllike$str_val))   # starting values for longitudinal models

    
  ########################## settings for Survival modeles
  (Sllike <-  cox_loglike(survObject, Jllike))
  Jlik2 <- Sllike$loglike  # log-likelihood function
  Spar <- Sllike$par   # fixed parameters in cox model
  Sp <- length(Spar)  # dimension of fixed parameters in Cox model
  SurvParValue <- Sllike$str_val
  cox.model <- survObject$cox.model
  status <- survObject$status
  
  ########################## settings for the data
  group <- long.data[ , idVar]  # grouping variable, e.g patient ID
  uniqueID <- unique(group)   
  n <- length(uniqueID)  # sample size
  ni <- table(group)   # number of repreated measurements  
  N <- nrow(long.data)

  ########################## initinal values
  Bi0 <-  as.data.frame(matrix(0, ncol=q, nrow=n)) # initial values of random effects
  B0 <-  as.data.frame(matrix(0, ncol=q, nrow=N)) 
  names(B0)=names(Bi0)=Jraneff 
  surv.data$h0 <- surv.data$Th0 <- 0
  
  if(is.null(SIGMA)) SIGMA <- diag(1,q,q)
  invSIGMA0 <- ginv(SIGMA)
  fixedest0 <- unlist(c(GlmeParValue, SurvParValue))
  sigma0 <- Vassign(Jdisp, rep(0.5, length(Jdisp)))
  adjSurvpar0 <- fixedest0[- (1:(p+Sp-q))]*sqrt(diag(SIGMA))
  
  new_h <- estBaseHazard(surv.data, Sllike, SurvParValue, Bi0,
                         status, nblock)
  surv.data <- new_h$ndat
  
  ##########################################
  ########  begining of iteration  #########
  ##########################################
  Diff=1
  m=1
  convergence=0
  
  
  while(Diff>iterVal & m<iterMax){

    # estimate random effects
    nB <- nBi <- c()      
    
    for(i in 1:n){
      subdat <- subset(long.data, long.data[, idVar]==uniqueID[i])
      subsurv.dat <- subset(surv.data, surv.data[, idVar]==uniqueID[i])
      
      bi <- estRaneff(RespLog=list(Jlik1, Jlik2), 
                    raneff=Jraneff, 
                    long.data=subdat, surv.dat=subsurv.dat, 
                    invSIGMA=invSIGMA0, sigma=sigma0, 
                    ParVal=fixedest0,
                    Silent=Silent)    
      nBi <- rbind(nBi, bi)
      nB <- rbind(nB, matrix(rep(bi, ni[i]), ncol=q, byrow=T))           
      if(Silent==F)  cat("i=",i, "\n")    
    }

    B <- as.data.frame(nB)
    Bi <- as.data.frame(nBi)
    row.names(Bi)=row.names(B)=c()
    names(B) <- names(Bi) <- Jraneff
    
    
    # estimate fixed parameters in longitudinal and survival model
    estResult <-  estFixeff(RespLog=list(Jlik1, Jlik2), 
                             fixed=c(Jfixed, Spar),     # pars to be estimated    
                             str_val=fixedest0,   # starting values of fixed pars
                             Dpars=Jraneff,   # pars taking derivatives to for H matrix
                             long.data, surv.data, # dataset
                             B, Bi,     # random effect values
                             invSIGMA=invSIGMA0, sigma=sigma0,    # dispersion pars
                             ParVal=NULL,   # other par values needed
                             Silent)
    
    fixedest <- estResult$gamma
    if(Silent==F) cat("fixed est=", round(fixedest, 2), '\n')
    

    # estimate dispersion parameters
    dispest <- estDisp(RespLog=list(Jlik1, Jlik2), 
                         Jdisp,      # pars to be estimated    
                         sigma0,   # starting values 
                         invSIGMA0=solve(cov(Bi)), # starting value
                         Dpars=c(Jraneff,Jfixed,Spar),   # pars taking derivatives to for H matrix
                         long.data, surv.data, # dataset
                         B, Bi,     # random effect values
                         ParVal=fixedest, # other par values, fixed
                         Silent)
    
    
    new_sigma <- as.list(dispest$sigma)
    new_invSIGMA <- dispest$invSIGMA
    
    print('cov(Bi)')
    print(cov(Bi))
    print('estimated cov(Bi)')
    print(solve(new_invSIGMA))
    
    # estimate baseline hazard function in Cox model
    new_h <- estBaseHazard(surv.data, Sllike, as.list(fixedest), Bi, 
                           status,nblock)
    surv.data <- new_h$ndat
    
    
    ####################################################    
    ################## update results ##################
    ####################################################
    adjSurvpar <- fixedest[-(1:(p+Sp-q))]*sqrt(diag(solve(new_invSIGMA)))
    
    Diff <- mean( 
      c(abs((fixedest-fixedest0)/fixedest0)[1:(p+Sp-q)],
        abs((adjSurvpar-adjSurvpar0)/adjSurvpar)))
    
    Diff0 <- mean(c(abs((fixedest-fixedest0)/fixedest0)[1:(p+Sp-q)]))
        
    if(Silent==F) cat("R.Abs.Diff:", Diff, '\n') 

    # calculate h-likelihood value
    log_val1 <- sum(with(long.data, with(B, with(as.list(fixedest), with(new_sigma,
                    eval(parse(text=Jlik1)))))))+
                sum(with(surv.data, with(Bi, with(
                  as.list(fixedest), with(new_sigma, eval(parse(text=Jlik2)))))))

    log_val2 <- sum(-diag(as.matrix(Bi)%*%new_invSIGMA%*%t(as.matrix(Bi)))/2+
                        0.5*log(det(new_invSIGMA))-q/2*log(2*pi))  
    
    hloglike_value <- log_val1+log_val2

    # calculate approximated log marginal likelihood value
    Hmats <- getHmat(RespLog=list(Jlik1, Jlik2), pars=Jraneff)
    nH1 <- evalMat(Hmats[[1]], q, long.data, 
                   par.val=c(fixedest, new_sigma), raneff.val=B)   # longitudinal part
    nH2 <- evalMat(Hmats[[2]], q, surv.data, 
                   par.val=c(fixedest, new_sigma), raneff.val=Bi)  # survival part
    nH3 <- -new_invSIGMA*n
    Hval <-  -(nH1+nH2+nH3)
    loglike_value <- hloglike_value -0.5*log(det(Hval)/2/pi)

##############
    cat("############## Iteration:", m, "###############","\n")
    cat("fixed.par:", round(fixedest, 2), "\n")
    cat("Diff=",Diff,'\n')
    cat("Diff0=",Diff0,'\n')
    cat("sigma:", round(unlist(new_sigma),2), "\n")
    cat("h-like:", hloglike_value, "\n")
    cat("loglike:", loglike_value, "\n")
    cat("###################################","\n")

    fixedest0 <- fixedest
    invSIGMA0 <- new_invSIGMA
    sigma0 <- new_sigma
    adjSurvpar0 <- adjSurvpar

    m=m+1  
  }
  
  
    if(Diff>iterVal & m>=iterMax){
       warning(paste("Iteration limit reached without covergence. Diff=", Diff))
       convergence=1
     }

    # estimate sd's of parameter estimates  
      finalHmat <- getHmat(RespLog=list(Jlik1, Jlik2), 
                           pars=c(Jraneff, Jfixed, Spar))
      mat1 <- evalMat(finalHmat$negH_long, q=p+q+Sp, data = long.data, 
              par.val =c(sigma0, fixedest0), raneff.val = B) 
      mat2 <- evalMat(finalHmat$negH_surv, q=p+q+Sp, data = surv.data, 
                      par.val = c(sigma0, fixedest0), 
                      raneff.val = Bi)  
      mat3 <- bdiag(-invSIGMA0*n, diag(0, p+Sp))
      Hval <-  as.matrix(-(mat1+mat2+mat3))
      covMat <- ginv(Hval)
      se2 <- diag(covMat)[-(1:q)]
      se2[(p+Sp+1):(p+Sp+q)] <- se2[(p+Sp+1):(p+Sp+q)]*diag(solve(invSIGMA0))
      
        return(list(fixedest=fixedest0, 
                    varest=se2,
                    raneff=Bi, 
                    B=B,
                    invSIGMA=invSIGMA0, sigma=sigma0,
                #    surv.data=surv.data,
                    convergence=convergence,
                    loglike_value=loglike_value,
                    hloglike_value=hloglike_value))
      
  }


