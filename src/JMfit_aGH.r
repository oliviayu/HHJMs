## newly updated on Feb 9, 2017
# This function is for joint modelling of multiple 
# (generalized) linear mixed-effect models and a Cox model, 
# via adaptive GH method. It also addresses left-censoring 
# problem in a continuous longitudinal data, e.g. NAb.

JMfit_aGH <- function(
  #### arguments, must specified by users  ####
  glmeObject, 
  survObject,
  long.data, surv.data,  
  idVar, eventTime, survFit,

  # arguments by defualt
  itertol=1e-3, Ptol=1e-2, epsilon=10^{-6},
  iterMax=10,  ghsize=4, 
  Silent=T
){

  ########################## settings for Longitudinal modeles
  Jllike <- Jglme_loglike(glmeObject)  
  Jraneff <- unique(na.omit(unlist(Jllike$raneff))) # random effects in glme's
  Jfixed <- unique(unlist(Jllike$fixed))   # fixed parameters in glme's
  Jdisp <- na.omit(unlist(Jllike$sigma))  # dispersion parameters
  Jlik1 <- paste("(", unlist(Jllike$loglike), ")", collapse="+") 
  p <- length(Jfixed)  #  dimension of fixed parameters in glme's
  q <- length(Jraneff) # diemnsion of random effects
  GlmeParValue <- as.list(unlist(Jllike$str_val))   # starting values for longitudinal models
  
  ########################## settings for Survival modeles
  if(class(survFit)=="coxph"){
    nblock=max(surv.data[, Sllike$resp])
    Haz = basehaz(survFit)
    names(Haz) <- c("Th0", eventTime)
    Haz$h0 <- c(1e-20, pmax(1e-20, diff(basehaz(survFit)[,1],1)))
    Nsurv <- merge(surv.data, Haz, by.all=eventTime)
    Nsurv <- Nsurv[order(Nsurv[,idVar]) ,]
    surv.data=Nsurv
    weibPar = NULL
  } else if(class(survFit)=="survreg" & survFit$dist=="weibull"){
    weibPar = c( -summary(survFit)$coeff[1]/summary(survFit)$scale,
                 1/summary(survFit)$scale)
  } 
  
  Sllike <-  cox_loglike(survObject, weibPar)
  Jlik2 <- Sllike$loglike  # log-likelihood function
  Spar <- Sllike$par   # fixed parameters in cox model
  Sp <- length(Spar)  # dimension of fixed parameters in Cox model
  SurvParValue <- Sllike$str_val
  
  ########################## settings for the data
  group <- long.data[ , idVar]  # grouping variable, e.g patient ID
  uniqueID <- unique(group)   
  n <- length(uniqueID)  # sample size
  ni <- table(group)   # number of repreated measurements  
  N <- nrow(long.data)
  
  ########################## initinal values
  fixedest0 <- unlist(c(GlmeParValue, SurvParValue))
  sigma0 <- Vassign(Jdisp, rep(0.5, length(Jdisp)))
  invSIGMA0 <- SIGMA <- diag(1,q,q)
  
  ## adjust fixed and disp parameters
  exDisp <- unlist(Jllike$disp)
  Jfixed <- unique(c(Jfixed[!Jfixed%in%exDisp], Spar[!Spar%in%survObject$disp]))
  Jdisp <- unique(c(Jdisp, exDisp, survObject$disp))
  fixedest00 <- fixedest0[which(names(fixedest0)%in%Jfixed)]
  sigma00 <- unlist(c(sigma0, fixedest0[which(!names(fixedest0)%in%Jfixed)]))
  fixedest0=fixedest00; sigma0=sigma00
  p <- length(fixedest0)
  Lval0 = NULL
  
  
  # add restrictions on the dispersion parameters
  lower=c(0, c(unlist(Jllike$lower), survObject$lower))
  upper=c(Inf, c(unlist(Jllike$upper), survObject$upper))
  names(lower) <- names(upper) <- Jdisp
  RespLog=list(Jlik1, Jlik2)
  if(survObject$distribution=="weibull"){
    fixedlower=c(rep(-Inf, p-1), 0)
  } else {
    fixedlower=rep(-Inf, p)
  } 
  
  Alower = c(rep(-Inf, p), lower) 
  Aupper = c(rep(Inf,p), upper)


  ##########################################
  ########  begining of iteration  #########
  ##########################################
  likDiff=Diff=1
  convergence=1
  m=1
  
  fixedest0 = c(fixedest0,sigma0)
  SIGMA0 =solve(invSIGMA0)
  
  GHzsamp = mgauss.hermite(n=ghsize, mu=rep(0,q), sigma=diag(1/2, q))
  
  while(likDiff > itertol & Diff>Ptol & m<iterMax){

    # estimate bi's
    output <- get_raneff(RespLog, long.data,surv.data, idVar, uniqueID,
                         Jraneff, invSIGMA0, sigma0, n, ni, q, N,
                         fixedest0, Silent)
    Bi <- output$Bi
    B <- output$B
    # print("estimate random effects --- done.")
    
    # get cov(bi) for each subject i
    idsigma = cholHi_adaptGH(RespLog=list(Jlik1, Jlik2), 
                           Jraneff, long.data, surv.data,
                           Bi,  invSIGMA0, sigma=NULL,  
                           ParVal=fixedest0,
                           idVar, uniqueID, Silent)
    
    # generate GH samples by subject
    GHsample <- as.list(rep(NA,n))
    for(i in 1:n){
      GHsample[[i]] = mgauss.hermite(n=ghsize, mu=as.numeric(Bi[i,]), sigma=idsigma[[i]])
    }

    ## estimate all the parameters
   est_result = estFixeff_adptGH(RespLog=list(Jlik1, Jlik2), 
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
                                 lower=Alower,
                                 upper = Aupper,
                                 Silent)
    
    fixedest = est_result$gamma[1:length(fixedest0)]
    new_invSIGMA = est_result$invSIGMA
    ####################################################    
    ################## update results ##################
    ####################################################
    # calculate approximated log marginal likelihood value
    loglike_value <- calc_llike(long.data, surv.data, RespLog,
                                Bi=Bi, B=B, fixedest, new_sigma=NULL, 
                                new_invSIGMA, Jraneff,q,n)
    
    if(m==1){
      likDiff = 1
    } else{
      likDiff <- abs(loglike_value-loglike_value0)/abs(loglike_value0)
    }
    
    # calcuate relative changes in mean parameters
    Diff <- mean(c(abs((fixedest-fixedest0)/fixedest0)))
    
    
    ############## print
    cat("############## Iteration:", m, "###############","\n")
    cat("fixed.par:", round(fixedest, 2), "\n")
    cat("FixedParDiff=", Diff, '\n')
    cat("likDiff=", likDiff,'\n')
    cat("loglike:", loglike_value, "\n")
    cat("###################################","\n")
    
    fixedest0 <- fixedest
    invSIGMA0 <- new_invSIGMA
    loglike_value0 <- loglike_value
    
    SIGMA0 = solve(invSIGMA0)
    SIGMA0 = diag(1/sqrt(diag(SIGMA0)))%*%SIGMA0%*%diag(1/sqrt(diag(SIGMA0)))
    invSIGMA0 = solve(SIGMA0)
    
    m=m+1  
  }
  
  
  ## messages about convergence success or failure
  if((likDiff > itertol  & Diff>Ptol ) ){
    warning("Iteration limit reached without covergence.")
    convergence=1
  }
  
  if(likDiff <= itertol & likDiff >= 0){
    message("Successful convergence. Iteration stops because likDiff <= itertol.")
    convergence=0
  }
  
  if(Diff<=Ptol){
    message("Successful convergence. Iteration stops because FixedParDiff <= Ptol.")
    convergence=0
  }
  
  
  # estimate sd's of parameter estimates  
  sd = try(get_aGH_sd2(long.data, surv.data,
                          fixedest0=c(fixedest0),
                          RespLog, p, q, GHzsamp, 
                          GHsample,
                          Dpars=Jraneff, ghsize,
                          uniqueID, idVar, 
                          invSIGMA=invSIGMA0, epsilon,
                          otherval=NULL), silent = T)
  
  return(list(fixedest=fixedest0, 
              fixedsd=sd,
              Bi=Bi, 
              B=B,
              covBi=solve(invSIGMA0), 
              sigma=sigma0,
              convergence=convergence,
              loglike_value=loglike_value,
              long.data=long.data,
              surv.data=surv.data,
              RespLog=RespLog,
              idVar=idVar, uniqueID=uniqueID,
              Jraneff=Jraneff))
}


