#' This function is for joint modelling of multiple (generalized) linear mixed-effect 
#' models and a survival model (Cox PH or Weibull), via h-likelihood method. It also 
#' addresses left-truncation of a continuous longitudinal data fitted by a LME model.
#' 
JMfit_HL <- function(glmeObject, survObject, long.data, surv.data, idVar, eventTime, survFit,
                     itertol=1e-3, Ptol=1e-2, iterMax=10, Silent=T){
  
  ########################## settings for Longitudinal models
  Jllike <- Jglme_loglike(glmeObject)  
  Jraneff <- unique(na.omit(unlist(Jllike$raneff))) # random effects 
  Jfixed <- unique(unlist(Jllike$fixed))   # fixed parameters 
  Jdisp <- na.omit(unlist(Jllike$sigma))  # dispersion parameters
  Jlik1 <- paste("(", unlist(Jllike$loglike), ")", collapse="+") # log likelihood conditional on random effects
  p <- length(Jfixed)  #  dimension of fixed parameters 
  q <- length(Jraneff) # diemnsion of random effects
  GlmeParValue <- as.list(unlist(Jllike$str_val))   # starting values for fixed parameters
  
  ########################## settings for Survival model
  if(class(survFit) == "coxph"){
    Sllike <- cox_loglike(survObject) 
    nblock <- max(surv.data[, Sllike$resp])
    Haz <- basehaz(survFit)
    names(Haz) <- c("Th0", eventTime)
    Haz$h0 <- c(1e-20, pmax(1e-20, diff(basehaz(survFit)[,1],1)))
    Nsurv <- merge(surv.data, Haz, by.all=eventTime)
    Nsurv <- Nsurv[order(Nsurv[,idVar]) ,]
    surv.data <- Nsurv
  } else if(class(survFit) == "survreg" & survFit$dist == "weibull"){
    weibPar <- c( -summary(survFit)$coeff[1]/summary(survFit)$scale,
                  1/summary(survFit)$scale)
    Sllike <- cox_loglike(survObject, weibPar) 
  }
  
  Jlik2 <- Sllike$loglike  # log-likelihood function
  Spar <- Sllike$par   # fixed parameters in survival model
  Sp <- length(Spar)  # dimension of fixed parameters in survival model
  SurvParValue <- Sllike$str_val  # starting values
  
  ########################## settings for the data
  group <- long.data[ , idVar]  # grouping variable, e.g patient ID
  uniqueID <- unique(group)   
  n <- length(uniqueID)  # sample size
  ni <- table(group)   # number of repreated measurements for each subject 
  N <- nrow(long.data)
  
  ########################## initial values
  fixedest0 <- unlist(c(GlmeParValue, SurvParValue))
  invSIGMA0 <- SIGMA <- diag(1,q,q)
  sigma0 <- Vassign(Jdisp, rep(0.5, length(Jdisp)))
  
  ## adjust fixed and disp parameters
  exDisp <- unlist(Jllike$disp)
  Jfixed <- unique(c(Jfixed[!Jfixed%in%exDisp], Spar[!Spar%in%survObject$disp]))
  Jdisp <- unique(c(Jdisp, exDisp, survObject$disp))
  fixedest00 <- fixedest0[which(names(fixedest0)%in%Jfixed)]
  sigma00 <- unlist(c(sigma0, fixedest0[which(!names(fixedest0)%in%Jfixed)]))
  fixedest0 <- fixedest00
  sigma0 <- sigma00
  p <- length(fixedest0)
  Lval0 <- NULL
  
  # add bounds on the dispersion parameters
  lower <- c(0, c(unlist(Jllike$lower), survObject$lower))
  upper <- c(Inf, c(unlist(Jllike$upper), survObject$upper))
  names(lower) <- names(upper) <- Jdisp
  RespLog <- list(Jlik1, Jlik2)
  
  if(class(survFit) == "coxph"){
    fixedlower <- NULL
    # fixedlower = -Inf
  }  else if(class(survFit) == "survreg" & survFit$dist == "weibull"){
    fixedlower <- c(rep(-Inf, p-1), 0)
  }
  
  likDiff <- Diff <- 1
  convergence <- 1
  m <- 1
  
  while(likDiff > itertol & Diff > Ptol & m < iterMax){
    
    # estimate random effects
    output <- get_raneff(RespLog, long.data,surv.data, idVar, uniqueID,
                         Jraneff, invSIGMA0, sigma0, n, ni, q, N,
                         fixedest0, Silent)
    Bi <- output$Bi
    B <- output$B
    print("estimate random effects --- done.")
    
    # estimate fixed parameters in longitudinal and survival model
    estResult <- estFixeff(RespLog, 
                           fixed = Jfixed,
                           str_val = fixedest0,   # starting values of fixed pars
                           Dpars = Jraneff,   # pars taking derivatives to for H matrix
                           long.data, surv.data, # dataset
                           B, Bi,     # random effect values
                           invSIGMA = invSIGMA0, sigma = sigma0,    # dispersion pars
                           ParVal = NULL,   # other par values needed
                           lower = fixedlower,
                           Silent)
    fixedest <- estResult$gamma
    print("estimate fixed parameters --- done.")

    # estimate dispersion parameters
    dispest <- estDisp(RespLog, 
                       Jdisp,      # pars to be estimated    
                       sigma0,   # starting values 
                       invSIGMA0, # starting value
                       Dpars = c(Jraneff,Jfixed),   # pars taking derivatives to for H matrix
                       long.data, surv.data, # dataset
                       B, Bi,     # random effect values
                       ParVal = fixedest, # other par values, fixed
                       Silent, Lval0,
                       lower, upper)
    
    new_sigma <- dispest$sigma
    new_invSIGMA <- dispest$invSIGMA
    Lval0 <- dispest$Lval
    print("estimate dispersion parameters --- done.")
    
    
    # estimate baseline hazard function in Cox model
    if(is.null(survObject$distribution)){
      new_h <- estBaseHazard(surv.data, Sllike, c(as.list(fixedest),new_sigma), Bi, 
                             survObject$event, nblock)
      surv.data <- new_h$ndat
    }
    
    ####################################################    
    ################## update results ##################
    ####################################################
    # calculate approximated log marginal likelihood value
    loglike_value <- calc_llike(long.data, surv.data, RespLog,
                                Bi, B, fixedest, new_sigma, 
                                new_invSIGMA,Jraneff,q,n)
    
    if(m == 1){
      likDiff <- 1
    } else{
      likDiff <- abs(loglike_value-loglike_value0)/abs(loglike_value0)
    }
    
    # calcuate relative changes in mean parameters
    Diff <- mean(c(abs((fixedest - fixedest0)/(fixedest0 + 1e-6))))
  
    ############## print
    cat("############## Iteration:", m, "###############","\n")
    cat("fixed.par:", round(fixedest, 2), "\n")
    cat("FixedParDiff = ", Diff, '\n')
    cat("likDiff = ", likDiff, '\n')
    if(!is.null(unlist(new_sigma))){
      cat("sigma:", round(unlist(new_sigma), 2), "\n")
    }
    cat("loglike:", loglike_value, "\n")
    cat("##########################################","\n")
    fixedest0 <- fixedest
    invSIGMA0 <- new_invSIGMA
    sigma0 <- new_sigma
    loglike_value0 <- loglike_value
    m <- m+1  
  }
  
  
  ## messages about convergence success or failure
  if((likDiff > itertol  & Diff>Ptol ) ){
    warning("Iteration limit reached without covergence.")
    convergence <- 1
  }
  if(likDiff <= itertol & likDiff >= 0){
    message("Successful convergence. Iteration stops because likDiff <= itertol.")
    convergence <- 0
  }
  if(Diff <= Ptol){
    message("Successful convergence. Iteration stops because FixedParDiff <= Ptol.")
    convergence <- 0
  }
  
  # estimate sd's of parameter estimates  
  sd <- get_sd(long.data, surv.data, Bi, B,
               fixedest0, sigma0, invSIGMA0,
               RespLog, Jfixed, Jraneff, p, q, n)
  
  
  return(list(fixedest = fixedest0, 
              fixedsd = sd,
              Bi = Bi, 
              B = B,
              covBi = solve(invSIGMA0), 
              sigma = sigma0,
              convergence = convergence,
              loglike_value = loglike_value,
              long.data = long.data,
              surv.data = surv.data,
              RespLog = RespLog,
              idVar = idVar, uniqueID = uniqueID,
              Jraneff = Jraneff))
}


