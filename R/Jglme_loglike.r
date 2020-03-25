
# This function returns joint log-likelihood function of 
# multiple longitudinal variables modeled by GLME models 
# (including LME model as it's a special case of GLME).
# It also incorporates the log-likelihood of cneosring model.

# condition on: glme_loglike() 

Jglme_loglike <- function(glmeObject, std=T){
  
  k <- length(glmeObject)  # number of glme models
  lik <- fixed <- raneff <- rvX <- resp <- rvZ <- linear_pred <- 
    str_val <- sigma <- lower <- upper <- 
    disp <- as.list(rep(NA,k))
   
  for(i in 1:k){
      glmeObject_i <- glmeObject[[i]]
      glmmReturn <- glme_loglike(glmeObject_i, std)
      
      if(is.null(glmeObject_i$CenObject) | 
         !is.list(glmeObject_i$CenObject)){  # NO CENSOR
        
        fixed[[i]] <- glmmReturn$fixed  # fixed parameters
        raneff[[i]] <- glmmReturn$raneff # random effects
        rvX[[i]] <- glmmReturn$rvX # covariates for fixed pars
        rvZ[[i]] <- glmmReturn$rvZ  # covariates for random effects
        resp[[i]] <- glmmReturn$resp  # responses
        linear_pred[[i]] <- glmmReturn$linear_pred
        str_val[[i]] <- glmeObject_i$str_val
        names(str_val[[i]]) <- fixed[[i]]
        sigma[[i]] <- glmmReturn$sigma
        lower[[i]] <- glmeObject_i$lower
        upper[[i]] <- glmeObject_i$upper
        disp[[i]] <- glmeObject_i$disp
        
        if(is.null(glmeObject_i$CenObject)){
          lik[[i]] <- glmmReturn$loglike  # log-likelihood 
        } else {
          
          Cresp <- glmeObject_i$CenObject[2]
          delim_val <- glmeObject_i$CenObject[3]
          CstdNmu <- paste("(",delim_val,"-(",
                           glmmReturn$linear_pred,"))/",
                           glmmReturn$sigma)
          leftint <- paste("1/(1+exp(-1.702*", CstdNmu,"))")
          
          if(glmeObject_i$CenObject[1]=="tobit"){
            lik[[i]] <- paste( "(", glmmReturn$loglike, 
                               ")*(1-", Cresp, ") +log(",
                               leftint, 
                               ")*",Cresp )
          }
          
          if(glmeObject_i$CenObject[1]=="truncated"){
            lik[[i]] <- paste("(", glmmReturn$loglike, 
                                       "-log(1-", leftint,
                                       "))*(1-", Cresp, ")")
          }
          
          
        }

      } else if(is.list(glmeObject_i$CenObject)){ 
        CenObject <- glmeObject_i$CenObject
        Creturn <- glme_loglike(CenObject, std)
        Cresp <- Creturn$resp
        fixed[[i]] <- c(glmmReturn$fixed, Creturn$fixed) 
        raneff[[i]] <- c(glmmReturn$raneff, Creturn$raneff)
        rvZ[[i]] <- c(glmmReturn$rvZ, Creturn$rvZ)
        rvX[[i]] <- c(glmmReturn$rvX, Creturn$rvX)   
        linear_pred[[i]] <- c(glmmReturn$linear_pred, Creturn$linear_pred)
        resp[[i]] <- c(glmmReturn$resp, Creturn$resp)
        str_val[[i]] <- c(glmeObject_i$str_val, glmeObject_i$CenObject$str_val)
        names(str_val[[i]]) <- fixed[[i]]
        sigma[[i]] <- glmmReturn$sigma
        lower[[i]] <- c(glmeObject_i$lower, glmeObject_i$CenObject$lower)
        upper[[i]] <- c(glmeObject_i$upper, glmeObject_i$CenObject$upper)
        disp[[i]] <- c(glmeObject_i$disp, glmeObject_i$CenObject$disp)
        # Assume that there is only one regime in the censored data 
        # and the uncensored data follows a normal distribution
        if(CenObject$Cregime==1 & CenObject$truncated==F){ 
            lik[[i]] <- paste( "(", glmmReturn$loglike, 
                              ")*(1-", Cresp, ") +(",
                              Creturn$linear_pred, 
                              ")*",Cresp, "-log(1+exp(", 
                              Creturn$linear_pred, "))")
      
        } 
        
        # Assume that there is only one regime in the censored data 
        # and the uncensored data follows a truncated normal distribution
        if(CenObject$Cregime==1 & CenObject$truncated==T){
          CstdNmu <- paste("(",CenObject$delim_val,"-(",
                           glmmReturn$linear_pred,"))/",
                           glmmReturn$sigma)
          
          leftint <- paste(
            "1/(1+exp(-1.702*", CstdNmu,"))")
          
          lik[[i]] <- paste("(", glmmReturn$loglike, 
                              "-log(1-", leftint,
                            "))*(1-", Cresp, ")+(",
                              Creturn$linear_pred, 
                              ")*",Cresp,"-log(1+exp(", 
                              Creturn$linear_pred, "))")
        }
        
        # Assume that there are two regimes in the censored data, one from 
        # normal distribution and one from point mass, and the uncensored 
        # data follows a normal distribution
        
        if(CenObject$Cregime==2 ){
          CstdNmu <- paste("(",CenObject$delim_val,"-(",
                           glmmReturn$linear_pred,"))/",
                           glmmReturn$sigma)
          
          leftint <- paste(
            "1/(1+exp(-1.702*", CstdNmu,"))"
          )
          
            lik[[i]] <- paste("(",glmmReturn$loglike,"+",
                              Creturn$linear_pred,")*(1-",Cresp,")+
                              log(1+",leftint,"*exp(",Creturn$linear_pred,
                              "))*", Cresp,"-log(1+exp(", Creturn$linear_pred, "))")
        }
      }
      
  }
  
  result <- list(loglike=lik, fixed=fixed, raneff=raneff, 
                 linear_pred=linear_pred,
                 rvX=rvX, rvZ=rvZ,
                 resp=resp, str_val=str_val, sigma=sigma,
                 lower=lower, upper=upper,
                 disp=disp)
  
  return(result)
}

