
# This function returns joint log-likelihood function of 
# multiple longitudinal variables modeled by GLME models 
# (including LME model as it's a special case of GLME).
# It also incorporates the log-likelihood of cneosring model.

# condition on: glme_loglike() 

Jglme_loglike <- function(glmeObject){
  
  k <- length(glmeObject)  # number of glme models
  lik <- fixed <- raneff <- rvX <- resp <- rvZ <- linear_pred <- 
    str_val <- sigma <- as.list(rep(NA,k))
   
  for(i in 1:k){
      glmeObject_i <- glmeObject[[i]]
      glmmReturn <- glme_loglike(glmeObject_i)
      
      if(is.null(glmeObject_i$CenObject)){  # NO CENSOR
        
        lik[[i]] <- glmmReturn$loglike  # log-likelihood 
        fixed[[i]] <- glmmReturn$fixed  # fixed parameters
        raneff[[i]] <- glmmReturn$raneff # random effects
        rvX[[i]] <- glmmReturn$rvX # covariates for fixed pars
        rvZ[[i]] <- glmmReturn$rvZ  # covariates for random effects
        resp[[i]] <- glmmReturn$resp  # responses
        linear_pred[[i]] <- glmmReturn$linear_pred
        str_val[[i]] <- glmeObject_i$str_val
        names(str_val[[i]]) <- fixed[[i]]
        sigma[[i]] <- glmmReturn$sigma

      } else { # with censoring 
        
        CenObject <- glmeObject_i$CenObject
        Creturn <- glme_loglike(CenObject)
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
            lik[[i]] <- paste("(", glmmReturn$loglike, 
                              "-log(1-leftint))*(1-", Cresp, ")+(",
                              Creturn$linear_pred, 
                              ")*",Cresp,"-log(1+exp(", 
                              Creturn$linear_pred, "))")
          } 
        
        # Assume that there are two regimes in the censored data, one from 
        # normal distribution and one from point mass, and the uncensored 
        # data follows a normal distribution
        
        if(CenObject$Cregime==2){
            lik[[i]] <- paste("log(",Cresp, "*leftint+(1-",Cresp,")*exp(",
                              glmmReturn$loglike,")/(1-leftint)+ exp(", 
                              Creturn$linear_pred, 
                              "))-log(1+exp(", 
                              Creturn$linear_pred, "))")
        }
      
      }
      
  }
  
  result <- list(loglike=lik, fixed=fixed, raneff=raneff, 
                 linear_pred=linear_pred,
                 rvX=rvX, rvZ=rvZ,
                 resp=resp, str_val=str_val, sigma=sigma)
  
  return(result)
}

