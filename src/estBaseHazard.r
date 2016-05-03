
# This function estimates the baseline hazard function in the 
# Cox model as a step function.

estBaseHazard <- function(surv.data, Sllike, 
                    SurvParValue, Bi0, status,
                    nblock){
  
      Etime <- surv.data[, Sllike$resp]
      Twindow <- seq(min(Etime)-0.5, max(Etime)+0.5, length.out=nblock)
      delta1 <- Twindow[1:(nblock-1)] # lower window
      delta2 <- Twindow[2:nblock]  # upper window
      Harz <- data.frame(delta1, delta2, h0=0, Th0=0)
      linear <- parse(text=Sllike$linear_pred)
      
      
      for(i in 1:nrow(Harz)){
        
        subdat <- surv.data[Etime<delta2[i] & Etime>=delta1[i], ]
        # subdat of participants who survived up to delta1[i]
        indx <- which(Etime >= delta1[i]) 
        subdat2 <- surv.data[indx,]
        subBi0 <- Bi0[indx,]
        
        numer <- sum(subdat[, status])
        denom <-  sum( with(subdat2, with(as.list(SurvParValue), 
              with(subBi0, exp( eval(linear) )  ) ) ))    
        
        Harz$h0[i] <- numer/denom          
        Harz$Th0[i] <- sum(Harz$h0[1:i])
      }
      
      for(i in 1:nrow(surv.data)){
        (kk <- which(delta1< Etime[i] & delta2>= Etime[i]))
        (surv.data$h0[i] <- Harz$h0[kk])
        (surv.data$Th0[i] <- Harz$Th0[kk])
      }
      
      surv.data$h0[which(surv.data$h0==0)] <- 10^(-16)
      
      return(list(ndat=surv.data, H=Harz))
}