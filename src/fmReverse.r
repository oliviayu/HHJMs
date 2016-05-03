
# Given a model, this function returns the response, 
# covariates of fixed parameters, and covariates of random effects.

# Example: 
# fm <- y ~ 1+x1+x2+(1|id)
# formulaReverse(fm)

fmReverse <- function(fm){
  sp1 <- strsplit(as.character(fm), "~",  fixed=T)
  response <- sp1[[2]]   # returns response variable
  
  if( "(" %in% strsplit(as.character(sp1[[3]]), "")[[1]] ){  # if model includes random effects 
    sp2 <- as.character(strsplit(as.character(sp1[[3]]), "(",  fixed=T)[[1]][1])
    rvX <- strsplit(sp2,"+",fixed=T)[[1]] 
    rvX <- rvX[-length(rvX)]   # returns covariates of fixed parameters
    
    sp3 <- as.character(strsplit(as.character(sp1[[3]]), "(",  fixed=T)[[1]][2])
    # returns covariates of random effects
    rvZ <- strsplit(as.character(strsplit(sp3,"|",fixed=T)[[1]][1]), "+", fixed=T)[[1]]           
    
##  sp4 <- as.character(strsplit(as.character(sp1[[3]]), "|",  fixed=T)[[1]][2])
##  sp5 <- strsplit(sp4, ")", fixed=T)[[1]]
##  group <- paste(strsplit(sp5,"")[[1]][-1], collapse="")
    
  } else {
    rvX <- strsplit(as.character(sp1[[3]]), "+", fixed=T)[[1]]
    rvZ <- NULL
    # group <- NULL
  }
  
  return(list(resp=response, rvX=rvX, rvZ=rvZ)) 
}

