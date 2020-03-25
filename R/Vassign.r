
# assign values to a vector of names

Vassign <- function(name, value){
  dat <- data.frame(as.list(value))
  names(dat) <- name
  return(dat)
}