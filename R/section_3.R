#------------------------------
# Structural intensity function
# ~ from FDA 
#------------------------------


library("not")
library("pbapply")


structural.intensity <- function(X,h){
  #' Returns Gasser & Kneip function
  #'
  #' @param X, array
  #' @param h, bandwidth for 
  
  m <- ncol(X)
  Y <- pblapply(1:ncol(X), function(i) features(not(X[,i], contrast = "pcwsLinMean"))$cpt)
  Y <- Y[-which(sapply(Y, anyNA) == TRUE)]
  
  intensity <- function(t) {
    
    epan.kern <- function(x){
      u <- (x-t)/h
      if (abs(u) > 1){
        return(0)
      } else {
        return((3/4)*(1-u**2))
      }
    } 
    
    M <- sapply(Y, function(j) sum(sapply(j, function(i) epan.kern(i))))
    return((1/(m*h))*sum(M))
  }
  
  return(intensity)
  
}



#----------
# Load data 
#----------


# check if data has been collected, load / collect
if (file.exists("data/snp_covid.csv")){
  
  snp.raw.covid <- read.csv("data/snp_covid.csv", row.names = 1)
  
} else {
  
  source("R/build_datasets.R")
  
  # collect 
  date.from <- as.Date("2020/01/01")
  date.to <- date.from + 100
  snp.raw.covid <- SnP500.data(date.from,date.to)
  
  # check missing values
  if (anyNA(snp.raw.covid)){
    drop <- which(apply(snp.raw.covid,2,anyNA) == T)
    snp.raw.covid <- snp.raw.covid[,-drop]
  }
  
  # write out
  write.csv(snp.raw.covid,"data/snp_covid.csv")
}



#--------------------------
# Plot structural intensity
#--------------------------




toy.intensity <- structural.intensity(snp.raw.covid,5)
plot.ts(sapply(1:69,toy.intensity), main = "Structural intensity for changepoints")



