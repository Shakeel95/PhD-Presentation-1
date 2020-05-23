library("tseries") # fit GARCH easily
library("docstring") # EZ documentation
library("pbapply") # lovely progress bar


garch.ar <- function(x,n = 50 ,orders.list = list(c(1,1),c(0,1),c(2,1),c(2,1))) {
  #'Fit GARCH process, fint AR(\infty) expansion
  #'
  #'@param x vector
  #'@param n int, where to truncate infinite expansion
  #'@param orders.list list, orders of GARCH models to fit
  
  p <- 0; q <- 0; garch.AIC <- Inf
  phi <- 0; theta <- 0
  
  # find best fit stationar GARCH
  for (ord in orders.list) {
    
    m <- garch(x, ord, trace = FALSE)
    
    if ((sum(m$coef[-1])<1) & (garch.AIC > (2*sum(ord) - 2*m$n.likeli))) {
      
      garch.AIC <- 2*sum(ord) - 2*m$n.likeli
      
      if ((ord[1] == 0) & (ord[2] > 0 )) {
        
        q <- 0; p <- ord[2] 
        phi <- coef(m)[-1]
        theta <- 0
        
      } else if ((ord[2] == 0) & (ord[1] > 0 )){
        
        p <- 0; q <- 0; phi <- 0; theta <- 0
        
      } else {
        
        q <- ord[1]; p <- max(ord)
        theta <- (-1)*coef(m)[-1][-(1:ord[1])]
        phi <- coef(m)[-1][1:ord[1]] - theta
        
      } 
    }
  }
  
  # if best fit model is constant
  if (garch.AIC == Inf) return(integer(n))
  
  
  # else perform AR expansion
  phi <- c(phi, numeric(n))
  theta <- c(theta, numeric(n))
  pie <- c(numeric(q),1,numeric(n))
  
  for (j in 1:n) pie[j + q + 1] = -phi[j] - sum(theta * pie[ifelse(q==0,0,(q:1)) +j])
  
  return(pie[(0:n) + q + 1])
  
}



ar.coef.similarity <- function(X){
  #'AR expansion similarity 
  #'
  #'Similarity matrix for a panel of time series based on l2 distance between AR(\infty) representations
  #'
  #'@param X matrix or dataframe of dimension (T x p)
  #'@references 
  
  
  ar.dist <- function(idxs, expansions){
    # Lil' function, finds l2 distance between truncated expansions
    
    if (idxs[1] == idxs[2]) return(0.0)
    
    P <- expansions[[idxs[1]]]
    Q <- expansions[[idxs[2]]]
    
    return(sqrt(sum(P-Q)**2))
  }
  
  expansions <- pblapply(1:ncol(X),function(i) garch.ar(X[,i]))
  
  d <- pbapply(combn(ncol(X),2),2,ar.dist, expansions = expansions)
  
  attr(d, "Size") <- ncol(X)
  xnames <- colnames(X)
  if (!is.null(xnames)) {
    attr(d, "Labels") <- xnames
  }
  attr(d, "Diag") <- FALSE
  attr(d, "Upper") <- FALSE
  class(d) <- "dist"
  
  return(d)
}
