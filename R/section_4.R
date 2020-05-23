#-----------------------
# Load required packages
#-----------------------

library("not") # narrowest over threshold
library("pbapply") # show progress#
library("zoo") # utils



#---------------------------
# Load data if not loaded...
#---------------------------


if (!(("snp.raw.covid" %in% ls()) & ("snp.LR.covid" %in% ls()))) {
  
  source("build_datasets.R")
  
  date.from <- as.Date("2020/01/01");date.to <- date.from + 100
  snp.raw.covid <- SnP500.data(date.from,date.to)
  
  snp.LR.covid <- apply(log(snp.raw.covid),2,diff)
  
}



#-----------------------------------------------------------
# Show fator model is appropriate ts driven by common shocks
#-----------------------------------------------------------


# get changepoints with NOT algorithm
TT <- nrow(snp.raw.covid); n <- ncol(snp.raw.covid)
cpts.covid <- pbsapply(1:n, function(i) 
  features(not(snp.raw.covid[,i], contrast = "pcwsLinMean"))$cpt,
  simplify = "array")


X <- matrix(0,TT,n)
for (i in 1:n) X[cpts.covid[[i]],i] <- 1
colnames(X) <- names(snp.raw.covid)
rownames(X) <- rownames(snp.raw.covid)

heatmap(t(X),Colv = NA, Rowv = NA, main = "changepoint locations from NOT algorithm")



#----------------------------
# Load low volatility period
# plot...
#----------------------------


# check if data has been collected, load / collect
if (file.exists("data/snp_low_volatility.csv")){
  
  low.vol <- read.csv("data/snp_low_volatility.csv", row.names = 1)
  
} else {
  
  # collect 
  date.to <- as.Date("2005-05-01")
  date.from <- date.to-100 
  low.vol <- SnP500.data(date.from,date.to)
  
  # check missing values
  if (anyNA(low.vol)) {
    drop <- which(apply(low.vol,2,anyNA)==T)
    low.vol <- low.vol[,-drop] 
  }
  
  # write out
  write.csv(low.vol,"data/snp_covid.csv")
}


# plot log returns / show low volatility
snp.LR.lowvol <- apply(log(low.vol),2,diff)
matplot(snp.LR.lowvol, type = "l", xaxt = "n",
        main = "SnP500 (log-returns)",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.LR.lowvol), rownames(snp.LR.lowvol))


# get changepoints with NOT algorithm
TT <- nrow(low.vol); n <- ncol(low.vol)
low.vol.cpts <- pbsapply(1:n, function(i) 
  features(not(low.vol[,i], contrast = "pcwsLinMean"))$cpt,
  simplify = "array")


X <- matrix(0,TT,n)
for (i in 1:n) X[low.vol.cpts[[i]],i] <- 1
colnames(X) <- names(low.vol)
rownames(X) <- rownames(low.vol)

heatmap(t(X),Colv = NA, Rowv = NA, main = "changepoint locations from NOT algorithm")



#------------------------------
# Find groups maximising median
#------------------------------


# smooth changepoints with rolling sum
Y <- apply(X,2,function(i)rollapply(i,3,max))
Z <- Y 

groups <- list()
times <- 1:TT
considered <- c()
i <- 1

while(ncol(Z) > 1) {
  
  # base case start 
  if (length(times) == TT){
    
    tt <- which(apply(Y,1,sum) == median(apply(Y,1,sum)))[1]
    
    groups[[i]] <- which(Y[tt,] == 1)
    i <- i+1 
    times <- times[!times == tt]
    considered <- c(considered,tt)
    
    
    # else reduced search space  
  } else {
    
    if ((length(considered) > TT - 1) | (length(unlist(groups))) > n-1) break 
    
    Z <- Y[-considered,-unlist(groups)]
    tt  <- which.min(abs(apply(Z,1,sum) - median(apply(Z,1,sum))))[1]
    print(tt)
    holder <- names(which(Z[tt,] == 1))
    groups[[i]] <- which(colnames(Y) %in% holder)
    i <- i+1
    times <- times[!times == tt]
    considered <- c(considered,tt)
  }
}


# add last value
new.order <- c(unlist(groups), setdiff(1:n,unlist(groups)))

heatmap(t(X[,new.order]),Colv = NA, Rowv = NA, main = "ordered (median clusters)")

W <-  apply(X,2,function(i)rollapply(i,2,max))
heatmap(t(W[,new.order]),Colv = NA, Rowv = NA, main = "ordered (median clusters)")




#---------------------------
# Visualise cluster fomation
#---------------------------


if (!dir.exists("plots/snp_perm_gif")) dir.create("./plots/snp_perm_gif")
j <- 1

perm.order <- 1:n
exclude <- NULL


for (i in 1:n) {
  
  if (!(i %in% exclude)) {
    
    if (perm.order[i] != new.order[i]) {
      
      pos <- which(new.order == perm.order[i])
      val <- perm.order[i]
      
      perm.order[i] <- new.order[i]
      perm.order[pos] <- val
      
      exclude <- c(exclude,i,pos)
    }
   
    name <- paste("plots/snp_perm_gif/snp_perm_", j, ".png", sep = "")
    png(file = name); j <- j+1
    
    heatmap(t(X[,perm.order]),Colv = NA, Rowv = NA)
    
    dev.off()
    
  }
}

