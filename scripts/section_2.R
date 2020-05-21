#--------------------------------
# Make pcws Lin contrast function
#--------------------------------


contrast <- function(b, s, e){
 
 l <- e-s
 A <- sqrt(6/(l*(l**2-1)*(1+(e-b+1)*(b-s)+(e-b)*(b-s-1))))
 B <- sqrt(((e-b+1)*(e-b))/((b-s-1)*(b-s)))
 
 contrast <- numeric(l)
 
 contrast[(s+1):b] <- (A*B)*((3*(b-s)+(e-b)-1)*c((s+1):b) - (b*(e-s-1)+2*(s+1)*(b-s)))
 
 contrast[(b+1):e] <- (-1)*(A/B)*((3*(e-b)+(b-s)+1)*c((b+1):e) - (b*(e-s-1)+2*e*(e-b+1)))
 
 return(contrast)
 
}



#----------------------------
# LOCALISED GREEDY SEARCH 
#
# set up signal + noise model
#----------------------------


# mk model 
set.seed(123)
n <- 150  
signal.clean <- c(0.5*(1:50),rep(25,50),25-0.5*(1:50))
signal.noisy <- signal.clean + rnorm(n, sd = 2)


# catch CUSUM 
CUSUM.clean <- numeric(n)
CUSUM.noisy <- numeric(n)



#---------------------------
# LOCALISED GREEDY SEARCH
#
# Loop over frames, make gif
#---------------------------


if (!dir.exists("plots/pscwLin_CUSUM_gif")) dir.create("./plots/pscwLin_CUSUM_gif")

for (i in 2:(n-2)){
  
  # set plotting options, save png 
  name <- paste("plots/pscwLin_CUSUM_gif/CUSUM_", i, ".png", sep = "")
  png(file = name)
  layout(matrix(c(1,1,1,3,1,1,1,3,2,2,2,3), nrow = 3, byrow = T))
  
  # plot signal + noise 
  plot(1:150, signal.clean,
       type = "l", col = "red", lwd = 2, lty = 2, 
       ylim = c(0,30),
       xlab = "Time",
       ylab = "Value",
       main = "Signal + Noise"
  )
  lines(signal.noisy, col = "grey", lwd = 2, main = 1)
  abline(v = i, col = "blue")
  
  # get CUSUM
  x <- contrast(i,0,n)
  CUSUM.noisy[i] <- abs(sum(x*signal.noisy))
  CUSUM.clean[i] <- abs(sum(x*signal.clean))
  
  # plot CUSUM
  plot(1:n, CUSUM.clean, type = "l", col = "red", lwd = 2,
       ylim = c(0,100),
       xlab = "Time",
       ylab = "CUSUM",
       main = "CUSUM Value"
  )
  lines(1:n, CUSUM.noisy, type = "l", col = "grey", lwd = 2)
  abline(v = i, col = "blue")
  
  # plot contrast function 
  plot(1:n, x, col = "blue", type = "l", 
       xlab = "Time",
       ylab = "Contrast", 
       main = "Contrast"
  )
  
  dev.off()
}



#-------------------------
# ILLUSTRATING TGUW
#
# Set up model, transform
#-------------------------

# mk model
signal.clean <- c(0.5*(1:25), 12.5 - 0.5*(1:25))
n <- length(signal.clean)
signal.noisy <- signal.clean + rnorm(n)


# 
signal.transformed <- signal.noisy
detail.coef <- numeric(n)


# transform
tguwfit <- TGUW(signal.noisy, p = 0)



#---------------------------
# LOCALISED GREEDY SEARCH
#
# Loop over frames, make gif
#---------------------------


if (!dir.exists("plots/TGUW_gif")) dir.create("./plots/TGUW_gif")

j <- 1 

for (i in 1:(n-2)) {
  
  
  # set plotting options, save png
  name <- paste("plots/TGUW_gif/TGUW_", j, ".png", sep = "")
  png(file = name)
  layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,3), nrow = 3, byrow = TRUE))
  j <- j+1
  
  
  # extract itter transform
  mg <- tguwfit$merging.hist[,,i]
  s <- mg[1,1]; b <- mg[1,2]; e <- mg[1,3]
  l <- seq(mg[3,2], mg[3,3], length.out = e-s+1)
  
  
  plot(1:n, signal.transformed,
       type = "l", col = "grey", lwd = 4,
       main = "TGUW Transform",
       xlab = "Time", 
       ylab = "Value"
  )
  abline(v = s, col = "blue")
  abline(v = e, col = "blue")
  lines(1:n,signal.noisy, col = "red", type = "l")
  
  
  barplot(detail.coef, main = "Detail Coefs.", horiz=TRUE)

    
  plot(1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "max", col = "white")
  text(x = 1, y = 1, paste(round(detail.coef, 1)[which.max(abs(detail.coef))]),cex = 3)
  
  dev.off()

  
  # update detail and smooth coefs
  signal.transformed[s:b] <- l[1:(b-s+1)]
  signal.transformed[(b+1):e] <- l[(b-s+2):(e-b+2)]
  detail.coef[s:e] <- mg[3,1]


  # set plotting options, save png
  name <- paste("plots/TGUW_gif/TGUW_", j, ".png", sep = "")
  png(file = name)
  layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,3), nrow = 3, byrow = TRUE))
  j <- j+1


  plot(1:n, signal.transformed,
       type = "l", col = "grey", lwd = 4,
       main = "TGUW Transform",
       xlab = "Time",
       ylab = "Value"
  )
  abline(v = s, col = "blue")
  abline(v = e, col = "blue")
  lines(1:n,signal.noisy, col = "red", type = "l", horiz=TRUE)


  barplot(detail.coef, main = "Detail Coefs.", horiz=TRUE)
  
  plot(1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "max", col = "white")
  text(x = 1, y = 1, paste(round(detail.coef, 1)[which.max(abs(detail.coef))]),cex = 3)
  
  dev.off()
  
}


