#R Functions

#After finishing the first version of the collection of functions,
#pack the functions into one R package for more convenient usage. 

# Store and organize all the functions we need

## Interval-Data Class
#Define a class called "intvl" (short for Interval data)

#test the existing interval function in R 
#library(lubridate) #used to deal with date
#library(survival) #in survival analysis 

#load the required packages

require(ggplot2)
require(reshape2)

#So far I want to design my own class of intvl data for general purpose. 
#it also allows me to impose more flexible design. 


#Create my own interval-data class: 
  
#- So far I may directly treat it like a data frame;
#- The main reason is that I need the two end (left and right of the data). 
#- Data frame preserves the matrix operations and also allow more direct visualization under ggplot2. 


#create my own interval-data class
#intvl: interval
#create an interval object
#so far I may directly treat it like a data frame
#data frame preserves the matrix operations and also allow more

intvl.dataframe <- function(left, right, intvl.lr.check = FALSE){
  #intvl.check: whether to check if left <= right
  #in our setup we do not require left <= right
  #the check (length(left) != length(right)) has already been done by data.frame and also by the left>right 
  if (length(left) != length(right)) {
    stop("Left and right vector differ in length.")
  }
  if (intvl.lr.check){
    if (any(left>right)){
      message("(This check is optional to users.) \n","Some left points are greather than right points.")
    }
    re <- data.frame(left=left, right=right)
  }
  re <- data.frame(left=left, right=right)
  return(re)
}


#Since the data frame is not so flexible as matrix, 
#if some operations are not allow for data frame, 
#I may need to switch to matrix. 
#To accomodate this uncertainty, 
#I will stay with matrix form. 

#avoid pre-optimization

#Operations I may need: 
  
#  - transpose of a matrix; 

#- inverse of a matrix. 

#- multidimensional matrix: array. Matrix can better interact with general arrays (since matrix is a 2-dim array). 

#Notes: 
#  - Although it may look not so direct as data.frame, but it is simple. 
#  - e.g. the classical simulation will create a 1-dim arrary, here we create a 2-dim arrary. 


#matrix form
intvl <- function(left, right, lr.check = FALSE){
  #matrix form
  #lr.check: whether to check if left <= right
  if (length(left) != length(right)) {
    stop("Left and right vector differ in length.")
  }
  if (lr.check){
    if (any(left>right)){
      # warning("(This check is optional to users.) \n", "There exist intervals with its left end greater than the right end.")
      warning("(This check is optional to users.) \n", "There exist intervals with: left end > right end.")
    }
    mat <- cbind(left, right)
  }
  mat <- cbind(left, right)
  colnames(mat) <- c("left", "right")
  return(mat)
}

as.intvl <- function(vec){
  #vec is a bivariate vector
  #force a vector becomes an interval with left<right
  re <- c(min(vec), max(vec))
  names(re) <- c("left", "right")
  return(re)
}

mat2intvl <- function(x, lr.check = FALSE){
  #x is a matrix
  if (lr.check){
    if (any(left>right)){
      #warning("(This check is optional to users.) \n", "There exist intervals with its left end greater than the right end.")
      warning("(This check is optional to users.) \n", "There exist intervals with: left end > right end.")
    }
    x.new <- x
    colnames(x.new) <- c("left", "right")
  }
  x.new <- x
  colnames(x.new) <- c("left", "right")
  return(x.new)
}



#test the function
#x <- intvl(1:3, 2:1, intvl.lr.check = FALSE) #return error
#x <- intvl(1:3, 3:1) #no error
#x <- intvl(1:3, 3:1, lr.check = TRUE)  #warning required by user


## Operate on Intvl
#Check the interval arithmetic 

###Center and Range


center.intvl <- function(intvl) {(intvl[1]+intvl[2])/2}
range.intvl <- function(intvl) {intvl[2]-intvl[1]}
radius.intvl <- function(intvl) {(intvl[2]-intvl[1])/2}

#range and radius are allowed to be negative

###Interval Conjugate

intvlconj <- function(intvl) c(intvl[2], intvl[1])


###Summation 


#preparation
cummean <- function(x) cumsum(x)/seq_along(x)
#for checking LLN
cumsqrtn <- function(x) cumsum(x)/sqrt(seq_along(x))
#for checking CLT


#We can directly use "+" and "-" for "intvl" objects. 

# check
# x <- intvl(1:3, 3:1)
# y <- intvl(2:5, 5:2)
# z <- intvl(2:4, 4:2)
# x+y #return error
# x+z #no error
# x-z



## Products

#We need to consider the sequential independence here. 

#- Product rules for intervals; 
#- Product rules under the expectation operation.




###Apply Functions

#general function
#we need to do the optimal
#apply to the interval data
#re <- optimize(function(x) x^2, interval = c(-1,2), maximum = TRUE)
#re$objective

apply.intvl0 <- function(dat, FUN){
  #directly apply the function 
  re.intvl <- apply(dat, 1, FUN) 
  t(re.intvl)
}

apply.intvl <- function(dat, FUN, ind.optim = FALSE, 
                        tol = .Machine$double.eps^0.25){
  #do we need to do the optimization to produce the results
  #tol: tolerance for the optimization 
  #it should also work for dat with only one row
  if (ind.optim) {
    FUN.opt <- function(x){
      opt.re.max <- optimize(FUN, interval = as.intvl(x), 
                             maximum = TRUE, tol = tol)
      opt.re.min <- optimize(FUN, interval = as.intvl(x), 
                             maximum = FALSE, tol = tol)
      min.re <- opt.re.min$objective
      max.re <- opt.re.max$objective
      c(min.re, max.re)
    }
    re.intvl <- apply(dat, 1, FUN.opt)
    re <- t(re.intvl)
    
  } else {
    re.intvl <- apply(dat, 1, function(x) as.intvl(FUN(x))) 
    #vector will be treated as column by default. 
    re <- t(re.intvl)
  }
  return(re)
}

#check 
#z1 <- intvl(1,2)
#as.intvl(varphi(z1))
#re.1 <- apply(z1, 1, function(x) c(x[1]^2, x[2]^2))
#re.2 <- t(re.1)
#re.2

cumsum.intvl <- function(dat){
  #dat is interval data 
  #n <- nrow(dat)
  #re <- sapply(seq_len(n), function(k) apply())
  apply(dat, 2, cumsum)
}

#re <- cumsum.intvl(w2.semiGnorm)

#re1 <- cumsum(w2.semiGnorm[,1])
#re <- apply.intvl(w.semiGnorm, function(x) x^2)

#I may put it into ggplot context

#create an interval-class object
#so that it allows me to do calculation and visualization accordingly 

#write my own function at this stage 

#check the existing R package to see whether they are useful



# apply.intvl <- function(dat, FUN, ind.optim = FALSE, 
#                        tol = .Machine$double.eps^0.25){
#   #do we need to consider the optimal function
#   if (ind.optim) {
#     FUN.opt <- function(x){
#       opt.re.max <- optimize(FUN, interval = as.intvl(x), 
#                              maximum = TRUE)
#       opt.re.min <- optimize(FUN, interval = as.intvl(x), 
#                              maximum = FALSE)
#       min.re <- opt.re.min$objective
#       max.re <- opt.re.max$objective
#       c(min.re, max.re)
#     }
#     re.intvl <- apply(dat, 1, FUN.opt)
#     re <- t(re.intvl)
#     
#   } else {
#     re.intvl <- apply(dat, 1, function(x) as.intvl(FUN(x))) 
#     re <- t(re.intvl)
#   }
#   return(re)
# }

#transform to a interval with left<right
# as.intvl.lr <- function(vec){
#   #vec is a bivariate vector
#   re <- c(min(vec), max(vec))
#   names(re) <- c("low", "high")
#   return(re)
# }
# 
# mat2intvl <- function(x, intvl.lr = FALSE){
#   #x is a matrix
#   #if intvl.lr = TRUE: user require the left point <= right point
#   if (intvl.lr) {
#     x.new <- apply.stl
#   }
# }



as.propintvl <- function(dat.intvl){
  #transform the interval into proper/observable interval
  apply.intvl(dat.intvl, function(x) x)
}


##Distributions

###Univariate

NA.intvl <- function(n){
  #initialization of an interval object
  #n is the number of observations
  #null.seq <- rep(NULL, n) #it will create an empty array without dimension
  NA.seq <- rep(NA,n)
  re <- intvl(NA.seq, NA.seq)
  return(re)
}

#check
#x <- NA.intvl(50)

#without considering the G-version distributions
rconstant.intvl <- function(n, c.intvl = c(0,1)){
  #for constant interval, it is allowed to be an improper one. 
  re <- intvl(left = rep(c.intvl[1], n), 
              right = rep(c.intvl[2], n))
  #re <- replicate(n, as.intvl(intvl)) 
  return(re)
}

rintvlnormal <- function(n, par){
  meanl <- par[1]
  meanr <- par[2]
  sdl <- par[3]
  sdr <- par[4]
  mean.intvl <- c(meanl, meanr)
  sd.intvl <- c(sdl, sdr)
  mu.intvl <- rconstant.intvl(n, c.intvl = mean.intvl)
  ep <- rnorm(n)
  noise.intvl <- intvl(left = sdl*ep, right = sdr*ep)
  x.intvl <- mu.intvl + noise.intvl
  x.intvl.obs <- as.propintvl(x.intvl)
  list(x.intvl = x.intvl, 
       x.intvl.obs = x.intvl.obs, 
       ep = ep)
}

#the ambiguity in the interval direction of a constant interval 
#is characerized by the maximal distribution

  
rmaximal.intvl <- function(n, mean.intvl = c(0,1)){
  #we does require the mean interval of the maximal distribution to be a proper one
  if (mean.intvl[1] > mean.intvl[2]) stop("The mean interval of maximal distribution should be proper.")
  #generate sequential independent maximal
  re <- intvl(left = rep(mean.intvl[1], n), 
              right = rep(mean.intvl[2], n))
  #re <- replicate(n, as.intvl(intvl)) 
  return(re) 
}

rsemiGnorm.intvl <- function(n, sd.intvl = c(1,2)){
  z.intvl <- rmaximal.intvl(n, mean.intvl = sd.intvl)
  y <- rnorm(n)
  y.intvl <- intvl(y,y)
  re <- NA.intvl(n)
  #key point here: Z-->Y 
  for (i in seq_len(n)){
    re[i,] <- z.intvl[i,]*y.intvl[i,]
  }
  return(re)
}

rsemiGbern.intvl <- function(n, sd.intvl = c(1,2)){
  z.intvl <- rmaximal.intvl(n, mean.intvl = sd.intvl)
  y <- 2*rbinom(n, 1, 1/2) - 1
  y.intvl <- intvl(y,y)
  re <- NA.intvl(n)
  #key point here: Z-->Y 
  for (i in seq_len(n)){
    re[i,] <- z.intvl[i,]*y.intvl[i,]
  }
  return(re)
}

rsemiGdistn.intvl <- function(n, sd.intvl = c(1,2), rdistn.std){
  #generate semi-G-version distribution based on a standardized classical distribution
  #standardized distn: mean=0, var=1. 
  #one example of rdistn.std
  #rdistn.std = function(n) rexp(n, rate = 1) - 1
  z.intvl <- rmaximal.intvl(n, mean.intvl = sd.intvl)
  y <- rdistn.std(n)
  y.intvl <- intvl(y,y)
  re <- NA.intvl(n)
  #key point here: Z-->Y 
  for (i in seq_len(n)){
    re[i,] <- z.intvl[i,]*y.intvl[i,]
  }
  return(re)
}

### Distribution Functions or DF
#cdf in theory 
pnorm.intvlr <- function(x, mu=1, sd.intvl=c(1,2)){
  sdl <- sd.intvl[1]
  sdr <- sd.intvl[2]
  v1 <- pnorm(x, mu, sdl)
  v2 <- pnorm(x, mu, sdr)
  #(x<mu)*v1*2 + (x>=mu)*(v1+v2)
  (x<mu)*v1 + (x>=mu)*v2
}

dnorm.intvlr <- function(x, mu=1, sd.intvl=c(1,2)){
  sdl <- sd.intvl[1]
  sdr <- sd.intvl[2]
  v1 <- dnorm(x, mu, sdl)
  v2 <- dnorm(x, mu, sdr)
  #(x<mu)*v1*2 + (x>=mu)*(v1+v2)
  (x<mu)*v1 + (x>=mu)*v2
}
#pnorm.intvlr(c(1,2))

pnorm.intvll <- function(x, mu=1, sd.intvl=c(1,2)){
  sdl <- sd.intvl[1]
  sdr <- sd.intvl[2]
  v1 <- pnorm(x, mu, sdr)
  v2 <- pnorm(x, mu, sdl)
  #(x<mu)*v1*2 + (x>=mu)*(v1+v2)
  (x<mu)*v1 + (x>=mu)*v2
}

dnorm.intvll <- function(x, mu=1, sd.intvl=c(1,2)){
  sdl <- sd.intvl[1]
  sdr <- sd.intvl[2]
  v1 <- dnorm(x, mu, sdr)
  v2 <- dnorm(x, mu, sdl)
  #(x<mu)*v1*2 + (x>=mu)*(v1+v2)
  (x<mu)*v1 + (x>=mu)*v2
}

### G-normal DF 

pGnorm.upper <- function(z, sdl=1, sdr=2){
  C1 <- 2*sdr/(sdl+sdr)
  C2 <- 2*sdl/(sdl+sdr)
  C1*pnorm(z/sdr)*(z<=0) + (1-C2*pnorm(-z/sdl))*(z>0)
}

dGnorm.upper <- function(z, sdl=1, sdr=2){
  C <- sqrt(2)/(sqrt(pi)*(sdl+sdr))
  C*(exp(-z^2/(2*sdr^2))*(z<=0) + exp(-z^2/(2*sdl^2))*(z>0))
}

pGnomr.lower <- function(z, sdl=1, sdr=2){
  1-pGnorm.upper(-z, sdl = sdl, sdr = sdr)
}

dGnorm.lower <- function(z, sdl=1, sdr=2){
  dGnorm.upper(-z, sdl = sdl, sdr = sdr)
}

#check <- rsemiGbern.intvl(1e2)



### Bivariate 






## Independence 


#Implement the G-EM procedure for interval data 



## Visualizations


#cite the ggplot2
#here we use ggplot2 to do the visualization
#library(ggplot2)



#A colorblind-friendly palette
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
#  scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
#  scale_colour_manual(values=cbPalette)





#get the variable name 
#ref: https://stackoverflow.com/questions/10520772/in-r-how-to-get-an-objects-name-after-it-is-sent-to-a-function

plot.intvl <- function(x.intvl, alpha = 0.3, main = NULL, 
                       linesplot = TRUE, fatten = 0.02){
  #dat should be intvl type
  #alpha: transprancy of the ribbon in between
  #add more graphic options
  if (is.null(main)) {
    main1 <- deparse(substitute(x.intvl))
    main <- paste("Interval Plot of ", main1)
  }
  dat <- x.intvl
  dat.new <- data.frame(index = seq_len(nrow(dat)), 
                        left = dat[,1], right = dat[,2], 
                        center = (dat[,2] + dat[,1])/2, 
                        range = dat[,2]-dat[,1])
  if (linesplot) {
    re <- ggplot(dat.new, aes(x = index, y = center)) +  geom_ribbon(aes(ymin=left, ymax=right), alpha = alpha, colour = "grey50") + geom_line(colour = "grey60") + geom_line(aes(y=left), colour = "blue") + geom_line(aes(y=right), colour = "red") + ylab("Interval Values") + xlab("Obs Index") + ggtitle(main)
  } else {
    #only want the interval plot
    re <- ggplot(dat.new, aes(x = index, y = center)) +  geom_pointrange(aes(ymin=left, ymax=right), fatten = fatten)  + ylab("Interval Values") + xlab("Obs Index") + ggtitle(main) + geom_line(colour = "grey60")
  }
  print(re)
}

#ggplot
#matlines

#examples
#test.seq <- rsemiGnorm.intvl(n)
#test.seq.obs <- as.propintvl(test.seq)
#plot.intvl(test.seq)
#plot.intvl(test.seq.obs, linesplot = TRUE)
#test.seq.uobs <- apply.intvl0(test.seq.obs, intvlconj)
#plot.intvl(test.seq.uobs, linesplot = TRUE)



#semi-G-normal 
#plot.intvl(test.seq, linesplot = TRUE)

#conjugate of semi-G-normal
#test.seq.conj <- apply.intvl0(test.seq, intvlconj)
#plot.intvl(test.seq.conj, linesplot = TRUE)
#switch between them 
#a mixture model 
#s.seq <- rbinom(n, 1, 1/2)
#test.seq.mix <- s.seq*test.seq + (1-s.seq)*test.seq.conj
#plot.intvl(test.seq.mix, linesplot = TRUE)



#from internet
#library(reshape2)

#version 1
hist.intvl1 <- function(dat.intvl, true.par = NULL, 
                        par.col=2, par.pch=2){
  #give a name 
  plot(dat.intvl[,1], dat.intvl[,2])
  points(true.par[1], true.par[2], col=par.col, pch=par.pch)
}

hist.intvl <- function(dat.intvl, true.par=NULL,
                       show.bias=TRUE){
  #visualization 
  #dat.intvl <- mean.estmat.intvl
  #use reshape2 library
  if (show.bias){
    dat1 <- apply.intvl(dat.intvl, function(x) x-true.par)
    dat <- melt(dat1)
    colnames(dat) <- c("Index", "EndType", "Value")
    re <- ggplot(data = dat, aes(x=Value, fill = EndType)) + geom_density(alpha=0.3) + geom_vline(xintercept = c(0,0))
    print(re)
  } else {
    dat <- melt(dat.intvl)
    colnames(dat) <- c("Index", "EndType", "Value")
    re <- ggplot(data = dat, aes(x=Value, fill = EndType)) + geom_density(alpha=0.3) + geom_vline(xintercept = true.par)
    #re2 <- ggplot(data = subset(dat, EndType=="left", Value), aes(x=Value)) + geom_density(alpha=0.3) + geom_vline(xintercept = true.par[1])
    print(re)
  }
  #geom_histogram(position = "identity", alpha = 0.4)
}
#true.par <- mu.intvl.true


## Estimations
library(Rfast)

intvlnormal.ctrrng <- function(intvl.obs, sig.level=0.95, 
                               meanCertainty=TRUE){
  rng <- apply(intvl.obs,1,range)
  rad <- apply(intvl.obs,1,radius)
  ctr <- apply(intvl.obs,1,center)
  if(meanCertainty){
    mu.est <- mean(ctr)
    a.est <- 2*sd(ctr)
    b.est <- sqrt(pi/2)*mean(rng)
    sdr.est <- (a.est+b.est)/2
    sdl.est <- (a.est-b.est)/2
    par.est <- c(mu.est, sdl.est, sdr.est)
    #get the confidence interval
    alpha <- 1-sig.level
    n <- nrow(intvl.obs)
    me <- qt(1-alpha/2, n-1)*sqrt(a.est/(2*n))
    confint.est <- mu.est + c(-1,1)*me
    names(par.est) <- c("mean.est", "sdl.est", "sdr.est")
    re <- list(par.est=par.est, 
               confint.est = confint.est)
  } else {
    normfold.est <- Rfast::foldnorm.mle(rad)
    #modify the maximization method if needed
    #Error in while (sum(abs(anew - aold)) > tol) { : 
    #missing value where TRUE/FALSE needed
    est.par <- unname(normfold.est$param)
    mu.rad <- est.par[1]
    sig.rad <- sqrt(est.par[2])
    mu.cen <- mean(ctr)
    sig.cen <- sd(ctr)
    mu.low <- mu.cen - mu.rad
    mu.up <- mu.cen + mu.rad
    sig.low <- sig.cen-sig.rad
    sig.up <- sig.cen+sig.rad
    par.est <- c(mu.low,mu.up,sig.low,sig.up)
    names(par.est) <- c("meanl.est","meanr.est", 
                        "sdl.est", "sdr.est")
    re <- par.est
  }
  return(re)
}

intvlnormal.adj <- function(intvl.obs){
  x.cen <- apply(intvl.obs, 1, center)
  x.cen.bar <- mean(x.cen)
  ind.seq <- as.numeric(x.cen-x.cen.bar<0)
  x.conj <- apply.intvl0(x.obs, intvlconj)
  #adjust.intvl <- function(a) (intvlconj(a)-a)*ind.seq + a
  x.adj <- ind.seq*x.conj + (1-ind.seq)*x.obs
  mean.est <- apply(x.adj, 2, mean)
  mean.est.cen <- mean(x.cen)
  sd.est <- apply(x.adj, 2, sd)
  par.est <- c(mean.est[1], 
               mean.est.cen, 
               mean.est[2], 
               sd.est)
  names(par.est) <- c("meanl.adj", "mean.cen.adj", "meanr.adj", "sdl.est","sdr.est")
  re <- list(x.adj = x.adj, par.est = par.est)
  re
}


## Hypothesis Tests

get.meanseq.novlp <- function(y, n.guess){
  m <- floor(length(y)/n.guess)
  y <- y[1:(m*n.guess)]  
  y.mat <- matrix(y, ncol = n.guess, byrow = TRUE)
  mean.seq <- apply(y.mat, 1, mean)
  mean.seq
}

get.meanseq.ovlp <- function(y, n.guess, step = 2){
  N <- length(y)
  n <- n.guess
  m <- floor((N-n)/step)
  ind.mat1 <- matrix(rep(1:n, m+1), ncol = n, byrow = TRUE)
  ind.mat2 <- matrix(rep((0:m)*step, n), ncol = n, byrow = FALSE)
  ind.mat <- ind.mat1 + ind.mat2
  #y.new <- y[1:(n+m*step)]
  y.mat <- matrix(NA, ncol = n, nrow = m+1)
  for (i in seq_len(m+1)){
    y.mat[i,] <- y[ind.mat[i,]]
  }
  mean.seq <- apply(y.mat, 1, mean)
  mean.seq
}

require(zoo)
#moving overlapping window
#n.guess = 3; step = 2;
get.groupfn.seq <- function(y, n.guess, novlp.ind = TRUE, 
                            step = NULL, groupfn = mean){
  y.ts <- zoo(y)
  if(novlp.ind){
    step = n.guess
  }
  re.temp <- rollapply(y.ts, width = n.guess, by = step, 
                       FUN = groupfn, align="left")
  re <- as.numeric(re.temp)
  re
}
#non-ovlp, let n.guess=step
#y.seq <- 1:6
#get.groupfn.seq(y = y.seq, n.guess = 3, step = 2, groupfn = mean)
#get.groupfn.seq(y = y.seq, n.guess = 3, step = 3, groupfn = mean)


