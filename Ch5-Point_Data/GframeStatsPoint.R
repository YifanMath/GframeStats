#organize our existing codes in the G-frame Stats

library(ggplot2)
library(MASS)
library(reshape2)
library(quantmod)

#user specify 
#this kind of information is unknown to the data analysts
#block length
#switching rule
#rmaximal
rmaximal <- function(n, min = 1, max = 2, gr.len.mean = 1e2, 
                     #len.gen.f = function(len) rpois(1, len), 
                     len.gen.f = function(len) rnorm(1, mean = len, sd = len/20), 
                     level.gen.f = function(k){
                       #return level in [0,1]
                       #alpha <- runif(1, 0,50)
                       #beta <- runif(1, 0,50)
                       #alpha <- beta <- 0.3
                       #lev <- rbeta(k, alpha, beta) 
                       lev <- sample(c(0,0.3,1), k, prob = c(3,2,1)/6)
                       lev
                     }) {
  a <- min; b <- max
  result <- numeric(n)
  count <- 0
  while (count < n) {
    #random decide the length
    #lambda <- sample(1:1e3, 1, prob = rep(1,1e3)/1e3)
    #len <- rpois(1, lambda)
    len <- len.gen.f(gr.len.mean) #more suitable for max-mean estimation
    lev <- level.gen.f(1)
    c <- lev*(b-a) + a
    #dirac distributions (constant)
    newvals <- numeric(len)+c
    result[count + seq_along(newvals)] <- newvals 
    #make some space for the newvals
    count <- count + length(newvals)
  }
  result[seq_len(n)]
}

#when block length = 1, it becomes a mixture distribution (which is also a single distribution)
#when there is only one block, it becomes one single distribution
#one can also use one Markov chain to mimic this kind of pattern 
#(but in this case, the block length follows a geometric distribution characterized by the diagonal elements in the transition probability matrix)
#pack it into a R package 

#to be shared and discuss with the community (especially, Prof. Peng, Prof. Defei Zhang)

rmaximal.extrem <- function(n, min = 1, max = 2, gr.len.mean = 1e2, 
                            alpha = 0.3, beta = 0.3) {
  a <- min; b <- max
  result <- numeric(n)
  count <- 0
  while (count < n) {
    #random decide the length
    ##version 1
    #lambda <- floor(1/runif(1)) #lambda has no upper bound
    #len <- rpois(1, lambda = lambda)
    #but the length is usually too small 
    ##version 2
    #gr.len.mean
    L <- gr.len.mean*10
    lambda <- sample(1:L, 1, prob = rep(1,L)/L)
    len <- rpois(1, lambda)
    #alpha <- 0.1
    #beta <- 10
    #beta <- 0.1
    c <- rbeta(1, alpha, beta) *(b-a) + a
    #dirac distributions (constant)
    newvals <- numeric(len)+c
    result[count + seq_along(newvals)] <- newvals 
    #make some space for the newvals
    count <- count + length(newvals)
  }
  result[seq_len(n)]
}

rmaximal.extrem.eq <- function(n, min = 1, max = 2, gr.len.mean = 1e2, 
                               alpha = 0.3, beta = 0.3) {
  a <- min; b <- max
  result <- numeric(n)
  count <- 0
  while (count < n) {
    #random decide the length
    ##version 1
    #lambda <- floor(1/runif(1)) #lambda has no upper bound
    #len <- rpois(1, lambda = lambda)
    #but the length is usually too small 
    
    ##version 2
    #lambda <- sample(1:1e3, 1, prob = rep(1,1e3)/1e3)
    len <- gr.len.mean
    #alpha <- 0.1
    #beta <- 10
    #beta <- 0.1
    c <- rbeta(1, alpha, beta) *(b-a) + a
    #dirac distributions (constant)
    newvals <- numeric(len)+c
    result[count + seq_along(newvals)] <- newvals 
    #make some space for the newvals
    count <- count + length(newvals)
  }
  result[seq_len(n)]
}

rmaximal.unif.eq <- function(n, min = 1, max = 2, gr.len.mean = 1e2) {
  a <- min; b <- max
  result <- numeric(n)
  count <- 0
  while (count < n) {
    #random decide the length
    ##version 1
    #lambda <- floor(1/runif(1)) #lambda has no upper bound
    #len <- rpois(1, lambda = lambda)
    #but the length is usually too small 
    ##version 2
    #lambda <- sample(1:1e3, 1, prob = rep(1,1e3)/1e3)
    len <- gr.len.mean
    #alpha <- 0.1
    #beta <- 10
    #beta <- 0.1
    c <- runif(1, min = a, max = b)
    #dirac distributions (constant)
    newvals <- numeric(len)+c
    result[count + seq_along(newvals)] <- newvals 
    #make some space for the newvals
    count <- count + length(newvals)
  }
  result[seq_len(n)]
}

rmaximal.incre.eq <- function(n, min = 1, max = 2, gr.len.mean = 1e2) {
  a <- min; b <- max
  m <- ceiling(n/gr.len.mean)
  s.seq <- seq(a,b,length.out = m)
  result <- rep(s.seq, each = gr.len.mean)
  result[seq_len(n)]
}

rmaximal.decre.eq <- function(n, min = 1, max = 2, gr.len.mean = 1e2) {
  a <- min; b <- max
  m <- ceiling(n/gr.len.mean)
  s.seq <- seq(b,a,length.out = m)
  result <- rep(s.seq, each = gr.len.mean)
  result[seq_len(n)]
}

#plot(rmaximal.unif.eq(1e3), type = "l")

rmaximal1 <- function(n, min = 1, max = 2) {
  a <- min; b <- max
  result <- numeric(n)
  count <- 0
  while (count < n) {
    #random decide the length
    lambda <- sample(1:1e3, 1, prob = rep(1,1e3)/1e3)
    len <- rpois(1, lambda)
    alpha <- runif(1, 0,50)
    beta <- runif(1, 0,50)
    c <- rbeta(1, alpha, beta) *(b-a) + a
    #dirac distributions (constant)
    newvals <- numeric(len)+c
    result[count + seq_along(newvals)] <- newvals 
    #make some space for the newvals
    count <- count + length(newvals)
  }
  result[seq_len(n)]
}

#almost no blocking 
rmaximal2 <- function(n, min = 1, max = 2) {
  a <- min; b <- max
  #linear measures: beta distn
  result <- numeric(n)
  count <- 0
  while (count < n) {
    #random decide the length
    lambda <- sample(1:1e3, 1, prob = rep(1,1e3)/1e3)
    len <- rpois(1, lambda)
    alpha <- runif(1, 0,50)
    beta <- runif(1, 0,50)
    newvals <- rbeta(len, alpha, beta) *(b-a) + a
    result[count + seq_along(newvals)] <- newvals 
    #make some space for the newvals
    count <- count + length(newvals)
  }
  result[seq_len(n)]
}

rmaximal3 <- function(n, min = 1, max = 2) {
  a <- min; b <- max
  result <- numeric(n)
  count <- 0
  while (count < n) {
    #random decide the length
    lambda <- floor(1/runif(1)) 
    #lambda has no upper bound
    len <- rpois(1, lambda = lambda)
    alpha <- runif(1, 0,50)
    beta <- runif(1, 0,50)
    c <- rbeta(1, alpha, beta) *(b-a) + a
    #dirac distributions (constant)
    newvals <- numeric(len)+c
    result[count + seq_along(newvals)] <- newvals 
    #make some space for the newvals
    count <- count + length(newvals)
  }
  result[seq_len(n)]
}

rf1 <- function(n) rnorm(n, mean = b2, sd = 2)
rf2 <- function(n) rexp(n, rate = 2/(a2+b2))
rf3 <- function(n) runif(n, min = a2, max = b2)
rmaximal.list <- function(n, min = 1, max = 4, 
                          density.set = list(rf1, rf2, rf3)) {
  a <- min; b <- max
  result <- numeric(n)
  count <- 0
  L <- length(density.set)
  #size of density.set
  while (count < n) {
    ind <- sample(1:L, 1)
    #random decide the length
    lambda <- floor(1/runif(1)) 
    #lambda has no upper bound
    len <- rpois(1, lambda = lambda)
    newvals <- density.set[[ind]](len)
    newvals <- newvals[a < newvals & newvals < b]
    result[count + seq_along(newvals)] <- newvals 
    #make some space for the newvals
    count <- count + length(newvals)
  }
  result[seq_len(n)]
}

#rsemiGnorm
rsemiGnorm <- function(n, sig.low=1, sig.up=2, gr.len.mean = 1e2, 
                       rmaximal.k=rmaximal){
  #choose the rmaximal function 
  #(pseudo sim of nl.iid maximal distn) 
  z.seq <- rmaximal.k(n, min=sig.low, max=sig.up, gr.len.mean = gr.len.mean)
  list(w.seq = rnorm(n) * z.seq, z.seq=z.seq)
}

#rbinom.maximal
rbinom.maximal <- function(n, size=1, rmaximal = rmaximal.extrem, p.par = c(.4,.6), varphi = function(x) 2*x-1, gr.len = 1e2){
x.seq <- numeric(n)
p.seq <- rmaximal(n, min = p.par[1], max = p.par[2], gr.len.mean = gr.len)
for (i in seq_along(p.seq)){
  p <- p.seq[i]
  x.seq[i] <- varphi(rbinom(1, size = size, prob = p))
}
list(p.seq=p.seq, x.seq=x.seq)
#x.seq
}

#from nl.CLT, generate in approximated sense
  

#rbisemiGnorm.seqind
#turn to 02-SemiGNormal-Part2-DataEx-fn

#simulate the bivariate distribution
#allow user to change the switching rule, 
#later in the LLN check, regardless of the switching rule
#we will see that it will always be sufficiently covered by the G-expectation.

##Estimation of variance uncertainty 
###Max-mean estimation

####function: max.mean (time order)

max.mean <- function(y, n.guess){
  #always choose m=N/n
  m <- floor(length(y)/n.guess)
  y <- y[1:(m*n.guess)]  
  y.mat <- matrix(y, ncol = n.guess, byrow = TRUE)
  mean.seq <- apply(y.mat, 1, mean)
  mu.low.est <- min(mean.seq)
  mu.up.est <- max(mean.seq)
  c(mu.low.est, mu.up.est)
}
 
#N is the length(y)
max.mean.new <- function(y, n.guess=floor(sqrt(N))){
  #always choose m=n 
  m <- n.guess
  y <- y[1:(m*n.guess)]  
  y.mat <- matrix(y, ncol = n.guess, byrow = TRUE)
  mean.seq <- apply(y.mat, 1, mean)
  mu.low.est <- min(mean.seq)
  mu.up.est <- max(mean.seq)
  c(mu.low.est, mu.up.est)
}
 

####function: max.mean2 (value order)

 
max.mean2 <- function(y, n.guess){
  m <- floor(length(y)/n.guess)
  y <- y[1:(m*n.guess)]
  y.mat <- matrix(y[order(y)], ncol = n.guess, byrow = TRUE)
  mean.seq <- apply(y.mat, 1, mean)
  mu.low.est <- min(mean.seq)
  mu.up.est <- max(mean.seq)
  c(mu.low.est, mu.up.est)
}

###Improvement: how to choose the group size
####function: CentralGroup.est

  
CentralGroup.est <- function(y, n.seq = seq(2,4e2,2),
                             plot.ind=TRUE, par.true=c(1,4)){
  center <- mean(y)
  N <- length(y)
  #y.ord <- y[order(y)]
  #y.ord[N/2]
  #y.ord[N/2+1]
  #median(y.ord)
  est.all.mat <- matrix(nrow = length(n.seq), ncol = 4)
  #i <- 100
  for (i in seq_along(n.seq)){
    n <- n.seq[i]
    m <- floor(N/n)
    y.mat <- matrix(y[1:(m*n)], ncol = n, byrow = TRUE)
    mean.seq <- apply(y.mat, 1, mean)
    min.seq <- apply(y.mat, 1, min)
    max.seq <- apply(y.mat, 1, max)
    mu.low.est <- min(mean.seq)
    mu.up.est <- max(mean.seq)
    est.all.mat[i,1:2] <- c(mu.low.est, mu.up.est)
    
    l.gr <- max(min.seq[min.seq<=center])
    r.gr <- min(max.seq[max.seq>=center])
    est.all.mat[i,3:4] <- c(l.gr, r.gr)
  }
  L <- est.all.mat[,1]
  R <- est.all.mat[,2]
  l <- est.all.mat[,3]
  r <- est.all.mat[,4]
  ind1 <- min(which(R-r <= 0))
  ind2 <- min(which(l-L <= 0))
  n1 <- n.seq[ind1]; n2 <- n.seq[ind2]
  a2.est <- L[ind2]
  b2.est <- R[ind1]
  if(plot.ind){
    matplot(n.seq, est.all.mat, type = "l",
            ylab = "values", xlab = "group size n", 
            main = "CentralGroup Estimation")
    abline(h=par.true, col="brown", lty=2)
    legend("top", 
           c("GroupMean.min", "GroupMean.max", 
             "CentralGroup.min", "CentralGroup.max",
             "par.true"), 
           col=c(1:4,"brown"), lty=c(1:4,2), 
           lwd=2, cex=0.4, box.lty = 2, box.col = "grey", 
           pch = 25)
  }
  return(list(est.all.mat = est.all.mat, ab.est = c(a2.est, b2.est),
              n.seq = n.seq))
}
 

####function: TimValOrd.est

  
TimValOrd.est <- function(y, step=10, plot.ind=TRUE, 
                          min.n=1, min.m=2, 
                          choprule="both.n", 
                          par.true=c(0,5),
                          start=1, len.frac=1,
                          y.scale=1){
  #time ordering 
  #c("both.n", "min.n", "max.n","fix.n")
  #min.n 
  #min.m = min number of groups
  y <- y*y.scale
  N <- length(y)
  n.seq <- seq(min.n, N/min.m, step)
  est.mat <- matrix(0, ncol=2, nrow=length(n.seq))
  for (i in seq_along(n.seq)){
    est.mat[i,] <- max.mean(y, n.seq[i])
  }
  est.mat2 <- matrix(0, ncol=2, nrow=length(n.seq))
  #i <- 2
  #?sample
  for (i in seq_along(n.seq)){
    est.mat2[i,] <- max.mean2(y, n.seq[i])
  }
  if(plot.ind){
    len <- floor(length(n.seq)*len.frac)
    ind <- start+seq_len(len)-1
    n.seq1 <- n.seq[ind]
    sum2 <- n.seq1[1]+n.seq1[length(n.seq1)]
    n.lab <- seq(0, ceiling(N/(min.m*100))*100, step*10)
    n.lab[1] <- 1
    n.lab1 <- unique(floor(N/seq(min.m, N, 1)))
    m.lab <- floor(N/n.lab1)
    matplot(n.seq[ind],cbind(est.mat[rev(ind),], est.mat2[ind,]), 
            type = "l",
            col = c(2,2,4,4),
            ylab = paste("estimation *", y.scale), xlab = "",
            main="Est for different group sizes",
            xaxt='n')
    #axis for the n.seq[rev(ind)]
    legend("top", 
           c("TimOrd.min  ", "TimOrd.max  ", 
             "ValOrd.min  ", "ValOrd.max  ",
             "par.true"), 
           col=c(2,2,4,4,"brown"), lty=1:5, 
           lwd=2, cex=0.4)
    axis(1, at = sum2 - n.lab, labels = n.lab, line = 1, 
         col = 2, col.ticks = 2, col.axis = 2)
    mtext("n.rev", 1, line = 1, at = -40, col = 2)
    #axis for the n.seq[ind]
    axis(1, at = n.lab, line = 3,
         col=4, col.ticks=4, col.axis=4)
    mtext("n", 1, line=3, at=-40, col=4)
    #par.true
    abline(h=par.true[1], col="brown", lty = 5)
    abline(h=par.true[2], col="brown", lty = 5)
  }
  n1 <- n.seq[sum(est.mat2[,2] - est.mat[rev(seq_along(n.seq)),2] > 0)]
  n2 <- n.seq[sum(est.mat[rev(seq_along(n.seq)),1] - est.mat2[,1] > 0)]
  n12 <- c(n1,n2); n12.ord <- n12[order(n12)]
  min.n <- n12.ord[1]; max.n <- n12.ord[2]
  n.fix <- n.seq[length(n.seq)]
  if(choprule=="min.n"){
    return(max.mean2(y, min.n)/y.scale)
  } else if(choprule == "max.n"){
    return(max.mean2(y, max.n)/y.scale)
  } else if(choprule == "both.n"){
    return(rbind(c(max.mean2(y, n1)/y.scale, n1), 
                 c(max.mean2(y, n2)/y.scale, n2)))
  } else if(choprule == "fix.n"){
    return(max.mean2(y, n.fix)/y.scale)
  }
}
 
###Distance Measure 
####function: dist
  
dist <- function(est){
  #est is the estimation
  #par is the true value
  sqrt(sum((est-par.true)^2))
}
 
##functions: get.meanseq

  
#non-overlapping groups
get.meanseq.novlp <- function(y, n.guess){
  m <- floor(length(y)/n.guess)
  y <- y[1:(m*n.guess)]  
  y.mat <- matrix(y, ncol = n.guess, byrow = TRUE)
  mean.seq <- apply(y.mat, 1, mean)
  mean.seq
}
 

#check the max-mean curve based on LIL

 
#overlapping groups
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
 

###function: get.groupfn.seq
  
require(zoo)
#moving overlapping window
#n.guess = 3; step = 2;
get.groupfn.seq <- function(y, n.guess, novlp.ind = FALSE, step = NULL, groupfn = mean){
  y.ts <- zoo(y)
  if(novlp.ind){
    step = n.guess
  }
  re.temp <- rollapply(y.ts, width = n.guess, by = step, FUN = groupfn, align="left")
  re <- as.numeric(re.temp)
  re
}
#non-ovlp, let n.guess=step
#y.seq <- 1:6
#get.groupfn.seq(y = y.seq, n.guess = 3, step = 2, groupfn = mean)
#get.groupfn.seq(y = y.seq, n.guess = 3, step = 3, groupfn = mean)
 

  
get.meanseq.sum <- function(y, n.guess, step=5){
  N <- length(y)
  n <- n.guess
  m.novlp <- floor(N/n)
  m.ovlp <- floor((N-n)/step) 
  if(m.novlp < m.ovlp){
    #choose ovlp
    m <- m.ovlp
    ind.mat1 <- matrix(rep(1:n, m+1), ncol = n, byrow = TRUE)
    ind.mat2 <- matrix(rep((0:m)*step, n), ncol = n, byrow = FALSE)
    ind.mat <- ind.mat1 + ind.mat2
    #y.new <- y[1:(n+m*step)]
    y.mat <- matrix(NA, ncol = n, nrow = m+1)
    for (i in seq_len(m+1)){
      y.mat[i,] <- y[ind.mat[i,]]
    }
    mean.seq <- apply(y.mat, 1, mean)
  } else {
    #choose novlp
    m <- m.novlp
    y <- y[1:(m*n)]  
    y.mat <- matrix(y, ncol = n.guess, byrow = TRUE)
    mean.seq <- apply(y.mat, 1, mean)
  }
  mean.seq
}
 

##function:G-norm CDF
  
#ref from Prof.Yang's paper
pGnorm <- function(z, sdl=.5, sdr=1){
  C1 <- 2*sdr/(sdl+sdr)
  C2 <- 2*sdl/(sdl+sdr)
  C1*pnorm(z/sdr)*(z<=0) + (1-C2*pnorm(-z/sdl))*(z>0)
}

#plot()

dGnorm <- function(z, sdl=.5, sdr=1){
  C <- sqrt(2)/(sqrt(pi)*(sdl+sdr))
  C*(exp(-z^2/(2*sdr^2))*(z<=0) + exp(-z^2/(2*sdl^2))*(z>0))
}
 
qGnorm <- function(p, sdl=.5, sdr=1){
  K <- sdr/(sdl+sdr)
  if(p<=K){
    sdr*qnorm(p/(2*K))
  } else {
    -sdl*qnorm((1-p)/(2*(1-K)))
  }
  # C1 <- sdr*qnorm(p/(2*K))
  # C2 <- sdl*qnorm((1-p)/(2*(1-K))) #this will produce error
  # C1*(p<=K)-C2*(p>K)
}

#qGnorm(0.5)
#qGnorm(0.7)
#sapply(seq(0,1,.1), function(p) pGnorm(qGnorm(p)))
#pGnorm(qGnorm(0.2))

pSemiGnorm <- function(z, sdl=.5, sdr=1){
  pnorm(z/sdl)*(z>=0) + pnorm(z/sdr)*(z<0)
}

dSemiGnorm <- function(z, sdl=.5, sdr=1){
  dnorm(z/sdl)*(z>=0) + dnorm(z/sdr)*(z<0)
}
##function:uppvalue.maxmean
  
uppvalue.maxmean <- function(c, m, sdl=0.5, sdr=1, pf = pGnorm){
  #upprob.less.c <- (pGnorm(c))^m 
  upprob.more.c <- 1- (1-pf(-c, sdl = sdl, sdr = sdr))^m
  #uppvalue.final <- upprob.less.c*(c<=0) + upprob.more.c*(c>0)
  #uppvalue.final
  upprob.more.c
}
 

  
cdf.maxmean <- function(c, m, sdl=0.5, sdr=1, pf = pGnorm){
  upprob.less.c <- (pf(c, sdl = sdl, sdr = sdr))^m
  upprob.less.c
}
#we can symbolically compute its pdf
pdf.maxmean <- function(c, m, sdl=0.5, sdr=1, pf = pGnorm, df = dGnorm){
  m*(pGnorm(c, sdl = sdl, sdr = sdr))^(m-1) * dGnorm(c, sdl = sdl, sdr = sdr)
}
 

##function:uppvalue.minmean
  
#min(Z_i)
uppvalue.minmean <- function(c, m, sdl=0.5, sdr=1){
  upprob.less.c <- 1 - (1 - pGnorm(c, sdl = sdl, sdr = sdr))^m
  upprob.less.c
}
 

  
#p.vaule may become 
#min(uppvalue.min(a), uppvalue.max(b))
#change m 
 

pdf.minmean <- function(c, m, sdl=0.5, sdr=1){
  m*(1 - pGnorm(c, sdl = sdl, sdr = sdr))^(m-1)*dGnorm(c, sdl = sdl, sdr = sdr)
}
 
#upper p-values
uppvalue.minmean <- function(c, m, sdl=0.5, sdr=1){
  upprob.less.c <- 1 - (1 - pGnorm(c, sdl = sdl, sdr = sdr))^m
  upprob.less.c
}

#combeind
uppvalue.comb <- function(a,b,m, sdl=0.5,sdr=1){
  #compute the upper p-values that Z_(1)<a or Z_(m)>b
  1-(1-pGnorm(-b,sdl=sdl,sdr=sdr)-pGnorm(a,sdl=sdl,sdr=sdr))^m
}

#cdf
#pdf

#find the critical value
#z.seq <- seq()
#pGnorm(-.2)
#pGnorm(0)
#pGnorm(.2)

#curve(pGnorm, from = -5, to = 5, n=1e3)
#curve(Vectorize(pGnorm), from = -10, to = -10)
#curve(dGnorm, from = -5, to = 5, n=1e3)

#critical value

#create plots

plotl <- function(x){
  #quickly draw a line plot 
  plot(x, type = "l")
}

  
get.name <- function(v1) {
  deparse(substitute(v1))
}

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#to combine legend
#https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplot

#deal with legend in ggplot2 
#https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels/

##function: testUncertainty

#colorblind-friendly palette
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbPalette)

#so far, it is the testUncertainty.plot
#we may also include the est and the numerical results later 
testUncertainty <- function(dat.seq, ref.seq=NULL,
                            moment.num=NULL, 
                            varphi = identity, 
                            par.true=NULL,
                            special.groupfn = FALSE,
                            mov.gr.step = 100, group.function = mean, 
                            initial.plot=FALSE, group.size.ini = 100, prop.ini = 1/4,
                            mov.plot=FALSE, n.seq = NULL,
                            n.min=20, n.max.prop=1/10, n.step=5,
                            test.plot=TRUE, test.iid.re=FALSE,
                            test.new.re=TRUE, n.check = NULL, 
                            test.new.plot = TRUE,
                            smooth.plot.test = FALSE, 
                            print.plot = TRUE, 
                            n.check.seq = NULL,  n.consistent.ind = TRUE, 
                            M.rep=5e2,
                            hist.bin=35, 
                            return.ind=FALSE, all.plot = TRUE, plot.re = FALSE, 
                            sig.level = .05, 
                            sd.vec.true = NULL){
  
  #test.plot = TRUE: includes the four visualizations of extrem-mean curves
  #test.new.plot = TRUE: include the p-values under each group size
  #1/n.max.prop is the number of the groups with largest size, we do not want to make it too small because it will not be representative (especially for the maximal distribution requiring many sample points to show its boundaries). 
  
  #smooth.plot.test: whether to add smooth lines on V-value curves
  
  N <- length(dat.seq)
  dat.name <- get.name(dat.seq)
  alpha <- sig.level
  
  #set a reminder to the user that test.new.re needs to be TRUE before making test.new.plot be TRUE
  
  #mov.gr.step 
  #consider overlapping group 
  #need to have special.groupfn = TRUE
  #to consider mov.gr.step and group.function
  if(!special.groupfn){
    get.meanseq=get.meanseq.novlp #if no special group function, use it as default
  } else {
    get.meanseq= function(y, n.guess){
      get.groupfn.seq(y, n.guess, step = mov.gr.step, groupfn = group.function)
    }
  }
  
  if(!is.null(moment.num)){
    varphi = function(x) x^moment.num
  }
  
  if(is.null(n.seq)){
    n.max <- floor(N*n.max.prop)
    n.seq <- seq(n.min,n.max,n.step)
  }
  
  n.len <- length(n.seq)
  if (is.null(ref.seq)){
    ref.seq <- sample(dat.seq, N, replace = TRUE)
    #consider bootstrap with blocking here
    #block.len = 5
    #increase the block.len
    #only use it for afterwards analysis on the inference part
    #it may be connected with Type I error
  }
  
  dat2.seq <- varphi(dat.seq) 
  ref2.seq <- varphi(ref.seq)
  dat2.meanseq.list <- ref2.meanseq.list <- vector("list", length = n.len)
  dat2.minmax.mat <- ref2.minmax.mat <- matrix(NA, nrow = n.len, ncol=2)
  for (i in seq_along(n.seq)){
    group.size <- n.seq[i]
    dat2.meanseq <- get.meanseq(dat2.seq, n.guess = group.size)
    ref2.meanseq <- get.meanseq(ref2.seq, n.guess = group.size)
    dat2.meanseq.list[[i]] <- dat2.meanseq
    ref2.meanseq.list[[i]] <- ref2.meanseq
    dat2.minmax.mat[i,] <- c(min(dat2.meanseq), max(dat2.meanseq))
    ref2.minmax.mat[i,] <- c(min(ref2.meanseq), max(ref2.meanseq))
  }
  #make the curve of critical value
  D.dat2 <- min(dat2.minmax.mat[,2]) - max(dat2.minmax.mat[,1]) #min(max-mean) -max(min-mean)
  if (all.plot==FALSE){initial.plot <- mov.plot <- test.plot <- all.plot}
  if (initial.plot){
    plot(dat.seq[1:floor(N*prop.ini)], type="l")
    hist(dat.seq, breaks = "scott")
    N1 <- floor(N/group.size.ini)*group.size.ini
    dat.seq.mat <- matrix(dat.seq[seq_len(N1)], nrow = group.size.ini, byrow = FALSE)
    boxplot(dat.seq.mat, ylab="W.seq", main=paste("Dat.seq with group size = ", group.size.ini))
    ref.seq.mat <- matrix(ref.seq[seq_len(N1)], nrow = group.size.ini, byrow = FALSE)
    boxplot(ref.seq.mat, ylab="Ref.seq", main=paste("Ref.seq with group size = ", group.size.ini))
  } 
  
  if (mov.plot){
    #design an automatic way to set the axis
    ##decide the x range
    p1 <- hist(dat2.meanseq.list[[1]], 
               freq = FALSE, plot = FALSE, warn.unused = FALSE) 
    p2 <- hist(dat2.meanseq.list[[1]], 
               freq = FALSE, plot = FALSE, warn.unused = FALSE) 
    breaks.tol <- c(p1$breaks, p2$breaks)
    x.lim <- c(min(breaks.tol), max(breaks.tol))
    ##decide the y range
    p1 <- hist(dat2.meanseq.list[[n.len]], 
               freq = FALSE, plot = FALSE, warn.unused = FALSE) 
    p2 <- hist(dat2.meanseq.list[[n.len]], 
               freq = FALSE, plot = FALSE, warn.unused = FALSE)
    #abline(v=par.true, lty=2, col=2, lwd=2)
    den.tol <- c(p1$density, p2$density)
    y.lim <- c(0, max(den.tol))
    for (i in seq_along(n.seq)){
      p1 <- hist(dat2.meanseq.list[[i]], plot = FALSE)
      p2 <- hist(ref2.meanseq.list[[i]], plot = FALSE) 
      
      plot(p1, col=rgb(0,0,1,1/4), freq = FALSE, xlim = x.lim, ylim = y.lim, border = rgb(0,0,1,1/4), main = paste("Group vars with n=", n.seq[i]))  
      #"Sampling distn of sample variance: M[1,2] vs U[1,2]"
      #, xlim=c(-1,1)*0.7
      plot(p2, col=rgb(1,0,0,1/4), freq = FALSE, border = rgb(1,0,0,1/4),
           add=TRUE) 
      abline(v=par.true, col="brown", lty=2, lwd=1.5)
    }
  }
  if (test.plot){
    points.dat1 <- data.frame(from="Our Data", maxmin="min", gr.var = dat2.minmax.mat[,1])
    points.dat2 <- data.frame(from="Our Data", maxmin="max", gr.var = dat2.minmax.mat[,2])
    points.dat3 <- data.frame(from="Ref Data", maxmin="min", gr.var = ref2.minmax.mat[,1])
    points.dat4 <- data.frame(from="Ref Data", maxmin="max", gr.var = ref2.minmax.mat[,2])
    
    #points.dat.our <- rbind(points.dat1, points.dat2)
    #points.dat.ref <- rbind(points.dat3, points.dat4)
    
    points.dat <- rbind(points.dat1, points.dat2, points.dat3, points.dat4)
    names <- paste(points.dat$from, points.dat$maxmin)
    points.dat.new <- cbind(points.dat, names, n.seq) 
    
    #maxmin-mean curve or 2M-mean curve
    #matplot(n.seq, cbind(dat2.minmax.mat, ref2.minmax.mat), 
    #           type = "l", col = c(1,1,2,2))
    
    ##print(ggplot(points.dat.new, aes(x=n.seq, y=gr.var)) + geom_point(aes(colour=names)) + geom_smooth(method="loess",formula=y~x,aes(colour=names)))
    plot.2Mmean.obj <- ggplot(points.dat.new, aes(x=n.seq, y=gr.var)) + geom_point(aes(colour=names)) + geom_smooth(method="loess",formula=y~x,aes(colour=names))
    #+geom_line
    ##print(plot.2Mmean.obj)
    #tips from Prof. Kulperger, it will be much better to have a smooth curve to show the trend when the max-mean points involves too much noise (waggling around)
    
    #scatter plot of maxmin-mean points of Our Data vs Ref Data
    plot.2Mpoints.obj <- ggplot(points.dat.new, aes(x=from, y=gr.var, colour=names)) + geom_point()
    ##print(plot.2Mpoints.obj)
    #qplot(factor(points.dat$from), points.dat$gr.var, colour=factor(names))
    
    #histogram of maxmin-mean points of Our Data vs Ref Data
    #if(is.null(par.true)){
    #  ##print(ggplot(points.dat.our, aes(x=gr.var, fill=maxmin)) + geom_histogram(position = "identity", alpha=0.8, bins = hist.bin) + labs("Our Data: Histogram of maxmean and minmean points"))
    #} else {
    #  ##print(ggplot(points.dat.our, aes(x=gr.var, fill=maxmin)) + geom_histogram(position = "identity", alpha=0.8, bins = hist.bin) + geom_line(xintercept=par.true, colour="brown") + labs("Our Data: Histogram of maxmean and minmean points"))
    #}
    
    ##print(ggplot(points.dat.ref, aes(x=gr.var, fill=maxmin)) + geom_histogram(position = "identity", alpha=0.8, bins = hist.bin) + geom_vline(xintercept=mean(ref2.seq), colour="brown") + labs("Ref Data: Histogram of maxmean and minmean points"))
    
    #put the histograms together 
    if(is.null(par.true)){
      plot.2Mhist.obj <- ggplot(points.dat, aes(x=gr.var, fill=from)) + geom_histogram(position = "identity", alpha=0.5, bins=hist.bin)
      #print(plot.2Mhist.obj)
    } else{
      plot.2Mhist.obj <- ggplot(points.dat, aes(x=gr.var, fill=from)) + geom_histogram(position = "identity", alpha=0.5, bins=hist.bin) + geom_vline(xintercept = c(mean(ref2.seq), par.true), colour=c("brown", "red", "red"))
      #print(plot.2Mhist.obj)
    }
    
    #maxmin-clustering plot
    minmax.mat.sum1 <- data.frame(min.gr.var = dat2.minmax.mat[,1], max.gr.var= dat2.minmax.mat[,2], names=factor("Our Data"))
    minmax.mat.sum2 <- data.frame(min.gr.var = ref2.minmax.mat[,1], max.gr.var= ref2.minmax.mat[,2], names=factor("Ref Data"))
    
    minmax.mat.sum <- rbind(minmax.mat.sum1, minmax.mat.sum2)
    
    #plot the points with different sizes
    #qplot(min.gr.var, max.gr.var, data = minmax.mat.sum, colour=names)
    plot.2Mcluster.obj <- ggplot(minmax.mat.sum,aes(min.gr.var, max.gr.var, colour=names)) + geom_point(alpha=0.6) + geom_abline(intercept = 0, slope = 1, colour = "red", size = 1.25)
    #print(plot.2Mcluster.obj)
    if(print.plot){
      print(plot.2Mmean.obj)
      print(plot.2Mpoints.obj)
      print(plot.2Mhist.obj)
      print(plot.2Mcluster.obj)
    }
  }
  #test with null: iid F
  if (test.iid.re){
    D.ref2.seq <- numeric(M.rep)
    D.ref2.seq <- replicate(M.rep,{
      ref2.seq <- sample(dat2.seq, N, replace = TRUE)
      ref2.minmax.mat <- matrix(NA, nrow = n.len, ncol=2)
      for (i in seq_along(n.seq)){
        group.size <- n.seq[i]
        ref2.meanseq <- get.meanseq(ref2.seq, n.guess = group.size)
        ref2.minmax.mat[i,] <- c(min(ref2.meanseq), max(ref2.meanseq))
      }
      D.ref2 <- min(ref2.minmax.mat[,2]) - max(ref2.minmax.mat[,1]) #min(max-mean) -max(min-mean)
      D.ref2
    })
    ##print(ggplot(points.dat, aes(x=gr.var, fill=from)) + geom_histogram(position = "identity", alpha=0.5, bins=hist.bin))
    D.seq <- c(D.ref2.seq, D.dat2)
    hist(D.ref2.seq, xlim = c(min(D.seq)*0.9, max(D.seq)*1.1))
    abline(v=D.dat2, col=2, lty=2, lwd=1.2)
    
    cat(paste("Under the null, the relative location of our dataset is at the prob", mean(D.ref2.seq>D.dat2))) #p-value
  }
  if (test.new.re){
    #choose the group size (from data or background)
    if(is.null(n.check)){n.check <- 2e2}
    #n.check <- 150
    m.check <- floor(N/n.check)
    #m.check <- N/n.check #to make the V-values smoother
    dat2.seq.std <- dat2.seq-mean(dat2.seq) 
    dat2.meanseq.std <- get.meanseq(dat2.seq.std, n.guess = n.check)
    #estimate the sd bounds
    if (is.null(sd.vec.true)){
      dat4.seq.std <- dat2.seq.std^2
      dat4.meanseq.std <- get.meanseq(dat4.seq.std, n.guess = n.check)
      sdl.est <- sqrt(min(dat4.meanseq.std)) #lower sd
      sdr.est <- sqrt(max(dat4.meanseq.std)) #upper sd
    } else {
      sdl.est <- sd.vec.true[1]
      sdr.est <- sd.vec.true[2]
    }
    #give the upper p-value for order statistics
    p.maxmean <- uppvalue.maxmean(sqrt(n.check)*max(dat2.meanseq.std), sdl = sdl.est, sdr = sdr.est, m=m.check)
    p.minmean <- uppvalue.minmean(sqrt(n.check)*min(dat2.meanseq.std), sdl = sdl.est, sdr = sdr.est, m=m.check)
    cat(paste("Under the varphi-certainty null,","with group size = ", n.check,",", "\n", 
              "the upper p-value of the max-mean stats is", format(p.maxmean, digits = 6), ";","\n", 
              "the upper p-value of the min-mean stats is", format(p.minmean, digits = 6)),".", "\n")
    p.new.re <- c(p.maxmean, p.minmean)
    if (test.new.plot){
      #set a reminder to the reader that test.new.re needs to be TRUE
      #print a warning for this  
      if(is.null(n.check.seq)){
        if(n.consistent.ind){
          n.check.seq <- n.seq
        } else {
          n.check.seq <- seq(1e2, floor(N/10), 20)
        }
      }
      stats.max.seq <- stats.min.seq <- numeric(length(n.check.seq))
      p.max.seq <- p.min.seq <- p.comb.seq <- numeric(length(n.check.seq))
      a.alpha.seq <- b.alpha.seq <- numeric(length(n.check.seq))
      for (i in seq_along(n.check.seq)){
        n.check <- n.check.seq[i]
        m.check <- floor(N/n.check)
        #estimate the variance bounds
        dat2.seq.std <- dat2.seq-mean(dat2.seq)
        dat2.meanseq.std <- get.meanseq(dat2.seq.std, n.guess = n.check)
        
        if (is.null(sd.vec.true)){
          dat4.seq.std <- dat2.seq.std^2
          dat4.meanseq.std <- get.meanseq(dat4.seq.std, n.guess = n.check)
          sdl.est <- sqrt(min(dat4.meanseq.std)) #lower sd
          sdr.est <- sqrt(max(dat4.meanseq.std)) #upper sd
        } else {
          sdl.est <- sd.vec.true[1]
          sdr.est <- sd.vec.true[2]
        }
        
        #order statistics of the group average of the centralized value
        stats.max.seq[i] <- max(dat2.meanseq.std) 
        stats.min.seq[i] <- min(dat2.meanseq.std)
        a <- sqrt(n.check)*stats.min.seq[i]
        b <- sqrt(n.check)*stats.max.seq[i]
        #give the upper p-value for order statistics
        p.max.seq[i] <- uppvalue.maxmean(b, sdl = sdl.est, sdr = sdr.est, m=m.check) #change m.check to N/n.check
        p.min.seq[i] <- uppvalue.minmean(a, sdl = sdl.est, sdr = sdr.est, m=m.check)
        #compute the combined upper p-values
        p.comb.seq[i] <- uppvalue.comb(a,b, sdl = sdl.est, sdr = sdr.est, m=m.check)
        #also compute the critical value
        crit.value <- qGnorm(1-(1-alpha)^(n.check/N), sdl = sdl.est, sdr = sdr.est)
        a.alpha.seq[i] <- crit.value/sqrt(n.check)
        b.alpha.seq[i] <- -crit.value/sqrt(n.check)
      }
      p.maxmin.dat1 <- data.frame(groupsize.check = n.check.seq, 
                                  uppvalue.maxmean = p.max.seq, 
                                  uppvalue.minmean = p.min.seq)
                                  #uppvalue.comb = p.comb.seq)
      p.maxmin.dat <- melt(p.maxmin.dat1, id.vars = "groupsize.check")
      # plot.test.obj <- ggplot(data = p.maxmin.dat, aes(x=groupsize.check, y=value))  + 
      #   geom_point(aes(colour=variable)) + geom_smooth(method="loess",formula=y~x,aes(colour=variable)) + 
      #   geom_hline(yintercept = c(0.1,0.05), col = "brown", lty=2) + ggtitle("V-values")
      # 
      if (smooth.plot.test){
        plot.test.obj <- ggplot(data = p.maxmin.dat, aes(x=groupsize.check, y=value))  + 
          geom_point(aes(colour=variable))  + 
          geom_hline(yintercept = c(0.1,0.05), col = "brown", lty=2) + ggtitle("V-values") +
          #geom_line(aes(colour=variable))  + 
          geom_smooth(method="loess",formula=y~x,aes(colour=variable, linetype = variable), se = FALSE)
      } else {
        plot.test.obj <- ggplot(data = p.maxmin.dat, aes(x=groupsize.check, y=value))  + 
          #geom_point(aes(colour=variable))  + 
          geom_hline(yintercept = c(0.1,0.05), col = "brown", lty=2) + ggtitle("V-values") +
          geom_line(aes(colour=variable))
          #geom_smooth(method="loess",formula=y~x,aes(colour=variable, linetype = variable), se = FALSE)
      }
      
      #check it
      #print(plot.test.obj)
      # + geom_line(aes(lty=variable))  #check it
      critval.maxmin.dat1 <- data.frame(groupsize.check = n.check.seq, 
                                        value.maxmean = stats.max.seq, 
                                        value.minmean = stats.min.seq,
                                  critval.maxmean.adj = b.alpha.seq, #adjust by timing 1/sqrt(n.check)
                                  critval.minmean.adj = a.alpha.seq)
      #melt the data
      crit.val.dat <- melt(critval.maxmin.dat1, id.vars = "groupsize.check")
     # plot.crit.obj0 <- ggplot(data = crit.val.dat, aes(x=groupsize.check, y=value)) + geom_smooth(method="loess",formula=y~x) + scale_linetype_manual(values = c(1,1,2,2))
      #legend1 <- g_legend(plot.crit.obj0)
      plot.crit.obj <- ggplot(data = crit.val.dat, aes(x=groupsize.check, y=value)) +
        geom_point(data = subset(crit.val.dat, variable == "value.maxmean" | variable== "value.minmean"), aes(colour=variable)) +
        geom_smooth(method="loess",formula=y~x,aes(colour=variable, linetype = variable), se = FALSE) +
        scale_linetype_manual(values = c(1,1,2,2), guide = FALSE) +
        ggtitle("Critical values and test statistics")
      # plot.crit.obj <- ggplot(data = crit.val.dat, aes(x=groupsize.check, y=value)) + 
      #   geom_smooth(method="loess",formula=y~x, aes(colour=variable, linetype = "dashed"))
      #   
      #print(plot.crit.obj)
    }
    #return(p.new.re) #return this pvalue vector under n.check
    if(print.plot){
      print(plot.test.obj)
      print(plot.crit.obj)
    }
  }
  if(return.ind){
    if(plot.re){
      return(list(plot.test.obj = plot.test.obj,
                  plot.2Mmean.obj = plot.2Mmean.obj, 
                  plot.2Mpoints.obj = plot.2Mpoints.obj, 
                  plot.2Mhist.obj = plot.2Mhist.obj,
                  plot.2Mcluster.obj = plot.2Mcluster.obj,
                  plot.testCV.obj = plot.crit.obj))
    } else {
      return(list(dat2.meanseq.list = dat2.meanseq.list, 
                  ref2.meanseq.list = ref2.meanseq.list,
                  dat2.minmax.mat = dat2.minmax.mat,
                  ref2.minmax.mat = ref2.minmax.mat, 
                  p.new.re = p.new.re)
      )
    }
  }
}
 

StrechSmalls <- log1p
 

  
#write a generic function for test of uncertainty
#we may apply any general operation to each group


 
