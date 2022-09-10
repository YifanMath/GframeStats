# pseudo simulation of sequential indep for a small family of varphi

rbisemiGnorm.seqind <- function(n, m, set.uncertainty = "core", 
                                sdl = 1, sdr = 2){
  #set.uncertainty = "core": only go through the four important scenarios (cover the cases under varphi in the space of homogeneous polylomial functions)
  #set.uncertainty = "full": go through a dense subset of the set of scenarios (\bar L), which includes all bounded functions
  #n: the length of one realization of one possible scenario
  #m: the number of scenarios to be (repeatedly) realized
  #generate a matrix
  if (set.uncertainty == "core"){
    #set.sd <- matrix(c(sdl, sdl, 
    #                   sdl, sdr, 
    #                   sdr, sdl,
    #                   sdr, sdr), ncol = 2, byrow = TRUE)
    set.sd <- c(sdl, sdr)
    x.array <- replicate(m, {
      #switch to one scenario
      #to do the future scenario analysis, 
      #here we simply equally treat each scenario
      ind.scene <- sample(c(1,2), 3, replace = TRUE) #make this switching rule flexible
      #generate one scenario (a bivariate sequence)
      #it also help us think how to consider the set for three variables under sequential independence
      sd.vec <- set.sd[ind.scene]
      theta <- function(s, B1, switch.p = 0,
                        sd1 = sd.vec[1], 
                        sd2 = sd.vec[2], 
                        sd3 = sd.vec[3]){
        if(0 <= s & s <=1){
          re <- sd1
        } else {
          re <- ifelse(B1>switch.p, sd2, sd3)
        }
      }
      dB1 <- rnorm(n)
      dB2 <- rnorm(n)
      #dB1 <- rnorm(n)/sqrt(2)
      #dB2 <- rnorm(n)/sqrt(2)
      B1 <- dB1
      B2 <- dB1+dB2
      bi.seq <- matrix(NA, nrow = 2, ncol = n)
      for (i in seq_len(n)){
        bi.seq[, i] <- c(theta(0.5, B1[i])*B1[i], 
                         theta(1.5, B1[i])*(B2[i]-B1[i]))
      }
      bi.seq
    })
  } else if (set.uncertainty == "full") {
    #go through more elements in the representation
    #change the switch.p 
    #arbitrarily change p 
    set.sd <- c(sdl, sdr)
    x.array <- replicate(m, {
      #switch to one scenario
      #to do the future scenario analysis, 
      #here we simply equally treat each scenario
      ind.scene <- sample(c(1,2), 3, replace = TRUE) #make this switching rule flexible
      #generate one scenario (a bivariate sequence)
      #it also help us think how to consider the set for three variables under sequential independence
      sd.vec <- set.sd[ind.scene]
      #also change switch.p
      #p1 <- sample(c(-sdr,-sdl, 0, sdl,sdr), 1)
      p1 <- runif(1, -sdr, sdr)
      theta <- function(s, B1, switch.p = p1,
                        sd1 = sd.vec[1], 
                        sd2 = sd.vec[2], 
                        sd3 = sd.vec[3]){
        if(0 <= s & s <=1){
          re <- sd1
        } else {
          re <- ifelse(B1>switch.p, sd2, sd3)
        }
        #add the return
        re
      }
      dB1 <- rnorm(n)
      dB2 <- rnorm(n)
      #dB1 <- rnorm(n)/sqrt(2)
      #dB2 <- rnorm(n)/sqrt(2)
      B1 <- dB1
      B2 <- dB1+dB2
      bi.seq <- matrix(NA, nrow = 2, ncol = n)
      for (i in seq_len(n)){
        bi.seq[, i] <- c(theta(0.5, B1[i])*B1[i], 
                         theta(1.5, B1[i])*(B2[i]-B1[i]))
      }
      bi.seq
    })
  } else {
    #warning
    x.array <- NULL
  }
  
  return(x.array)
  #m groups
  #n is the group size 
  #in each group, one prob rule is used to generate a bivariate random vector
}

rbisemiGnorm.semiseqind <- function(n, m, set.uncertainty = "core", 
                                sdl = 1, sdr = 2){
  #set.uncertainty = "core": only go through the four important scenarios (cover the cases under varphi in the space of homogeneous polylomial functions)
  #set.uncertainty = "full": go through a dense subset of the set of scenarios (\bar L), which includes all bounded functions
  #n: the length of one realization of one possible scenario
  #m: the number of scenarios to be (repeatedly) realized
  #generate a matrix
  if (set.uncertainty == "core"){
    #set.sd <- matrix(c(sdl, sdl, 
    #                   sdl, sdr, 
    #                   sdr, sdl,
    #                   sdr, sdr), ncol = 2, byrow = TRUE)
    set.sd <- c(sdl, sdr)
    x.array <- replicate(m, {
      #switch to one scenario
      #to do the future scenario analysis, 
      #here we simply equally treat each scenario
      ind.scene <- sample(c(1,2), 2, replace = TRUE) #make this switching rule flexible
      #generate one scenario (a bivariate sequence)
      #it also help us think how to consider the set for three variables under sequential independence
      sd.vec <- set.sd[ind.scene]
      theta <- function(s, 
                        sd1 = sd.vec[1], 
                        sd2 = sd.vec[2]){
        if(0 <= s & s <=1){
          re <- sd1
        } else {
          re <- sd2
        }
        re
      }
      dB1 <- rnorm(n)
      dB2 <- rnorm(n)
      #dB1 <- rnorm(n)/sqrt(2)
      #dB2 <- rnorm(n)/sqrt(2)
      B1 <- dB1
      B2 <- dB1+dB2
      bi.seq <- matrix(NA, nrow = 2, ncol = n)
      for (i in seq_len(n)){
        bi.seq[, i] <- c(theta(0.5)*B1[i], 
                         theta(1.5)*(B2[i]-B1[i]))
      }
      bi.seq
    })
  } else if (set.uncertainty == "full") {
    #go through more elements in the representation
    #set.sd <- seq(sdl, sdr)
    x.array <- replicate(m, {
      #switch to one scenario
      #to do the future scenario analysis, 
      #here we simply equally treat each scenario
      #ind.scene <- sample(c(1,2), 2, replace = TRUE) #make this switching rule flexible
      #sd.vec <- set.sd[ind.scene]
      #generate one scenario (a bivariate sequence)
      #it also help us think how to consider the set for three variables under sequential independence
      sd.vec <- runif(2, sdl, sdr)
      theta <- function(s, 
                        sd1 = sd.vec[1], 
                        sd2 = sd.vec[2]){
        if(0 <= s & s <=1){
          re <- sd1
        } else {
          re <- sd2
        }
        re
      }
      dB1 <- rnorm(n)
      dB2 <- rnorm(n)
      #dB1 <- rnorm(n)/sqrt(2)
      #dB2 <- rnorm(n)/sqrt(2)
      B1 <- dB1
      B2 <- dB1+dB2
      bi.seq <- matrix(NA, nrow = 2, ncol = n)
      for (i in seq_len(n)){
        bi.seq[, i] <- c(theta(0.5)*B1[i], 
                         theta(1.5)*(B2[i]-B1[i]))
      }
      bi.seq
    })
  } else {
    #warning
    x.array <- NULL
  }
  
  return(x.array)
  #m groups
  #n is the group size 
  #in each group, one prob rule is used to generate a bivariate random vector
}

seq2mat.mov.ovlp <- function(x.seq, gr.size, step, gr.num = NULL){
  N <- length(x.seq)
  n <- gr.size
  m <- gr.num <- floor((N-n)/step)
  ind.mat1 <- matrix(rep(1:n, m+1), ncol = n, byrow = TRUE)
  ind.mat2 <- matrix(rep((0:m)*step, n), ncol = n, byrow = FALSE)
  ind.mat <- ind.mat1 + ind.mat2
  x.mat <- matrix(NA, ncol = n, nrow = m+1)
  for (i in seq_len(m+1)){
    x.mat[i,] <- x.seq[ind.mat[i,]]
  }
  return(x.mat)
}

cummean <- function(x) {cumsum(x)/seq_along(x)}

#illustrate SLLN
#used to illustrate the asymmetry of independence 
slln.check.f <- function(f, x.seq, gr.size, step, add.plot = FALSE, 
                         prop.ex = 1/25, main, 
                         par.true = NULL, 
                         maxmin.mean.plot = TRUE,
                         save.plot = FALSE, 
                         filename = FALSE, 
                         ylim = NULL, 
                         width = 840, 
                         height = 550, 
                         cex.lab = 1.5, 
                         cex.main = 1.5){
  #prop.ex: proportion to exclude at the beginning (to visualize)
  #main: main title
  #concentrate.bound: user needs to specify two bounds in theory in function form
  #maxmin.mean.plot: whether to draw the lines plot of max and min mean for each group size n
  #filename: name of the file to be saved, 
  #if user wants to add the file to a path, add the path at the beginning of the file name
  #such as "Figs/myfile"
  x.mat <- seq2mat.mov.ovlp(x.seq, 
                            gr.size = gr.size, step = step)
  cummean.f <- function(x) cummean(f(x))
  cummeanf.mat <- apply(x.mat, 1, cummean.f)
  ind.ex1 <- seq_len(floor(gr.size*prop.ex))
  n.seq.old <- seq_len(gr.size)
  maxmin.mean.mat <- apply(cummeanf.mat, 1, function(x) c(min(x),max(x)))
  maxmin.mean.mat <- t(maxmin.mean.mat)
  if(add.plot){
    cols <- 1:6 
    matplot(n.seq.old[-ind.ex1],
            cummeanf.mat[-ind.ex1, ], type = "l", 
            xlab = "sample size", ylab = "cumavg",
            col = alpha(cols, 0.42),
            main = main, 
            ylim = ylim)
    abline(h=par.true, col = "red", lwd = 1.8)
    if(maxmin.mean.plot){
      matplot(n.seq.old[-ind.ex1],
              maxmin.mean.mat[-ind.ex1, ], type = "l", 
              xlab = "sample size", ylab = "maxmin.mean",
              main = main, 
              ylim = ylim, 
              lty = 1)
      abline(h=par.true, col = "red", lwd = 1.5, lty = 2)
    }
  }
  if(save.plot){
    #pdf(file = paste0(filename,".pdf"), width = 7, height = 5)
    png(file = paste0(filename,".png"), width = width, height = height)
    cols <- 1:6 
    matplot(n.seq.old[-ind.ex1],
            cummeanf.mat[-ind.ex1, ], type = "l", 
            xlab = "sample size", ylab = "cumavg",
            col = alpha(cols, 0.42),
            main = main, 
            ylim = ylim, 
            cex.lab = cex.lab, 
            cex.main = cex.main)
    abline(h=par.true, col = "red", lwd = 1.8)
    dev.off()
    if(maxmin.mean.plot){
      png(file = paste0(filename,"-maxmin_mean.png"), width = width, height = height)
      matplot(n.seq.old[-ind.ex1],
              maxmin.mean.mat[-ind.ex1, ], type = "l", 
              xlab = "sample size", ylab = "maxmin.mean",
              main = main, 
              ylim = ylim, 
              lty = 1)
      abline(h=par.true, col = "red", lwd = 1.5, lty = 2)
      dev.off()
    }
  }
  return(list(n.seq = n.seq.old[-ind.ex1], 
              ind.ex1 = ind.ex1,
              par.true = par.true, 
              cummeanf.mat = cummeanf.mat, 
              maxmin.mean.mat = maxmin.mean.mat))
}

#long sequence version

number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

theta.choose1 <- function(x){
  floor(7*exp(x+1)) %% 4
}

#v1: threshold = 1, info.pre = y.seq[i-1]
#v2: threshold = 2, info.pre = sum(y.seq[i-1])/sqrt(i-1) 
#v2 does have blocking design feature

theta.choose2 <- function(x, threshold = 0.5){
  #v1: threshold=1
  #v2: threshold=0.5, 2
  #monotonically depend on the abs(x)
  #threshold>0
  a <- threshold
  #if the previous value is small
  #then the next stage will tend to make it jump into low-vol
  #if the previous value is large
  #then the next stage will tend to make it jump into high-vol
  #2*(x<=-a) + 1*(-a<x & x<=0) + 3*(0<x & x<=a) + 4*(x>a)
  1*(x<=-a) + 2*(-a<x & x<=0) + 0*(0<x & x<=a) + 3*(x>a)
}

#v3: threshold = 10, info.pre = S.seq[i-1] 
theta.choose3 <- function(x, threshold = 4){
  a <- threshold
  floor(abs(x)/a) %% 4
}

#theta.choose3(seq(-20,-80,-2))


# Functions

## set up the default rule here
#the current sig is given 
#choose the prob rule
#then generate the next sig
prob.sigRS.rule <- function(sig.now, ep.now, 
                            info.pre, 
                            theta.choose = theta.choose2, 
                            sdl = sdl0, sdr = sdr0,
                            switch.loc = 0){
  #info.now = c(sig.now, ep.now)
  #info.pre: the previous information that many have effects on the future theta
  #theta.nextperiod = theta.choose(info.pre)
  #info.pre could a subset of (sig.seq[s], ep.seq[s], s<=i) 
  #theta.choose: the theta(.) function characterizing the dependence
  #sdl = sdl0; sdr = sdr0; the default one is exactly the global one
  #switch.loc = 0 by default, user may change it according to their need
  y.now <- sig.now*ep.now
  #theta.nextperiod is the choice of model in the next period
  theta.nextperiod <- theta.choose(info.pre)
  ind.scene <- number2binary(theta.nextperiod, noBits = 2) + 1 
  set.sd <- c(sdl,sdr)
  sd.vec <- set.sd[ind.scene]
  sig.next <- ifelse(y.now>switch.loc, sd.vec[1], sd.vec[2])
  re <- list(theta.nextperiod = theta.nextperiod, sig.next = sig.next)
  return(re)
}


#VaR defines the capital to be prepared to support one's position
#it is a risk measure 
#it is positive
VaR <- function(x, alpha = .01){
  #x is the dataset of gain or return
  -quantile(x, probs = alpha)
}
#VaR(rnorm(1e3))
#also compute ES

rvol.seq.long <- function(N, set.sd){
  sdl0 <- set.sd[1]
  sdr0 <- set.sd[2]
  
  sig.seq <- y.seq <- S.seq <- numeric(N)
  theta.seq <- numeric(N-1) #theta.seq[i] is used from time i to i+1
  #Initialization 
  
  #time t=1
  sig.seq[1] <- sample(set.sd, 1)
  ep.seq <- rnorm(N) #our main focus is the sigma sequence 
  theta.seq[1] <- sample(c(0,1,2,3),1)
  y.seq[1] <- sig.seq[1]*ep.seq[1]
  S.seq[1] <- y.seq[1] 
  
  #from t=1 to t=2
  ind.scene <- number2binary(theta.seq[1], noBits = 2) + 1 
  #transfer 1 to (1,1), to check, run the code below
  #sapply(0:3, function(x) number2binary(x, noBits = 2) + 1)
  sd.vec <- set.sd[ind.scene] #choose from the four possible scenarios
  sig.seq[2] <- ifelse(y.seq[1]>0, sd.vec[1], sd.vec[2])
  y.seq[2] <- sig.seq[2]*ep.seq[2]
  S.seq[2] <- sum(y.seq[seq_len(2)])
  
  #from time 3 to N
  for (i in seq(2,N-1,1)){
    #get the sig at time i+1
    #we also compute the theta.seq[i] which is used from time i to i+1
    #from time i to i+1, or write it as [i,i+1]
    re <- prob.sigRS.rule(sig.now = sig.seq[i], 
                          ep.now = ep.seq[i], 
                          info.pre = y.seq[i-1], 
                          theta.choose = theta.choose2)
    theta.seq[i] <- re$theta.nextperiod #record the theta used in [i,i+1]
    sig.seq[i+1] <- re$sig.next
    y.seq[i+1] <- sig.seq[i+1]*ep.seq[i+1]
    S.seq[i+1] <- sum(y.seq[seq_len(i+1)])
  }
  
  #dat.history <- data.frame(y.seq = y.seq, sig.seq=sig.seq, 
  #                          S.seq=S.seq, ep.seq=ep.seq)
  
  list.history <- list(y.seq = y.seq, 
                       sig.seq=sig.seq, 
                       ep.seq=ep.seq, 
                       S.seq=S.seq, 
                       theta.seq=theta.seq)
  return(list.history)
}

rvol.seq.long2 <- function(N, set.sd){
  #version 2 of the generation: 
  #this version has stronger volatility-clustering feature 
  sdl0 <- set.sd[1]
  sdr0 <- set.sd[2]
  
  sig.seq <- y.seq <- S.seq <- numeric(N)
  theta.seq <- numeric(N-1) #theta.seq[i] is used from time i to i+1
  #Initialization 
  
  #time t=1
  sig.seq[1] <- sample(set.sd, 1)
  ep.seq <- rnorm(N) #our main focus is the sigma sequence 
  theta.seq[1] <- sample(c(0,1,2,3),1)
  y.seq[1] <- sig.seq[1]*ep.seq[1]
  S.seq[1] <- y.seq[1] 
  
  #from t=1 to t=2
  ind.scene <- number2binary(theta.seq[1], noBits = 2) + 1 
  #transfer 1 to (1,1), to check, run the code below
  #sapply(0:3, function(x) number2binary(x, noBits = 2) + 1)
  sd.vec <- set.sd[ind.scene] #choose from the four possible scenarios
  sig.seq[2] <- ifelse(y.seq[1]>0, sd.vec[1], sd.vec[2])
  y.seq[2] <- sig.seq[2]*ep.seq[2]
  S.seq[2] <- sum(y.seq[seq_len(2)])
  
  #from time 3 to N
  for (i in seq(2,N-1,1)){
    #get the sig at time i+1
    #we also compute the theta.seq[i] which is used from time i to i+1
    #from time i to i+1, or write it as [i,i+1]
    re <- prob.sigRS.rule(sig.now = sig.seq[i], 
                          ep.now = ep.seq[i], 
                          info.pre = S.seq[i-1], 
                          theta.choose = theta.choose3)
    theta.seq[i] <- re$theta.nextperiod #record the theta used in [i,i+1]
    sig.seq[i+1] <- re$sig.next
    y.seq[i+1] <- sig.seq[i+1]*ep.seq[i+1]
    S.seq[i+1] <- sum(y.seq[seq_len(i+1)])
  }
  
  #dat.history <- data.frame(y.seq = y.seq, sig.seq=sig.seq, 
  #                          S.seq=S.seq, ep.seq=ep.seq)
  
  list.history <- list(y.seq = y.seq, 
                       sig.seq=sig.seq, 
                       ep.seq=ep.seq, 
                       S.seq=S.seq, 
                       theta.seq=theta.seq)
  return(list.history)
}

rvol.seq.long.VaR <- function(N, set.sd){
  sdl0 <- set.sd[1]
  sdr0 <- set.sd[2]
  
  sig.seq <- y.seq <- S.seq <- numeric(N)
  theta.seq <- numeric(N-1) #theta.seq[i] is used from time i to i+1
  #Initialization 
  
  #time t=1
  sig.seq[1] <- sample(set.sd, 1)
  ep.seq <- rnorm(N) #our main focus is the sigma sequence 
  theta.seq[1] <- sample(c(0,1,2,3),1)
  y.seq[1] <- sig.seq[1]*ep.seq[1]
  S.seq[1] <- y.seq[1] 
  
}

par.est.VaR <- function(x, alpha = c(.01, .05, .1)){
  re.est <- VaR(x, alpha = alpha)
  names(re.est) <-  c("VaR.01", "VaR.05", "VaR.1")
  re.est
}

par.est.full <- function(x){
  #user specifies the functions they are interested in 
  re.est <- c(mean(x), mean(x^2), mean(x^3), 
              VaR(x, alpha = c(.01, .05, .1)))
  names(re.est) <- c("1st Moment", "2nd Moment", 
                     "3rd Moment", 
                     "VaR.01", "VaR.05", "VaR.1")
  re.est
}
#par.est.full(rnorm(1e3))


MC.function <- function(x) {
  re <- c(mean(x), sd(x))
  names(re) <- c("MC.mean", "MC.error")
  return(re)
}

#par.est(rnorm(1e3))
#re.est <- par.est(rnorm(1e3))
#cat(paste("The estimated", names(re.est), "is", re.est), sep = "\n")

#cat("The estimated", names(re.est), "is", re.est, "\n")

#True future 

summary.future <- function(list.history, 
                           sample.size = 1e3, 
                           MC.size = 5e2, 
                           plot.ind = TRUE,
                           print.ind = FALSE,
                           MC.ind = FALSE, 
                           theta.choose = theta.choose2, 
                           par.est = par.est.full, 
                           silent.ind = FALSE,
                           f.ymat = sum){
  if(silent.ind) plot.ind <- print.ind <- FALSE
  #par.est = par.est will give this error
  #promise already under evaluation: recursive default argument reference or earlier problems?
  #use different function name
  y.seq <- list.history$y.seq 
  N <- length(y.seq)
  theta.new <- theta.choose(y.seq[c(N-1,N)])
  #print(theta.choose2(y.seq[c(N-1,N)]))
  y.new.mat <- replicate(sample.size, {
    y.seq.ext <- c(list.history$y.seq, numeric(2)) 
    S.seq.ext <- c(list.history$S.seq, numeric(2)) 
    sig.seq.ext <- c(list.history$sig.seq, numeric(2)) 
    #theta.seq.ext <- c(list.history$theta.seq, numeric(2)) 
    ep.seq.ext <- c(list.history$ep.seq, rnorm(2)) 
    #use the ep.seq in the history
    for (i in c(N,N+1)){ 
      re <- prob.sigRS.rule(sig.now = sig.seq.ext[i], 
                            ep.now = ep.seq.ext[i], 
                            info.pre = y.seq.ext[i-1]) 
      sig.seq.ext[i+1] <- re$sig.next 
      #theta.seq.ext[i] <- re$theta.nextperiod 
      #record the theta used in [i,i+1]
      y.seq.ext[i+1] <- sig.seq.ext[i+1]*ep.seq.ext[i+1]
      S.seq.ext[i+1] <- sum(y.seq.ext[seq_len(i+1)])
    } 
    y.seq.new <- y.seq.ext[c(N+1, N+2)]
    #theta.new <- theta.seq.ext[c(N, N+1)]
    #print(theta.new)
    y.seq.new 
  })
  ## compute (y[N+1]+y[N+2])
  y.new.sum <- apply(y.new.mat, 2, f.ymat)
  ## draw a histogram
  if (plot.ind) {
    hist(y.new.sum, probability = TRUE, breaks = "scott")
    lines(density(y.new.sum), col = 2)
    }
  
  if (!MC.ind){
    ## compute the k-th moment, k=1,2,3
    moment.est <- par.est(y.new.sum)
    re <- list(y.new.mat = y.new.mat, y.new.sum = y.new.sum,
               moment.est = moment.est, theta.new = theta.new)
    
    if(print.ind){
      print(paste("Given the current path history, the true future is", theta.new[1], theta.new[2])) 
      cat(paste("The estimated", names(moment.est), 
                "is", format(moment.est, digits = 5)), sep = "\n")
    }
    
    return(re)
  } else {
    #if use Monte Carlo to estimate the moment
    #we need to re-generate the sample MC.size times
    y.now <- y.seq[N]
    switch.loc <- 0
    y.seq.new <- numeric(2)
    
    moment.est.mat <- replicate(MC.size, {
      y.new.mat <- replicate(sample.size, {
        ep.seq.new <- rnorm(2)
        for (i in c(1,2)){
          ind.scene <- number2binary(theta.new[i], noBits = 2) + 1
          sd.vec <- set.sd[ind.scene]
          sig.next <- ifelse(y.now>switch.loc, sd.vec[1], sd.vec[2])
          y.seq.new[i] <- sig.next*ep.seq.new[i]
          y.now <- y.seq.new[i]
        }
        y.seq.new
      })
      y.new.sum <- apply(y.new.mat, 2, sum)
      ## compute the k-th moment, k=1,2,3 and VaR
      moment.est <- par.est(y.new.sum)
      moment.est
    })
    
    moment.est.MC <- apply(moment.est.mat, 1, MC.function)
    re <- list(y.new.mat = y.new.mat, y.new.sum = y.new.sum,
               moment.est.MC = moment.est.MC)
    if(print.ind){
      moment.est.mat.cat <- format(moment.est.MC, digits = 5)
      cat("Given the current path history, the true future is", 
          theta.new[1], theta.new[2],"\n")
      print(moment.est.MC)
      # cat("Given the current path history, the true future is", 
      #     theta.new[1], theta.new[2], "\n", 
      #     "The estimated 1st Moment is", paste0(moment.est.mat.cat[1,1],"(",moment.est.mat.cat[2,1],")"), "\n",
      #     "The estimated 2nd Moment is", paste0(moment.est.mat.cat[1,2],"(",moment.est.mat.cat[2,2],")"), "\n",
      #     "The estimated 3rd Moment is", paste0(moment.est.mat.cat[1,3],"(",moment.est.mat.cat[2,3],")"), "\n",
      #     "The estimated VaR01 is",      paste0(moment.est.mat.cat[1,4],"(",moment.est.mat.cat[2,4],")")
      # )
    }
    
    return(re)
  }
}


#fit a HMM model (the most suitable one)
# consider a HMM setup, 


#HMM
#without blocking design 

my.sim.hmm <- function (n, initial, Pi, pm) {
  #generate state
  x <- y <- numeric(n)
  x[1] <- sample(c(1,2), 1, prob = Pi[1,])*as.numeric(initial==1) + sample(c(1,2), 1, prob = Pi[2,])*as.numeric(initial==2)
  for (i in seq_len(n-1)){
    x[i+1] <- sample(c(1,2), 1, prob = Pi[1,])*as.numeric(x[i]==1) + sample(c(1,2), 1, prob = Pi[2,])*as.numeric(x[i]==2)
  }
  #generate observation
  y <- rnorm(n, mean = pm$mean[1], sd=pm$sd[1])*as.numeric(x==1) + rnorm(n, mean = pm$mean[2], sd=pm$sd[2])*as.numeric(x==2) 
  #rname <- paste("r", distn, sep = "")
  return(list(x = x, y = y))
}

hist.den <- function(x){
  hist(x, probability = TRUE, breaks = "scott")
  lines(density(x), col=2)
}

#overlay two density using ggplot2
#den.compare <- function(x1.seq,x2.seq){
# 
#}

hmm.fit.predict <- function(dat.seq, plot.ind = TRUE, 
                            n.pred = 1e3){
  #[]also consider second order hmm 
  N <- length(dat.seq)
  dat <- data.frame(y=dat.seq)
  suppressMessages(mod.hmm <- depmix(y ~ 1, family = gaussian(), nstates = 2, 
                    data = data.frame(y=dat.seq)))
  suppressMessages(hmmfit <- fit(mod.hmm, verbose = FALSE))
  re.post <- posterior(hmmfit)
  #AIC(hmmfit)
  #BIC(hmmfit)
  #plot(re$state, type = "l")
  state.seq <- re.post$state
  #simulate from a given state
  re2 <- forwardbackward(hmmfit)
  #matplot(re2$gamma[seq_len(1e2),], type = "l")
  #sim.re <- simulate(mod.hmm, nsim = 1, 
  #                   times = as.logical(c(rep(0,N),1,1)))
  #sim.seq <- as.numeric(sim.re@response[[1]][[1]]@y)
  #sim.seq2 <- as.numeric(sim.re@response[[2]][[1]]@y) 
  #to compare
  #response.list <- sim.re@response
  #predict the future values 
  #a more straightforward way
  pars.fit <- getpars(hmmfit)
  mat.vec <- pars.fit[3:6]
  #mat.vec <- c(hmmfit@transition[[1]]@parameters$coefficients, 
  #             hmmfit@transition[[2]]@parameters$coefficients)
  trans.mat <- matrix(mat.vec, ncol = 2, byrow = TRUE)
  y.new.mat <- replicate(n.pred, {
  re.new <- my.sim.hmm(2, Pi = trans.mat, 
                          pm = list(mean = pars.fit[c(7,9)], 
                                   sd = pars.fit [c(8,10)]), 
                          initial = state.seq[N])
    
  re.new$y
  })
  y.new.mat <- t(y.new.mat)
  y.new.sum <- apply(y.new.mat, 1, sum)
  
  if(plot.ind) hist.den(y.new.sum)
  
  #make a summary of the future, 
  #compare with the true future 
  mat.vec <- pars.fit[3:6]
  trans.mat <- matrix(mat.vec, ncol = 2, byrow = TRUE)
  return(list(y.new.mat = y.new.mat, 
              y.new.sum = y.new.sum, 
              state.seq = state.seq, 
              gamma.mat = re2$gamma, 
              pars.fit = pars.fit, 
              trans.mat = trans.mat, 
              initial.sim = state.seq[N]))
}


summary.future.hmm <- function(dat.seq, 
                           sample.size = 1e3, 
                           MC.size = 5e2, 
                           plot.ind = TRUE, 
                           print.ind = TRUE, 
                           MC.ind = FALSE, 
                           #theta.choose = theta.choose2, 
                           par.est = par.est.full, 
                           silent.ind = FALSE){
  #par.est = par.est will give this error
  #promise already under evaluation: recursive default argument reference or earlier problems?
  #use different function name
  if(silent.ind) plot.ind <- print.ind <- FALSE
  
  y.seq <- dat.seq
  N <- length(y.seq)
  re <- hmm.fit.predict(y.seq, plot.ind = FALSE,
                               n.pred = sample.size)
  y.new.mat <- re$y.new.mat
  y.new.sum <- re$y.new.sum
  pars.fit <- re$pars.fit
  initial.sim <- re$initial.sim #state.seq[N]
  trans.mat <- re$trans.mat
  ## draw a histogram
  if (plot.ind) {
    hist(y.new.sum, probability = TRUE, breaks = "scott")
    lines(density(y.new.sum), col = 2) 
    }
  
  if (!MC.ind){
    ## compute the k-th moment, k=1,2,3
    moment.est <- par.est(y.new.sum)
    re <- list(y.new.mat = y.new.mat,
               y.new.sum = y.new.sum,
               moment.est = moment.est)
    if(print.ind){
      print("HMM univariate:")
      cat(paste("The estimated", names(moment.est), 
                "is", format(moment.est, digits = 5)), sep = "\n")
    }
    #print(paste("Given the current path history, the model future is", theta.new[1], theta.new[2])) 
    return(re)
  } else {
    #if use Monte Carlo to estimate the moment
    #we need to re-generate the sample MC.size times
    moment.est.mat <- replicate(MC.size, {
      y.new.sum <- replicate(sample.size, {
        re.new <- my.sim.hmm(2, Pi = trans.mat, 
                             pm = list(mean = pars.fit[c(7,9)], 
                                       sd = pars.fit [c(8,10)]), 
                             initial = initial.sim)
        
        sum(re.new$y)
      })
      ## compute the k-th moment, k=1,2,3 and VaR
      moment.est <- par.est(y.new.sum)
      moment.est
    })
    
    moment.est.MC <- apply(moment.est.mat, 1, MC.function)
    re <- list(y.new.mat = y.new.mat,
               y.new.sum = y.new.sum, 
               moment.est.MC = moment.est.MC)
    if(print.ind){
      moment.est.mat.cat <- format(moment.est.MC, digits = 5)
      print("HMM univariate:")
      print(moment.est.MC)
      # cat("The estimated 1st Moment is", paste0(moment.est.mat.cat[1,1],"(",moment.est.mat.cat[2,1],")"), "\n",
      #     "The estimated 2nd Moment is", paste0(moment.est.mat.cat[1,2],"(",moment.est.mat.cat[2,2],")"), "\n",
      #     "The estimated 3rd Moment is", paste0(moment.est.mat.cat[1,3],"(",moment.est.mat.cat[2,3],")"), "\n",
      #     "The estimated VaR01 is",      paste0(moment.est.mat.cat[1,4],"(",moment.est.mat.cat[2,4],")")
      # )
    }
   
    return(re)
  }
}

#fit GARCH(1,1)

my.garchfit.old <- function(x, garch.order = c(1,1)){
  #x is the input sequence
  #garch.order = c(2,2)
  #dat <- data.frame(x = x)
  dat <- matrix(x, nrow=length(x))
  #spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = TRUE))
  spec = ugarchspec(variance.model = list(model = "fGARCH", garchOrder = garch.order, submodel = "GARCH", external.regressors = NULL, variance.targeting = FALSE), 
             mean.model     = list(armaOrder = c(0, 0), external.regressors = NULL), distribution.model = "norm", start.pars = list(), fixed.pars = list())
  #spec = ugarchspec()
  #show(spec)
  mod.fit <- ugarchfit(spec, data = dat, solver = "hybrid")
  #return the fitted sigma.t and parameters
  est.par <- coef(mod.fit)
  #it treats it as white noise 
  sig.org <- sigma(mod.fit)
  res <- residuals(mod.fit)
  x.fit <- fitted(mod.fit)
  #plot(as.numeric(sig.org), type = "l")
  #plot(x.fit)
  #plot(x, type = "l")
}

#we can also write our own function to fit garch(1,1)
#ref: SS9864-practice

#log likelihood for garch(1,1)
LLik <- function(para, x){ 
  n <- length(x)
  omg <- para[1]
  alpha <- para[2]
  beta <- para[3]
  sigmasq <- rep(omg^2, n) 
  l <- 0
  for (t in 2:n) {
    sigmasq[t] <- omg + alpha * x[t-1]^2 + beta * sigmasq[t-1]
    l <- l + log(dnorm(x[t], mean = 0, sd = sqrt(sigmasq[t])))
  }
  return(-l) 
}

garchfit.my0 <- function(y){
  optim.re <- nlminb(start=c(0.1,0.02,0.6),function(p) LLik(para=p,x=c(y[1],y)), 
               lower = c(0,0,0), upper = c(1,1,1))
  optim.re2 <- nlminb(start=optim.re$par,function(p) LLik(para=p,x=c(y[1],y)), 
                      lower = c(0,0,0), upper = c(1,1,1))
  est.par <- optim.re$par
  est.par
}

#the result is close to the garchFit rather than the one from rugarch
#garchfit.my0 uses MLE

#garchFit uses quasi-MLE 

#to validate these functions 

garchfit.my1 <- function(y, include.mean = FALSE){
  #options(warnings=-1)
  dat <- data.frame(y = y)
  suppressWarnings(mod.garch <- garchFit(y~garch(1,1), data = dat,
                                         include.mean = include.mean,
                                         trace = FALSE))
  # suppressWarnings(mod.garch <- garchFit(y~garch(2,2), data = dat,
  #                                        include.mean = include.mean,
  #                                        trace = FALSE))
  #Warning message:
  #Using formula(x) is deprecated when x is a character vector of length > 1.
  #Consider formula(paste(x, collapse = " ")) instead. 
  sig.org <- mod.garch@sigma.t
  est.par <- mod.garch@fit$par
  #res.org <- mod.garch@residuals
  #this is the raw residuals
  #it is equal to y.seq - est.par[1] (est.par[1] is mu)
  #summary(mod.garch)
  #plot(sig.org, type="l")
  if (!include.mean) est.par <- c(0, est.par) #for the simulation later
  return(list(sig.org = sig.org, est.par = est.par, 
              fit.re = mod.garch))
}

#consider GJR-garch
#https://cran.r-project.org/web/packages/rugarch/vignettes/Introduction_to_the_rugarch_package.pdf
#https://quant.stackexchange.com/questions/44415/gjr-garch-model-using-garchfit-function

garchfit.my <- function(y, include.mean = TRUE, var.model = "gjrGARCH",
                        garch.order = c(2,2)){
  #options(warnings=-1)
  dat <- data.frame(y = y)
  mod.spec <- ugarchspec(variance.model = list(model = var.model, garchOrder = garch.order), 
                     mean.model = list(armaOrder = c(1,1), include.mean = include.mean), 
                     distribution.model = "norm")
  mod.re <- ugarchfit(dat, spec = mod.spec)
  
  #mod.re
  #show(mod.re)
  sig.org <- sigma(mod.re)
  est.par <- coef(mod.re)
  #res.org <- mod.garch@residuals
  #this is the raw residuals
  #it is equal to y.seq - est.par[1] (est.par[1] is mu)
  #summary(mod.garch)
  #plot(sig.org, type="l")
  #if (!include.mean) est.par <- c(0, est.par) #for the simulation later
  return(list(sig.org = sig.org, est.par = est.par, 
              mod.re = mod.re))
}

#decide the order of garch model
#use AIC or BIC
#check garch(2,2)
#garch(1,2),
#garch(2,1)

#also try the GJR-GARCH

#system.time(re0 <- garchfit.my0(y.seq0))
#system.time(re1 <- garchfit.my(y.seq0))

garch.fit.predict <- function(dat.seq, plot.ind = TRUE, 
                            pred.step = 2, 
                            boot.time = 1e3,
                            fit.re = NULL, 
                            garchfit = garchfit.my, 
                            CB.ind = TRUE){
  #rugarch already provides the bootstrap step for forecasting
  
  #CB.ind = TRUE: apply conditional bootstrap (CB)
  #the estimated parameter is fixed during the bootstrap
  #we also need to select the order
  #garch.order = c(1,1)
  #p <- garch.order[1]
  #q <- garch.order[2]
  #m <- max(p,q)
  N <- length(dat.seq)
  if(is.null(fit.re)){
    #if there is no fit result, then we need to fit by ourselves
    fit.re <- garchfit(dat.seq)
  } 
  sig.org <- fit.re$sig.org
  est.par <- fit.re$est.par
  fit <- fit.re$mod.re
  #bootstrap
  #ref: Pascual(2006)
  
  if(CB.ind){
    boot.re <- ugarchboot(fit, method = c("Partial", "Full")[1], 
                          n.ahead = pred.step, n.bootpred = boot.time)
    y.new.mat <- boot.re@fseries
  } else {
    boot.re <- ugarchboot(fit, method = c("Partial", "Full")[2], 
                          n.ahead = pred.step, n.bootpred = boot.time)
    y.new.mat <- boot.re@fseries
  }
  
  y.new.sum <- apply(y.new.mat, 1, sum)
  
  if(plot.ind) hist.den(y.new.sum)
  
  #make a summary of the future, 
  #compare with the true future 
  
  return(list(y.new.mat = y.new.mat, 
              y.new.sum = y.new.sum,
              sigma.t = sig.org,
              pars.fit = est.par))
}

summary.future.garch <-function(dat.seq, 
                                           sample.size = 1e3, 
                                           MC.size = 5e2, 
                                           plot.ind = TRUE, 
                                           print.ind = FALSE, 
                                           MC.ind = FALSE, 
                                           #theta.choose = theta.choose2, 
                                           par.est = par.est.full, 
                                silent.ind = FALSE, 
                                CB.ind = TRUE){
  #par.est = par.est will give this error
  #promise already under evaluation: recursive default argument reference or earlier problems?
  #use different function name
  if(silent.ind) plot.ind <- print.ind <- FALSE
  y.seq <- dat.seq
  N <- length(y.seq)
  re <- garch.fit.predict(y.seq, plot.ind = FALSE,
                        boot.time = sample.size, CB.ind = CB.ind)
  y.new.mat <- re$y.new.mat
  y.new.sum <- re$y.new.sum
  pars.fit <- re$pars.fit
  ## draw a histogram
  if (plot.ind) {
    hist(y.new.sum, probability = TRUE, breaks = "scott")
    lines(density(y.new.sum), col = 2) 
  }
  
  if (!MC.ind){
    ## compute the k-th moment, k=1,2,3
    moment.est <- par.est(y.new.sum)
    re <- list(y.new.mat = y.new.mat,
               y.new.sum = y.new.sum,
               moment.est = moment.est)
    if(print.ind){
      print("GARCH(1,1) bootstrap:")
      cat(paste("The estimated", names(moment.est), 
                "is", format(moment.est, digits = 5)), sep = "\n")
    }
    #print(paste("Given the current path history, the model future is", theta.new[1], theta.new[2])) 
    return(re)
  } else {
    #if use Monte Carlo to estimate the moment
    #we need to re-generate the sample MC.size times
    fit.re0 <- garchfit.my(y.seq)
    
    moment.est.mat <- replicate(MC.size, {
      re <- garch.fit.predict(y.seq, plot.ind = FALSE,
                              boot.time = sample.size, 
                              fit.re = fit.re0)
      y.new.sum <- re$y.new.sum
      moment.est <- par.est(y.new.sum)
      moment.est
    })
    
    moment.est.MC <- apply(moment.est.mat, 1, MC.function)
    re <- list(y.new.mat = y.new.mat,
               y.new.sum = y.new.sum, 
               moment.est.MC = moment.est.MC)
    if(print.ind){
      moment.est.mat.cat <- format(moment.est.MC, digits = 5)
      print("GARCH univariate:")
      print(moment.est.MC)
    }
    
    return(re)
  }
}


get.mov.seq <- function(dat.seq, window.size, step = 5){
  #N total length
  y <- dat.seq
  N <- length(y)
  n <- window.size
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
  } else {
    #choose novlp
    m <- m.novlp
    ind.mat <- matrix(1:(m*n), ncol = n.guess, byrow = TRUE)
    y <- y[1:(m*n)]  
    y.mat <- matrix(y, ncol = n.guess, byrow = TRUE)
  }
  
  return(list(y.mat = y.mat, 
              ind.mat = ind.mat))
}

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

#to compute the upcdf

#gives the complete solution for u(t,x)
Gheat.itn.sol <- function(varphi.org, n=2, 
                          T1=1, a=1, b=2,
                          x.step = .02,
                          x.range.scale = 5, 
                          v.opt.mat=NULL, varphi.list=NULL,
                          method.spline = "fmm",
                          csv.name="test", save.ind=FALSE){
  #MC.size = 2e3, 
  #"fmm" may not be good at extrapolation; but I checked that it works for x^n
  varphi.list <- vector("list", n+1)
  varphi.list[[1]] <- varphi.org
  x.min <- -x.range.scale*b; x.max <- x.range.scale*b; 
  x.seq <- seq(x.min, x.max, x.step)
  y.seq <- numeric(length(x.seq))
  varphi.mat<- matrix(NA,nrow = length(x.seq), ncol = n+1)
  varphi.mat[,1] <- varphi.list[[1]](x.seq) 
  #save each v.opt
  v.opt.mat <- matrix(NA, nrow = (n-1)+2, ncol = length(x.seq)+2) 
  #also save the optimal v for the last x and last iteration, 
  #although it may not be very useful
  v.opt.mat[,1] <- b
  for (i in seq_len(n)){
    #the iteration step i
    temp <- varphi.list[[i]]
    for (j in seq_along(x.seq)){
      #the jth value of x.seq 
      v.init <- v.opt.mat[i,j]
      x1 <- x.seq[j]
      #define the expt of N(x,v^2); x is fixed and v is the arg.
      expt.v.fn <- function(v){
        f <- function(m) temp(x1+m)*(sqrt(n)/v)*dnorm(sqrt(n)*m/v)
        -integrate(Vectorize(f), -Inf, Inf)[[1]]
      }
      opt.list <- optim(v.init, expt.v.fn, method = "L-BFGS-B", 
                        lower = a, upper = b)
      y.seq[j] <- opt.list$value * (-1)
      v.opt.mat[i,j]<- v.opt.mat[i,(j+1)] <- opt.list$par 
      #update the v.opt at (i,j), guess the same for (i,j+1)
    }
    varphi.mat[,(i+1)] <- y.seq
    varphi.list[[i+1]]<- splinefun(x.seq, y.seq, 
                                   method = method.spline)
  }
  if(save.ind){
    #file.name <- paste0(csv.name, format(Sys.time(), "%Y-%m%d-%H%M%S"), ".csv")
    file.name <- paste0(csv.name, ".csv")
    write.csv(cbind(x.seq, varphi.mat), file = file.name) 
  }
  u.tx <- varphi.list[[n+1]]
  value <- varphi.list[[n+1]](0)
  return(list(value=value,
              varphi.mat = varphi.mat, 
              varphi.list = varphi.list,
              v.opt.mat = v.opt.mat,
              x.seq=x.seq,
              u.tx = u.tx))
}


#noise sequence
psemiGnorm.bisum <- function(q, set.sd = c(0.5,2)){
  sdl <- set.sd[1]
  sdr <- set.sd[2]
  #only consider two points
  varphi1 <- function(x){
    max(pnorm(q, mean = x, sd = set.sd))
  }
  
  fl <- function(m) dnorm(m)*varphi1(sdl*m)
  fr <- function(m) dnorm(m)*varphi1(sdr*m)
  
  l <- integrate(Vectorize(fl), -Inf, Inf)[[1]]
  r <- integrate(Vectorize(fr), -Inf, Inf)[[1]]
  
  max(l,r)
}
 

qsemiGnorm.bisum <- function(p, set.sd = c(0.5,2)){
  sdl <- set.sd[1]
  sdr <- set.sd[2]
  re <- uniroot(function(q) psemiGnorm.bisum(q, set.sd = set.sd)-p, 
                interval = c(-10,10)*sdr)
  re$root
}

VaR.semiGnorm.bisum <- function(alpha, set.sd = c(0.5,2)){
  -qsemiGnorm.bisum(alpha, set.sd = set.sd)
}

#VaR.semiGnorm.bisum(.01)

#estimation method based on point-interval transformation 

point.intvl.est <- function(y, n.guess, 
                            intvl.transform = function(x) x^2){
  #transform a point-valued noise into interval-valued 
  #then perform the estimation
  #intvl.transform: transformation on the interval
  #x^2 means the goal is to estimate the variance interval
  #(simpler interpretation)
  #(similar performance compared with max-mean)
  m <- floor(length(y)/n.guess)
  y <- y[1:(m*n.guess)]  
  y.mat <- matrix(y, ncol = n.guess, byrow = TRUE)
  y.intvl <- apply(y.mat, 2, range)
  
  y.mat <- matrix(y, ncol = n.guess, byrow = TRUE)
  y.cen <- apply(y.mat, 2, mean)
  y.cen.std <- scale(y.cen)
  #y.cen.std <- (y.cen - mean(y.cen))/sd(y.cen)
  # mean(w.cen.std) #almost zero
  # check whether we remove the mean part or not
  # consider the randomness of the center part itself; or remove the bias of the center part itself
  #plotl(w.cen.std)
  sd.set <- apply(y.mat, 1, sd)
  sd.intvl <- c(min(sd.set), max(sd.set))
  e.intvl <- cbind(y.cen.std,y.cen.std) 
  v.intvl <- rconstant.intvl(nrow(e.intvl), sd.intvl)
  y.intvl <- v.intvl*e.intvl
  #plot.intvl(w.intvl) #one version of the interval-valued data
  y.intvl.sq <- apply.intvl0(y.intvl, function(x) x^2)
  return(list(y.intvl=y.intvl, est.intvl = apply(y.intvl.sq, 2, mean)))
}

# test of model structure ----

## choices of varphi ----
#they should be non-convex and non-concave
#they should be symmetric around origin
exponent <- function(x, pow) (abs(x)^pow)*sign(x)
#for x^(1/3)
f_x3 <- function(x) x^3
f_sinx <- function(x) sin(x) 
#this may not be a good choice due to its seasonality, so we need to do standardization
f_sigmoid <- function(x) sigmoid(x) - 0.5 

#mean(f_sigmoid(rnorm(1e5)))

## main functions ----

sd.test.f1 <- function(sd.intvl){
  #written for the third moment
  sdl <- sd.intvl[1]
  sdr <- sd.intvl[2]
  # sdl.test <- sqrt(8*15*sdl^6)
  # sdr.test <- sqrt(8*15*sdr^6)
  sdl.test <- sqrt(15*sdl^6)
  sdr.test <- sqrt(15*sdr^6)
  c(sdl.test, sdr.test)
}

test.modelstructure <- function(x,
                                sdl.test = NULL, 
                                sdr.test = NULL, 
                                gr.size = 2, 
                                n.seq = seq(10, floor(length(x)/10), 5),
                                groupfn = mean, 
                                varphi = function(x) x^3,
                                sd.test.f = sd.test.f1, 
                                sig.level = 0.05, 
                                ovlp.ind = TRUE, 
                                step = 1,
                                return.ind = TRUE, 
                                print.ind = FALSE, 
                                print.intvl = FALSE, 
                                asym.stats = "semi-G-normal", 
                                est.sd.ind = TRUE, 
                                std.ind = FALSE){
  #x: input dataset 
  #std.ind: whether to standardize the dataset
  #sig.level: significance level
  #sdl.test, sdr.test: the lower and upper sd of the test stats
  #n.seq: the sequence of group size to be chosen from
  #which depend on varphi
  #asym.stats: the choice of the asymtotic distribution used for the test stats
  #pSemiGnorm or pGnorm
  #est.sd.ind: if TRUE, use original input dataset to estimate sd.intvl 
  #then use sd.test.f to compute sd.test (recommended for the current version), 
  #otherwise, directly use the transformed data (B.seq) to estimate sd.test
  
  N <- length(x)
  n <- gr.size
  if(std.ind){
    x <- (x-mean(x))/sd(x)
  }
  if(ovlp.ind){
    A <- get.groupfn.seq(x, n.guess = gr.size, step = step, groupfn = groupfn)
    # A <- get.meanseq.ovlp(x, n.guess = gr.size, step = step)
    #the two lines above are equivalent (checked)
    #A <- x[-1] + x[-N]
    #consider different order of dependency
    #use the function get.meanseq.ovlp
    #it will hold once we introduce the m-dependent version of semi-G-CLT
  } else {
    A <- get.groupfn.seq(x, n.guess = gr.size, step = step, 
                          groupfn = groupfn, 
                          novlp.ind = TRUE)
    # A <- get.meanseq.novlp(x, n.guess = gr.size)
    #k <- seq_len(N/2)
    #A <- x[2*k-1] + x[2*k]
  }
  B <- varphi(sqrt(n)*A)
  #check the sequence B
  #then think about estimate sd.test from B 
  m <- length(B)
  test.stat <- sqrt(m)*mean(B)
  #test.stat
  # sdl.test <- sqrt(8*15*sdl^6)
  # sdr.test <- sqrt(8*15*sdr^6)
  is.sd.test.null <- is.null(sdl.test) & is.null(sdr.test)
  if(is.sd.test.null){
    if (est.sd.ind){
      #use the data-driven rule to choose group size
      #directly use x sequence
      re <- choose.gr.size.onepath(x, n.seq = n.seq,
                                   f.trans = function(x) x^2,
                                   plot.ind = FALSE,
                                   print.re.ind = FALSE)
      #choose.gr.size.onepath is used to choose group size for var.intvl estimation
      #here we still use x^2 for numerical purpose
      var.est <- max.mean(x^2, n.guess = re$n.opt)
      sd.est <- sqrt(var.est)
      sd.test.est <- sd.test.f(sd.est)
      sdl.test <- sd.test.est[1]
      sdr.test <- sd.test.est[2]
    } else {
      #directly use B sequence 
      #test whether the B sequence has zero mean
      # re <- choose.gr.size.onepath(B, n.seq = n.seq, 
      #                              f.trans = function(x) x, 
      #                              plot.ind = FALSE, 
      #                              print.re.ind = FALSE)
      #use A sequence to choose group size 
      A1 <- sqrt(n)*A
      re <- choose.gr.size.onepath(A1, n.seq = n.seq,
                                   f.trans = function(x) x^2,
                                   plot.ind = FALSE,
                                   print.re.ind = FALSE)
      var.est <- max.mean(B^2, n.guess = re$n.opt)
      sd.test.est <- sqrt(var.est)
      #sd.test.est <- sd.test.f(sd.est)
      sdl.test <- sd.test.est[1]
      sdr.test <- sd.test.est[2]
      sd.est <- NULL
      #if directly use B sequence to estimate the sd.test 
      #the choice of varphi is important, otherwise
      #the estimation error will be large (especially when it involves high-order moments)
    }
    
  }
  
  if (asym.stats == "G-normal"){
    pf.asym <- pGnorm
    V.value <- 2*(1-pf.asym(abs(test.stat), upper = FALSE, 
                            sdl = sdl.test, sdr = sdr.test))
    v.value <- 2*(1-pf.asym(abs(test.stat), upper = TRUE, 
                            sdl = sdl.test, sdr = sdr.test))
    #we use the one sided tail prob
    #to approximate for G-normal two-sided tail prob 
    #(as mentioned in Peng and Zhou)
    #but we need to be careful about the small value of test stat
    #improve the current G-normal tail probability
    V.value <- truncate.f(V.value, 0, 1)
    v.value <- truncate.f(v.value, 0, 1)
  } else if (asym.stats == "semi-G-normal"){
    pf.asym <- pSemiGnorm
    V.value <- 2*(1-pf.asym(abs(test.stat), upper = FALSE, 
                            sdl = sdl.test, sdr = sdr.test))
    v.value <- 2*(1-pf.asym(abs(test.stat), upper = TRUE, 
                            sdl = sdl.test, sdr = sdr.test))
    #for semi-G-normal
    #the expression above is accurate
  } else {
    stop("Please provide an applicable input for asym.stats: 
         \"G-normal\" or \"semi-G-normal\".")
  }
   
  v.value1 <- format(v.value, digits = 3)
  V.value1 <- format(V.value, digits = 3)
  if(print.ind){
    if (print.intvl){
      cat("Under all possible situations in H0, \n the p-value is in the interval",
          paste0("[", v.value1, ",", V.value1, "],"),"\n at significance level", 
          sig.level, ".\n")
    } else {
      cat("Under H0, the V-value is",
         V.value1,", at sig.level =", sig.level, ", we", 
         ifelse(V.value<sig.level, "reject", "accept"), "H0.\n")
    }
    
  }
  
  #The V value is
  #under all possible situations in H0, the p-value is at most
  #we can also provide the pair of p-value if needed
  if(return.ind){
    if(is.sd.test.null){
      return(list(test.stat = test.stat, 
                  upper.p.value=V.value, 
                  lower.p.value=v.value, 
                  p.value.intvl = c(v.value, V.value),
                  sd.test.est = sd.test.est, 
                  sd.est = sd.est, 
                  A.seq = A,
                  B.seq = B))
    } else {
      return(list(test.stat = test.stat, 
                  upper.p.value=V.value, 
                  lower.p.value=v.value, 
                  p.value.intvl = c(v.value, V.value)))
    }
  }
  
  
}

test.ms <- test.modelstructure
#short form of test.modelstructure

test.ms.sim <- function(x.gen = NULL, 
                        x.len = 5e3,
                        sdl.test = NULL, sdr.test = NULL,
                        gr.size = 2, 
                        step = 2, 
                        groupfn = mean, 
                        varphi = function(x) x^3,
                        sig.level = 0.05,
                        repeat.times = 1e3,
                        plot.ind = TRUE,
                        print.ind = TRUE,
                        intvl.plot = FALSE, 
                        return.ind = FALSE, 
                        main.model = NULL, 
                        est.sd.ind = TRUE, 
                        std.ind = FALSE){
  #this is a function used for simulation study
  #x.gen: a function of n to generate x, e.g. x.gen = function(n) rnorm(n)
  #plot.ind: whether to plot the result if repeat.ind is TRUE
  #x.len: the length of simulated data
  #main: the title of the plot (the name of the generation scheme)
  #...: arguments to be passed to test.ms
  V.val.mat <- replicate(repeat.times, {
    x <- x.gen(x.len)
    re <- test.ms(x, std.ind = std.ind, gr.size = gr.size, step = step,
                  sdl.test = sdl.test, sdr.test = sdr.test, 
                  groupfn = groupfn, varphi = varphi, 
                  sig.level = sig.level, 
                  print.ind = FALSE, est.sd.ind = est.sd.ind)
    re$p.value.intvl
  })
  V.val.mat <- t(V.val.mat)
  V.val.seq <- V.val.mat[,2]
  #if(is.null(main))
  if(plot.ind){
    hist(V.val.seq, main = paste("Histogram of V-values",
                                 main.model))
    if(intvl.plot){
      re <- plot.intvl(V.val.mat, 
                       main = paste("Plot of the lower (in blue) 
                                  and upper (in red) p-values with", main),
                       print.ind = FALSE, return.ind = TRUE)
      re2 <- re + geom_hline(yintercept = sig.level, colour = "brown", 
                             linetype = 'dashed')
      print(re2)
    }
  }
  
  if(print.ind){
    cat("The rejection rate = ", mean(V.val.seq < sig.level), "\n", 
        "The acceptance rate = ", mean(V.val.seq >= sig.level), "\n") 
  }
  
  if(return.ind){
    return(list(pvalue.lowup.mat=V.val.mat, 
                V.val.seq =  V.val.seq, 
                reject.rate = mean(V.val.seq < sig.level)))
  }
}

truncate.f <- function(x, a, b){
  #truncate x into an interval [a,b]
  a*(x<=a)+x*(a<x&x<=b)+b*(x>=b)
}

# f <- function(a,b) {a+b}
# g <- function(x, ...){
#   #pass ... to f
#   x + f(...)
# }

## noise generation ---- 

garch.gjr <- function(N, return.full = FALSE, 
                      phi = 0.1754, 
                      a0 = 0.0154, 
                      a1 = 0, 
                      b1 = 0.8974, 
                      ind.trunc = TRUE,
                      sdl = 1, 
                      sdr = 2){
  ep <- rnorm(N)
  v <- rep(sdl, N)
  y <- v*ep
  if(ind.trunc){
    for (i in seq_len(N-1)){
      v2 <- a0 + a1*y[i]^2 + b1*v[i]^2 + phi*y[i]^2*(y[i]>=0)
      v[i+1] <- truncate.f(sqrt(v2), sdl, sdr)
      y[i+1] <- v[i+1]*ep[i+1]
    }
  } else {
    for (i in seq_len(N-1)){
      v2 <- a0 + a1*y[i]^2 + b1*v[i]^2 + phi*y[i]^2*(y[i]>=0)
      v[i+1] <- sqrt(v2)
      y[i+1] <- v[i+1]*ep[i+1]
    }
  }
  
  if(return.full){
    return(list(sig=v, y=y))
  } else {
    return(y)
  }
}

garch.gjr2 <- function(N, return.full = FALSE, 
                      phi = 0.1754, 
                      a0 = 0.0154, 
                      a1 = 0, 
                      b1 = 0.8974, 
                      ind.trunc = TRUE,
                      sdl = 1, 
                      sdr = 2){
  ep <- rnorm(N)
  v <- rep(sdl, N)
  y <- v*ep
  if(ind.trunc){
    for (i in seq_len(N-1)){
      v2 <- a0 + a1*y[i]^2 + b1*v[i]^2 + phi*y[i]^2*(y[i]<0)
      v[i+1] <- truncate.f(sqrt(v2), sdl, sdr)
      y[i+1] <- v[i+1]*ep[i+1]
    }
  } else {
    for (i in seq_len(N-1)){
      v2 <- a0 + a1*y[i]^2 + b1*v[i]^2 + phi*y[i]^2*(y[i]<0)
      v[i+1] <- sqrt(v2)
      y[i+1] <- v[i+1]*ep[i+1]
    }
  }
  
  if(return.full){
    return(list(sig=v, y=y))
  } else {
    return(y)
  }
}


rthres.lag <- function(N, order.lag = 2){
  ep <- rnorm(N)
  y <- v <- numeric(N)
  k <- order.lag
  v[seq_len(k)] <- sample(c(sdl,sdr), k, replace = TRUE)
  y[seq_len(k)] <- v[seq_len(k)]*ep[seq_len(k)]
  #sw.loc <- seq(0, sdr*2, length.out = N-2)
  for (i in seq_len(N-k)){
    v[i+k] <- ifelse(y[i]>0, sdl, sdr)
    #v[i+k] <- ifelse(y[i]<0, sdl, sdr)
    #v[i+2] <- ifelse(y[i]>runif(1,-sdr,sdr), sdr, sdl)
    #v[i+2] <- ifelse(y[i]>sw.loc[i], sdr, sdl)
    y[i+k] <- v[i+k]*ep[i+k]
  }
  #plotl(v)
  y
  #for this structure
  #we should accept the null hypothesis for group size = 2
  #but reject the null hypothesis for group size > 2
}



