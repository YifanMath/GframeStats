#organize existing code
library(numDeriv)
library(MASS)

bias.SE <- function(x, par.true){
  bias <- mean(x - par.true)
  SE <- sd(x - par.true)
  re <- c(bias, SE)
  names(re) <- c("bias", "SE")
  re
}

#linear regression
#max-mean loss function with regularization

#try it if the sub.lm1 does not work well
f.gr <- function(theta, gr.size){
  #write the loss function 
  
}

#Sublinear Least Abosolute 
#change it to L2 if needed 

# lm.cl <- function(dat){
#   mod <- lm(y~x, data = dat)
#   
# }

lm.sub <- function(y.name, x.name, data, gr.size, initialPar = c(0,1), type = "L1"){
  #y is the column name of the data for response
  #x is the column names of the data for covariates
  #group the data
  #type: type of the regression, chosen from c("L1","L2")
  if(type == "L1"){
    loss.f <- function(x) abs(x)
  } else if (type == "L2"){
    loss.f <- function(x) x^2
  } else {
    stop("The type of regression has to be specified. The current option is L1 or L2.")
  }
  n <- gr.size 
  y <- data[,y.name]
  x <- data[,x.name]
  #equal group size
  m <- floor(length(y)/n) # num of groups
  N <- n*m #the length of applicable rows
  dat <- data[seq_len(N), ]
  #y.mat <- matrix(y, ncol = n, byrow = TRUE)
  ind <- seq_len(N)
  ind.mat <- matrix(ind, ncol = n, byrow = TRUE)
  H <- function(beta){
    beta.t <- matrix(beta, ncol = 1)
    SSE.seq <- numeric(m)
    for (i in seq_len(m)){
      group.dat <- dat[ind.mat[i,], ]
      # y.gr <- as.matrix(group.dat[,y.name])
      # x.gr <- as.matrix(group.dat[,x.name])
      y.gr <- group.dat[,y.name]
      x.gr <- group.dat[,x.name]
      SSE.seq[i] <- mean(loss.f(y.gr - cbind(1,x.gr) %*% beta.t))
      #SSE.seq[i] <- mean((y.gr - cbind(1,x.gr) %*% beta.t)^2)
    }
    max(SSE.seq)
  }
  re <- optim(initialPar, H) 
  #improve this optimization by adding a regularization term 
  beta.G <- re$par
  x.names <- c("Intercept", x.name)
  names(beta.G) <- x.names
  re <- beta.G #to be extended later
  return(re)
}

yw.cl <- function(x){
  re <- ar.yw(x, order.max = 1, aic = FALSE, demean = FALSE)
  re$ar
}

yw.G <- function(x, n.guess){
  #improved version
  #to be compared with the classical yw
  n <- n.guess
  N <- length(x)
  x2.0 <- x^2
  x2.1 <- x[-N]*x[-1]
  g0.itvl <- max.mean(x2.0, n)
  g1.itvl <- max.mean(x2.1, n)
  sum(g1.itvl)/sum(g0.itvl)
}



lm.sub1 <- function(y.name, x.name, data, gr.size, initialPar, level=.95, compareLS=FALSE, 
                    type = "L1"){
  #y is the column name of the data for response
  #x is the column names of the data for covariates
  #x <- cbind(1, x1, ..., xp)
  #n is the group size. 
  #group the data
  #type: type of the regression, chosen from c("L1","L2")
  if(type == "L1"){
    loss.f <- function(x) abs(x)
  } else if (type == "L2"){
    loss.f <- function(x) x^2
  } else {
    stop("The type of regression has to be specified. The current option is L1 or L2.")
  }
  
  n <- gr.size 
  y <- data[,y.name]
  x <- data[,x.name]
  #equal group size
  m <- floor(length(y)/n) # num of groups
  N <- n*m #the length of applicable rows
  dat <- data[seq_len(N), ]
  #y.mat <- matrix(y, ncol = n, byrow = TRUE)
  ind <- seq_len(N)
  ind.mat <- matrix(ind, ncol = n, byrow = TRUE)
  H <- function(beta){
    beta.t <- matrix(beta, ncol = 1)
    SSE.seq <- numeric(m)
    for (i in seq_len(m)){
      group.dat <- dat[ind.mat[i,], ]
      # y.gr <- as.matrix(group.dat[,y.name])
      # x.gr <- as.matrix(group.dat[,x.name])
      y.gr <- group.dat[,y.name]
      x.gr <- group.dat[,x.name]
      #SSE.seq[i] <- mean(loss.f(y.gr - cbind(1,x.gr) %*% beta.t))
      SSE.seq[i] <- mean((y.gr - cbind(1,x.gr) %*% beta.t)^2)
    }
    max(SSE.seq)
  }
  re <- optim(initialPar, H) 
  #improve this optimization by adding a regularization term 
  beta.G <- re$par
  
  #ref: regression under sublinear expectation
  #compute the sig.max
  
  beta.G.t <- matrix(beta.G, ncol = 1)
  sig.seq <- numeric(m)
  for (i in seq_len(m)){
    group.dat <- dat[ind.mat[i,], ]
    y.gr <- group.dat[,y.name]
    x.gr <- group.dat[,x.name]
    sig.seq[i] <- var(y.gr - cbind(1,x.gr) %*% beta.G.t)
  }
  sig.max <- max(sig.seq)
  xx.inv <- solve(t(x)%*%x)
  #z <- abs(qnorm((1-level)/2)) #when level=0.95, z=1.96
  #sd.error.seq <- sqrt(sig.max*diag(xx.inv))
  #CI.wid.G <- t(t(sd.error.seq))%*%t(c(-z,z))/sqrt(n)
  #C.I.mat <- CI.wid.G + cbind(beta.G, beta.G)
  x.names <- c("Intercept", x.name)
  names(beta.G) <- x.names
  #names(beta.ls) <- x.names
  if(compareLS){
    #if initialPar is the beta from least square
    formula.c <- paste(y.name,"~",x.name)
    mod.ls <- lm(as.formula(formula.c), data = data)
    beta.ls.t <- matrix(mod.ls$coefficients, ncol=1)
    sig.seq.ls <- numeric(m)
    for (i in seq_len(m)){
      group.dat <- dat[ind.mat[i,], ]
      y.gr <- group.dat[,1]
      x.gr <- group.dat[,-1]
      sig.seq.ls[i] <- var(y.gr - x.gr %*% beta.ls.t)
    }
    sig.ls <- var(y - x%*%beta.ls.t)
    sd.error.seq.ls <- sqrt(c(sig.ls)*diag(xx.inv))
    CI.wid.ls <- t(t(sd.error.seq.ls))%*%t(c(-z,z))/sqrt(N) #caution! this is N not n
    row.names(CI.wid.G) <- row.names(CI.wid.ls) <- x.names
    return(list(par.G=beta.G, par.ls=beta.ls,
                sig.G.max=sig.max, sig.ls.max=max(sig.seq.ls),
                sig.ls.global=c(unname(sig.ls)),
                mean.sig.G = mean(sig.seq),
                mean.sig.ls = mean(sig.seq.ls),
                CI.width.G=CI.wid.G,
                CI.width.ls=CI.wid.ls))
  } else {
    #names(sd.error.seq) <- x.names
    #rownames(C.I.mat) <- x.names
    return(list(par=beta.G, sig.max=sig.max))
  }
}

sub.lm1 <- function(y, x, gr.size, initialPar, level=.95, compareLS=FALSE, 
                    type = "L1"){
  #y is the vector of data for response variable.
  #x is the matrix of covariates (predictors).
  #x <- cbind(1, x1, ..., xp)
  #n is the group size. 
  #group the data
  #type: type of the regression, chosen from c("L1","L2")
  if(type == "L1"){
    loss.f <- function(x) abs(x)
  } else if (type == "L2"){
    loss.f <- function(x) x^2
  } else {
    stop("The type of regression has to be specified. The current option is L1 or L2.")
  }
  n <- gr.size 
  x <- t(t(x))
  m <- floor(length(y)/n)
  N <- n*m
  y <- y[seq_len(N)]
  x <- t(t(x[seq_len(N),]))
  dat <- cbind(y, x)
  #get the names
  if (all(x[,1]==1)){
    p <- ncol(x)-1
    x.names <- c("(Intercept)",paste0("x", 1:p)) #paste0:no space
  } else {
    p <- ncol(x)
    x.names <- c(paste0("x", 1:p)) #paste0:no space
  }
  #y.mat <- matrix(y, ncol = n, byrow = TRUE)
  ind <- seq_len(N)
  ind.mat <- matrix(ind, ncol = n, byrow = TRUE)
  H <- function(beta){
    beta.t <- matrix(beta, ncol = 1)
    SSE.seq <- numeric(m)
    for (i in seq_len(m)){
      group.dat <- dat[ind.mat[i,], ]
      y.gr <- group.dat[,1]
      x.gr <- t(t(group.dat[,-1]))
      SSE.seq[i] <- mean(abs(y.gr - cbind(1,x.gr) %*% beta.t))
    }
    max(SSE.seq)
  }
  re <- optim(initialPar, H) 
  #improve this optimization by adding a regularization term 
  beta.G <- re$par
  #compute the sig.max
  beta.G.t <- matrix(beta.G, ncol = 1)
  sig.seq <- numeric(m)
  for (i in seq_len(m)){
    group.dat <- dat[ind.mat[i,], ]
    y.gr <- group.dat[,1]
    x.gr <- group.dat[,-1]
    sig.seq[i] <- var(y.gr - x.gr %*% beta.G.t)
  }
  sig.max <- max(sig.seq)
  xx.inv <- solve(t(x)%*%x)
  #z <- abs(qnorm((1-level)/2)) #when level=0.95, z=1.96
  #sd.error.seq <- sqrt(sig.max*diag(xx.inv))
  #CI.wid.G <- t(t(sd.error.seq))%*%t(c(-z,z))/sqrt(n)
  #C.I.mat <- CI.wid.G + cbind(beta.G, beta.G)
  names(beta.G) <- x.names
  #names(beta.ls) <- x.names
  if(compareLS){
    #if initialPar is the beta from least square
    beta.ls.t <- matrix(initialPar, ncol = 1)
    sig.seq.ls <- numeric(m)
    for (i in seq_len(m)){
      group.dat <- dat[ind.mat[i,], ]
      y.gr <- group.dat[,1]
      x.gr <- group.dat[,-1]
      sig.seq.ls[i] <- var(y.gr - x.gr %*% beta.ls.t)
    }
    sig.ls <- var(y - x%*%beta.ls.t)
    sd.error.seq.ls <- sqrt(c(sig.ls)*diag(xx.inv))
    CI.wid.ls <- t(t(sd.error.seq.ls))%*%t(c(-z,z))/sqrt(N) #caution! this is N not n
    row.names(CI.wid.G) <- row.names(CI.wid.ls) <- x.names
    return(list(par.G=beta.G, par.ls=beta.ls,
                sig.G.max=sig.max, sig.ls.max=max(sig.seq.ls),
                sig.ls.global=c(unname(sig.ls)),
                mean.sig.G = mean(sig.seq),
                mean.sig.ls = mean(sig.seq.ls),
                CI.width.G=CI.wid.G,
                CI.width.ls=CI.wid.ls))
  } else {
    #names(sd.error.seq) <- x.names
    #rownames(C.I.mat) <- x.names
    return(list(par=beta.G, sig.max=sig.max))
  }
}

boot.mov.block <- function(N.boot, block.size, step, dat.seq){
  #generate the index for moving-block bootstrap
  N <- length(dat.seq)
  b <- block.size
  m <- ceiling(N.boot/b) # the number of blocks
  head.ind <- seq(1, N-b+1, step)
  head.ind.sub <- sample(head.ind, m, replace = TRUE)
  ind.mat <- sapply(head.ind.sub, function(x) x+seq_len(b)-1)
  re.seq <- numeric(m*b)
  for (j in seq_len(m)){
    #seq_len vs seq_along
    re.seq[(j-1)*b + seq_len(b)] <- dat.seq[ind.mat[,j]]
  }
  boot.seq <- re.seq[seq_len(N.boot)]
  #return(list(ind.boot.mat = t(ind.mat), 
  #            boot.seq = boot.seq))
  return(boot.seq)
}

#also check the setup for Yule walker estimator 

noise.gen.eqblock <- function(N, gr.size, 
                              sd.set = c(0.5, 1, 2), 
                              prob.set = c(3,2,1)/6,
                              z.ind = NULL){
  #z.ind should be consistent with sd.set
  m <- ceiling(N/gr.size)
  if (is.null(z.ind)){
    z.val <- sample(sd.set, m, replace = TRUE, prob = prob.set)
  } else {
    z.val <- sd.set[z.ind]
  }
  z.seq <- rep(z.val, each = gr.size)
  w.seq <- z.seq[seq_len(N)]*rnorm(N)
  return(w.seq)
}


ar1.gen <- function(noise, par){
  #generate AR(1)
  N <- length(noise)
  x.seq <- numeric(N-1)
  x.seq[1] <- noise[1]
  for (t in seq_len(N-1)){
    x.seq[t+1] <- x.seq[t]*par + noise[t+1]
  }
  return(x.seq)
}

ma1.gen <- function(noise, par){
  #generate AR(1)
  N <- length(noise)
  x.seq <- numeric(N)
  x.seq[1] <- noise[1]
  #or use vectorization
  for (t in seq_len(N-1)){
    x.seq[t+1] <- noise[t]*par + noise[t+1]
  }
  return(x.seq)
}

ma.gen <- function(noise, par){
  # par = c(theta1, theta2, ..., thetap)
  p <- length(par) #the order of MA
  N <- length(noise) 
  #the length of the time series depends on the order p
  #create the noise 
  noise.mat1 <- sapply(seq_len(p+1), function(start.ind) {
    ind.seq <- seq(from = start.ind, length.out = N-p, by = 1)
    noise[ind.seq]
    })
  noise.mat <- noise.mat1[,rev(seq_len(p+1))]
  par1 <- c(1,par)
  x.seq <- apply(noise.mat, 1, function(x) sum(x*par1))
  return(x.seq)
}
#to check
# noise1 <- 1:10
# par <- c(1,2)
# ma.gen(noise1, par)

ar1.gen.sigdepn <- function(N, gr.size = 2e2, sdl=0.5, sdr=2, b.true=0.7, 
                            plot.ind = FALSE, a = NULL){
  x.seq <- w.seq <- sig.seq <- numeric(N)
  e.seq <- rnorm(N)
  #intermediate vol level
  sdm <- sdl+(sdr-sdl)/3
  #a <- sqrt(sdm^2/(1-b.true^2))*1.5
  if (is.null(a)){
    a <- sdl
  }
  m <- N/gr.size
  #the index of switching location 
  t.sw <- seq(1, N, gr.size)
  sig.val <- numeric(m)
  sig.val[1] <- sample(c(sdl,sdm,sdr),1)
  sig.seq[seq_len(gr.size)] <- rep(sig.val[1],gr.size)
  w.seq[seq_len(gr.size)] <- sig.seq[seq_len(gr.size)]*e.seq[seq_len(gr.size)]
  x.seq[1] <- w.seq[1]
  for (t in seq_len(gr.size-1)){
    x.seq[t+1] <- x.seq[t]*b.true + w.seq[t+1]
  }
  #generate w.seq
  for (i in seq_len(m-1)){
    #generate sig.seq and w.seq in block (i+1)
    x.incre <- x.seq[i*gr.size]-x.seq[(i-1)*gr.size+1]
    #add a markov dependence if needed
    #use the model with leverage effect 
    #to create more deviation to the normalized sum
    sig.val[i+1] <- (x.incre<=-a)*sdm + (-a<x.incre&x.incre<=a)*sdl + (x.incre>a)*sdr
    ind.new <- i*gr.size+seq_len(gr.size)
    sig.seq[ind.new] <- rep(sig.val[i+1],gr.size)
    w.seq[ind.new] <- sig.seq[ind.new]*e.seq[ind.new]
    x.seq[ind.new[1]] <- w.seq[ind.new[1]]
    #ind.new2 <- i*gr.size+seq_len(gr.size-1)
    for (t in seq_len(gr.size-1)){
      #generate x.seq in block (i+1)
      x.seq[ind.new[t+1]] <- x.seq[ind.new[t]]*b.true + w.seq[ind.new[t+1]]
    }
  }
  if(plot.ind){
    t.seq <- seq_len(N)
    plot(t.seq, sig.seq, type = "l")
    plot(t.seq, w.seq, type = "l")
    plot(t.seq, x.seq, type = "l")
  }
  return(list(sig.val=sig.val, sig.seq=sig.seq, 
              w.seq=w.seq, x.seq=x.seq))
  #make an external variable if needed
  #make it hard to be modeled by a single classical SV model
  #we need to have a short bad period (the proportion should not be large)
  #consider each group as one period 
}


f.bounded <- function(x, min, max){
  min*(x<=min)+x*(min<x & x<=max)+x*(x>max)
}

ar1.gen.sigdepn.noblock <- function(N, sdl=0.5, sdr=2, 
                                    b.true=0.7, 
                                    rho = -0.6, sig = 1.5, phi = 0.5, b.sc = 0.1,
                                    plot.ind = FALSE, a = NULL, bounded.var = FALSE){
  #(rho, sig, phi, b.sc) 
  #these are the parameters for SV 
  #without leverage: rho=0
  #with leverage: rho is negative 
  #intermediate vol level
  #ep and eta
  Sigma1 <- matrix(c(1, rho*sig,rho*sig, sig^2), ncol = 2)
  ep.eta.mat <- mvrnorm(N, mu = c(0,0), Sigma = Sigma1)
  e.seq <- ep.eta.mat[,1]
  eta.seq <- ep.eta.mat[,2]
  w.seq <- x.seq <- g.seq <- numeric(N)
  for (t in seq_len(N-1)){
    g.seq[t+1] <- phi*g.seq[t] + eta.seq[t]
  }
  #check whether we need bounded variance here
  if (bounded.var) {
    sig.seq <- f.bounded(exp(g.seq/2)*b.sc, min=sdl, max=sdr)
  } else {
    sig.seq <- exp(g.seq/2)*b.sc
  }
  w.seq <- sig.seq*e.seq
  x.seq <- ar1.gen(w.seq, par = b.true)
  
  if(plot.ind){
    t.seq <- seq_len(N)
    plot(t.seq, sig.seq, type = "l")
    plot(t.seq, w.seq, type = "l")
    plot(t.seq, x.seq, type = "l")
  }
  return(list(sig.seq=sig.seq, 
              w.seq=w.seq, x.seq=x.seq))
  #make an external variable if needed
  #make it hard to be modeled by a single classical SV model
  #we need to have a short bad period (the proportion should not be large)
  #consider each group as one period 
}

#traditional Yule-walker method

#compare the performance of estimation

#check coverage
cover.check <- function(CI, par.true){
  k <- nrow(CI)
  cover.ind <- logical(k)
  for (i in seq_len(k)){
    cover.ind[i] <- (CI[i,1]<=par.true[i]) & (par.true[i] <= CI[i,2])
  }
  names(cover.ind) <- rownames(CI)
  cover.ind
}

CI.cover.comp <- function(N, ar1.gen, name.char, 
                          save.plot = TRUE, 
                          print.plot = FALSE, 
                          rep.time=2e3, 
                          boot.time = 1e3,
                          size.coverage = 1e3, 
                          alpha = 0.05, 
                          b.true = 0.7, 
                          gr.size = NULL, 
                          n.min = 10, n.max.prop=1/2, n.step=5, 
                          plot.n.choose = FALSE){
  #if gr.size = NULL, we use the newly developed criterion based on BIC to decide
  #n.max.prop*N should be an integer
  #n.seq = seq(n.min, n.max.prop*N, n.step) gives the sequence of group sizes to check
  re <- ar1.gen(N)
  x.seq <- re$x.seq
  dat.tx <- data.frame(t = seq_len(N), x = x.seq)
  p.x <- ggplot(data = dat.tx, aes(x=t, y=x)) + geom_line()
  if (is.null(gr.size)) {
    cover.mat <- replicate(rep.time, {
      re <- ar1.gen(N)
      x.seq <- re$x.seq
      y <- x.seq[-1]
      x <- x.seq[-N]
      dat.yx <- data.frame(y=y, x=x)
      #mod.cl.L2 <- lm(y~0+x, data=dat.yx)
      #CI.cl.L2 <- confint(mod.cl.L2, level = 0.95)
      #b.cl.L2 <- mod.cl.L2$coefficients
      b.cl.yw <- mean(y*x)/mean(x^2)
      b.cl <- b.cl.yw
      re <- acf(x.seq,1, type = "covariance", plot = FALSE)
      acov.seq <- as.matrix(re$acf)
      sig2.hat <- acov.seq[1,] - b.cl*acov.seq[2,]
      sig.hat <- sqrt(sig2.hat)
      #autocovariance matrix
      acov.mat <- cbind(acov.seq, rev(acov.seq))
      acov.mat.inv <- solve(acov.mat)
      acov.mat.inv.sqrt <- sqrtm(acov.mat.inv)
      ME.cl.L2 <- qnorm(1-alpha/2)*(sig.hat*acov.mat.inv.sqrt[1,1])/sqrt(N)
      CI.cl.L2 <- b.cl + c(-1,1)*ME.cl.L2
      cover.cl.L2 <- (CI.cl.L2[1] <= b.true) & (b.true <= CI.cl.L2[2])
      
      #write this procedure into a function
      res.seq <- y-b.cl*x
      #var.seq.est <- get.meanseq.ovlp(res.seq^2, n.guess = gr.size, step = 50)
      #var.seq.est <- get.meanseq.novlp(res.seq^2, n.guess = gr.size)
      #here gr.size = NULL, use a data-driven way as an adaptive way to decide it 
      n.opt.re <- choose.gr.size.onepath(res.seq, f.obj = f.BIC, 
                                           f.trans = function(x) x, 
                                           n.seq = seq(n.min, n.max.prop*N, n.step), 
                                         plot.ind = plot.n.choose)
      gr.size <- n.opt.re$n.opt
      var.intvl.est <- max.mean(res.seq^2, n.guess = gr.size)
      #var.intvl.est <- quantile(var.seq.est, probs = c(0.25, 0.75))
      sd.intvl.est <- sqrt(var.intvl.est)
      #one heuristic rule to get less robust but sharper results
      #not using the maximum sigma
      #may be the 75% quantile
      #semi-G-normal under semi-sequential indep 
      beta1.ME.semiG <- qnorm(1-alpha/2)*(sd.intvl.est[2]*acov.mat.inv.sqrt[1,1])/sqrt(N)
      beta1.CI.semiG <- b.cl + c(-1,1)*beta1.ME.semiG
      #semi-G-normal under sequential indep, then treat it as G-normal as asymptotic
      sc.G <- (1+sd.intvl.est[1]/sd.intvl.est[2])/2
      beta1.ME.G <- qnorm(1-sc.G*alpha/2)*(sd.intvl.est[2]*acov.mat.inv.sqrt[1,1])/sqrt(N)
      beta1.CI.G <- b.cl+ c(-1,1)*beta1.ME.G
      #cover.cl.L2 <- cover.check(CI.cl.L2, par.true = beta.true)
      cover.semiG <- (beta1.CI.semiG[1] <= b.true) & (b.true <= beta1.CI.semiG[2])
      cover.G <- (beta1.CI.G[1] <= b.true) & (b.true <= beta1.CI.G[2])
      cover.re <- c(cover.cl.L2, cover.semiG, cover.G)
      names(cover.re) <- c("cover.cl.L2", "cover.semiG", "cover.G")
      cover.re
    })
  } else {
    cover.mat <- replicate(rep.time, {
      re <- ar1.gen(N)
      x.seq <- re$x.seq
      y <- x.seq[-1]
      x <- x.seq[-N]
      dat.yx <- data.frame(y=y, x=x)
      #mod.cl.L2 <- lm(y~0+x, data=dat.yx)
      #CI.cl.L2 <- confint(mod.cl.L2, level = 0.95)
      #b.cl.L2 <- mod.cl.L2$coefficients
      b.cl.yw <- mean(y*x)/mean(x^2)
      b.cl <- b.cl.yw
      re <- acf(x.seq,1, type = "covariance", plot = FALSE)
      acov.seq <- as.matrix(re$acf)
      sig2.hat <- acov.seq[1,] - b.cl*acov.seq[2,]
      sig.hat <- sqrt(sig2.hat)
      #autocovariance matrix
      acov.mat <- cbind(acov.seq, rev(acov.seq))
      acov.mat.inv <- solve(acov.mat)
      acov.mat.inv.sqrt <- sqrtm(acov.mat.inv)
      ME.cl.L2 <- qnorm(1-alpha/2)*(sig.hat*acov.mat.inv.sqrt[1,1])/sqrt(N)
      CI.cl.L2 <- b.cl + c(-1,1)*ME.cl.L2
      cover.cl.L2 <- (CI.cl.L2[1] <= b.true) & (b.true <= CI.cl.L2[2])
      
      #write this procedure into a function
      res.seq <- y-b.cl*x
      #var.seq.est <- get.meanseq.ovlp(res.seq^2, n.guess = gr.size, step = 50)
      #var.seq.est <- get.meanseq.novlp(res.seq^2, n.guess = gr.size)
      #here gr.size is given by user
      var.intvl.est <- max.mean(res.seq^2, n.guess = gr.size)
      #var.intvl.est <- quantile(var.seq.est, probs = c(0.25, 0.75))
      sd.intvl.est <- sqrt(var.intvl.est)
      #one heuristic rule to get less robust but sharper results
      #not using the maximum sigma
      #may be the 75% quantile
      #semi-G-normal under semi-sequential indep 
      beta1.ME.semiG <- qnorm(1-alpha/2)*(sd.intvl.est[2]*acov.mat.inv.sqrt[1,1])/sqrt(N)
      beta1.CI.semiG <- b.cl + c(-1,1)*beta1.ME.semiG
      #semi-G-normal under sequential indep, then treat it as G-normal as asymptotic
      sc.G <- (1+sd.intvl.est[1]/sd.intvl.est[2])/2
      beta1.ME.G <- qnorm(1-sc.G*alpha/2)*(sd.intvl.est[2]*acov.mat.inv.sqrt[1,1])/sqrt(N)
      beta1.CI.G <- b.cl+ c(-1,1)*beta1.ME.G
      #cover.cl.L2 <- cover.check(CI.cl.L2, par.true = beta.true)
      cover.semiG <- (beta1.CI.semiG[1] <= b.true) & (b.true <= beta1.CI.semiG[2])
      cover.G <- (beta1.CI.G[1] <= b.true) & (b.true <= beta1.CI.G[2])
      cover.re <- c(cover.cl.L2, cover.semiG, cover.G)
      names(cover.re) <- c("cover.cl.L2", "cover.semiG", "cover.G")
      cover.re
    })
  }
  
  cover.mat <- t(cover.mat)
  
  L <- nrow(cover.mat)
  cover.rate.mat <- replicate(boot.time, {
    boot.ind <- sample(seq_len(L), size.coverage, replace = TRUE)
    cover.mat.boot <- cover.mat[boot.ind, ]
    #cover.seq.boot <- sample(cover.seq, 2e3, replace = TRUE)
    apply(cover.mat.boot, 2, mean)
  })
  cover.rate.mat <- t(cover.rate.mat)
  
  #histogram
  #sharper results
  dat.short <- as.data.frame(cover.rate.mat)
  dat.short$ID <- seq_along(dat.short[,1])
  dat.long <- melt(dat.short, variable.name = "method", 
                   value.name = "cover.rate.est", id.vars = "ID")
  
  p1 <- ggplot(dat.long, aes(x=cover.rate.est, fill = method)) + geom_histogram(position = "identity", alpha = 0.5, bins = 25) + geom_vline(xintercept = 0.95, colour = "blue", linetype = "dashed")
  
  p2 <- ggplot(dat.long, aes(x=cover.rate.est, fill = method)) + geom_density(alpha = 0.5) + geom_vline(xintercept = 0.95, colour = "blue", linetype = "dashed")
  
  p.boxplot <- ggplot(dat.long, aes(y=cover.rate.est, fill = method)) + geom_boxplot() + geom_hline(yintercept = 0.95, colour = "blue", linetype = "dashed")
  
  if (print.plot){
    print(p.x)
    print(p1)
    print(p2)
    print(p.boxplot)
  }
  
  if (save.plot){
    ggsave(paste0("ar1-x-plot-",name.char,".pdf"), p.x, device = "pdf", width = 6, height = 4)
    ggsave(paste0("ar1-yw-hist-comp-",name.char,".pdf"), p1, device = "pdf", width = 6, height = 4)
    ggsave(paste0("ar1-yw-den-comp-",name.char,".pdf"), p2, device = "pdf", width = 6, height = 4)
    ggsave(paste0("ar1-yw-boxplot-comp-",name.char,".pdf"), p.boxplot, device = "pdf", width = 6, height = 4)
  }
  
  return(list(cover.rate.mat = cover.rate.mat, 
              p.x=p.x, 
              p.hist=p1, p.den=p2, p.boxplot=p.boxplot))
}



