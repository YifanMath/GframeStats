#last update date 2021-07-17
choose.gr.size.onepath <- function(ep.seq, n.seq, f.obj = f.BIC, 
                                   f.trans = function(x) x,
                                   plot.ind = TRUE,
                                   print.re.ind = TRUE,
                                   save.ind = FALSE, name.char = NULL, 
                                   width = 6, height = 4, 
                                   novlp.ind = TRUE, 
                                   step = 10){
  #n.seq: the sequence of candidate group sizes
  #f.obj: the objective function we can compute under each group size n 
  #rep.time: the times of replications to draw n.seq vs val.seq
  #novlp.ind: whether or not we use non-overlapping groups
  #if novlp.ind=FALSE, we use ovlp group, then we need to use argument step
  #step: the step size of the moving group (to produce ovlp groups)
  
  N <- length(ep.seq)
  val.seq <- numeric(length(n.seq))
  
  y.seq <- f.trans(ep.seq)
  
  if(novlp.ind){
    for (i in seq_along(n.seq)){
      n.guess <- n.seq[i]
      m <- ceiling(N/n.guess)
      group.ind <- rep(seq_len(m), each=n.guess) #group 1 to m
      group.ind1 <- group.ind[seq_len(N)] #truncation
      dat <- data.frame(y = y.seq, group = as.factor(group.ind1))
      val.seq[i] <- f.obj(dat)
    }
  } else {
    y <- y.seq
    for (i in seq_along(n.seq)){
      #change it to check overlapping groups 
      n <- n.seq[i]
      m <- floor((N-n)/step)
      
      ind.mat1 <- matrix(rep(1:n, m+1), ncol = n, byrow = TRUE)
      ind.mat2 <- matrix(rep((0:m)*step, n), ncol = n, byrow = FALSE)
      ind.mat <- ind.mat1 + ind.mat2
      
      y.mat <- matrix(NA, ncol = n, nrow = m+1)
      for (j in seq_len(m+1)){ #use different index j
        y.mat[j,] <- y[ind.mat[j,]]
      }
      
      group.ind <- rep(seq_len(m+1), each=n) 
      dat <- data.frame(y = as.vector(t(y.mat)), group = as.factor(group.ind))
      val.seq[i] <- f.obj(dat)
    }
  }

  
  if(plot.ind){
    plot(n.seq, val.seq, type = "l")
    #scatter.smooth(n.seq, val.seq, lpars = list(col = "red", lty = 3, lwd = 3))
  }
  
  
  
  if (save.ind){
    names1 <- paste0("choose-gr-size-", name.char, "-onepath.pdf")
    
    pdf(file = names1, width = width, height = height)
    plot(n.seq, val.seq, type = "l")
    dev.off()
    
  }
  
  n.opt.ind <- which.min(val.seq)
  n.opt <- n.seq[n.opt.ind]
  
  if (print.re.ind){
    print(paste("The suggested group size by this algorithm is", n.opt))
  }
  
  return(list(val.seq=val.seq, n.opt=n.opt))
}


#also check the change point detection 

choose.gr.size <- function(f.noise.gen, N, 
                           n.seq, f.obj, rep.time = 5e2, 
                           f.trans = function(x) x,
                           plot.ind = TRUE, 
                           save.ind = FALSE, name.char = NULL, 
                           width = 6, height = 4, 
                           print.re.ind = TRUE, 
                           center.val.mat = TRUE, 
                           ylab = "centered BIC", 
                           novlp.ind = TRUE, 
                           step = 10){
  #f.noise.gen: the function to generate the noise sequence
  #N: the length of the noise
  #n.seq: the sequence of candidate group sizes
  #f.obj: the objective function we can compute under each group size n 
  #rep.time: the times of replications to draw n.seq vs val.seq
  #center.val.mat: whether to move the B(n) path start from the same level
  if(novlp.ind){
    val.mat <- replicate(rep.time, {
      ep.seq <- f.noise.gen(N)
      y.seq <- f.trans(ep.seq)
      # re <- ar1.gen.sigdepn.noblock(N, rho = 0)
      # y.seq <- log(abs(re$w.seq)) #box-cox transform
      val.seq <- numeric(length(n.seq))
      for (i in seq_along(n.seq)){
        n.guess <- n.seq[i]
        m <- ceiling(N/n.guess)
        group.ind <- rep(seq_len(m), each=n.guess) #group 1 to m
        group.ind1 <- group.ind[seq_len(N)] #truncation
        dat <- data.frame(y = y.seq, group = as.factor(group.ind1))
        val.seq[i] <- f.obj(dat)
        ##one-way univariate ANOVA (general)
        #re <- oneway.test(y~group, data = dat)
        #p.val <- re$p.value
        
        ##one-way univariate ANOVA
        # fit.aov <- aov(y~group, data = dat)
        # re.aov <- summary(fit.aov)[[1]]
        # pv.f.test <- re.aov$`Pr(>F)`
        # p.val <- pv.f.test[1]
        
        # p.seq[i] <- log(p.val)
      }
      val.seq
    })
  } else {
    val.mat <- replicate(rep.time, {
      ep.seq <- f.noise.gen(N)
      y.seq <- f.trans(ep.seq)
      # re <- ar1.gen.sigdepn.noblock(N, rho = 0)
      # y.seq <- log(abs(re$w.seq)) #box-cox transform
      val.seq <- numeric(length(n.seq))
      
      y <- y.seq
      for (i in seq_along(n.seq)){
        #change it to check overlapping groups 
        n <- n.seq[i]
        m <- floor((N-n)/step)
        
        ind.mat1 <- matrix(rep(1:n, m+1), ncol = n, byrow = TRUE)
        ind.mat2 <- matrix(rep((0:m)*step, n), ncol = n, byrow = FALSE)
        ind.mat <- ind.mat1 + ind.mat2
        
        y.mat <- matrix(NA, ncol = n, nrow = m+1)
        for (j in seq_len(m+1)){ #use different index j
          y.mat[j,] <- y[ind.mat[j,]]
        }
        
        group.ind <- rep(seq_len(m+1), each=n) 
        dat <- data.frame(y = as.vector(t(y.mat)), group = as.factor(group.ind))
        val.seq[i] <- f.obj(dat)
      }
      
      val.seq
    })
  }
  
  
  if(center.val.mat){
    val.mat <- apply(val.mat, 2, function(x) x-mean(x))
  }
  
  avg.seq <- apply(val.mat, 1, mean)
  #record n.opt for each path
  n.opt.seq <- apply(val.mat, 2, function(x) n.seq[which.min(x)])
  #compute the mean and sd of the n.opt.seq
  #n.opt.mean <- round(mean(n.opt.seq))
  n.opt.med <- round(median(n.opt.seq))
  #n.opt.sd <- round(sd(n.opt.seq))
  n.opt.mad <- round(mad(n.opt.seq))
  
  if (plot.ind){
    matplot(n.seq, val.mat, type = "l", ylab = ylab)
    #pdf(file = "choose-gr-size-BIC-2.pdf", width = 6, height = 4)
    matplot(n.seq, val.mat, type = "l",col = alpha("lightblue", 0.2), ylab = ylab)
    lines(n.seq, avg.seq, col=2, lwd = 1.5)
    #dev.off()
  }
  
  #save the plot if needed
  if (save.ind){
    names1 <- paste0("choose-gr-size-", name.char, "-1.pdf")
    names2 <- paste0("choose-gr-size-", name.char, "-2.pdf")
    
    pdf(file = names1, width = width, height = height)
    matplot(n.seq, val.mat, type = "l", ylab = ylab)
    dev.off()
    
    pdf(file = names2, width = width, height = height)
    matplot(n.seq, val.mat, type = "l",col = alpha("lightblue", 0.2), ylab = ylab)
    lines(n.seq, avg.seq, col=2, lwd = 1.5)
    dev.off()
  }
  
  if (print.re.ind){
    print(paste("The suggested group size (median) by this algorithm is", n.opt.med, "(MAD:",n.opt.mad,")."))
  #print(paste("The suggested group size by this algorithm is", n.opt.mean, "(sd:",n.opt.sd,")."))
  }
  
  return(list(val.mat = val.mat, avg.seq = avg.seq, n.opt.seq = n.opt.seq))
  
}

f.BIC <- function(dat, short.ind = FALSE){
  N <- nrow(dat)
  m <- length(unique(dat$group))
  if (short.ind){
    mod.fit <- lm(y~group, data = dat)
    #why the BIC curves seems to parallel with each other
    #the likelihood is fixed 
    re <- BIC(mod.fit)
  } else {
    dat.split <- split(dat$y, dat$group)
    var.gr <- sapply(dat.split, function(x) mean(x^2))
    #compute the likelihood 
    loglik.gr.vec <- sapply(seq_along(var.gr), function(k){
      v2 <- var.gr[k]
      v <- sqrt(var.gr[k])
      y.vec <- dat.split[[k]]
      #C <- (-1/2)*log(2*pi)
      vec <- - log(v) - (y.vec^2)/(2*v2)
      #sum(C+vec)
      sum(vec)
    })
    loglik.full <- sum(loglik.gr.vec)
    re <- -2*loglik.full + log(N)*m
  }
  return(re)
}

f.BIC.full <- function(dat, short.ind = FALSE){
  #also consider changing mean
  N <- nrow(dat)
  m <- length(unique(dat$group))
  if (short.ind){
    mod.fit <- lm(y~group, data = dat)
    #why the BIC curves seems to parallel with each other
    #the likelihood is fixed 
    re <- BIC(mod.fit)
  } else {
    dat.split <- split(dat$y, dat$group)
    mean.gr <- sapply(dat.split, mean)
    var.gr <- sapply(dat.split, function(x) mean((x-mean(x))^2))
    #compute the likelihood 
    loglik.gr.vec <- sapply(seq_along(var.gr), function(k){
      u <- mean.gr[k]
      v2 <- var.gr[k]
      v <- sqrt(var.gr[k])
      y.vec <- dat.split[[k]]
      #C <- (-1/2)*log(2*pi)
      vec <- - log(v) - ((y.vec-u)^2)/(2*v2)
      #sum(C+vec)
      sum(vec)
    })
    loglik.full <- sum(loglik.gr.vec)
    re <- -2*loglik.full + log(N)*m
  }
  return(re)
}


f.AIC <- function(dat, short.ind = FALSE){
  N <- nrow(dat)
  m <- length(unique(dat$group))
  if (short.ind){
    mod.fit <- lm(y~group, data = dat)
    #why the BIC curves seems to parallel with each other
    #the likelihood is fixed 
    re <- AIC(mod.fit)
  } else {
    dat.split <- split(dat$y, dat$group)
    var.gr <- sapply(dat.split, function(x) mean(x^2))
    #compute the likelihood 
    loglik.gr.vec <- sapply(seq_along(var.gr), function(k){
      v2 <- var.gr[k]
      v <- sqrt(var.gr[k])
      y.vec <- dat.split[[k]]
      #C <- (-1/2)*log(2*pi)
      vec <- - log(v) - (y.vec^2)/(2*v2)
      #sum(C+vec)
      sum(vec)
    })
    loglik.full <- sum(loglik.gr.vec)
    re <- -2*loglik.full + 2*m
  }
  return(re)
}

#f.ks.test

f.ks.unif <- function(dat, gr.f = function(x) sqrt(mean(x))){
  dat.split <- split(dat$y, dat$group)
  mean.gr <- sapply(dat.split, gr.f)
  num.gr <- length(mean.gr)
  re.test <- ks.test(mean.gr, punif)
  p.val <- re.test$p.value
  #val <- -log(p.val)
  D <- re.test$statistic
  #check the distance to uniform
  #to see how spread out the empirical distribution is
  #val <- D/sqrt(num.gr)
  val <- D
  val 
  #the goal is to obtain the largest p.value or the smallest val
  #also check overlapping groups
}

#check other test related to relative entropy 

#f.pv.F.test 

#f.adj.R

