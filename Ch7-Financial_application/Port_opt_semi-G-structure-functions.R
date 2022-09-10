library(matrixcalc)
library(CVXR)

library(tidyverse)
library(plotly)
library(rgl)
library(plot3D)

A.f <- function(sd.mat, rho.intvl, par.weight, i=1){
  #sd.mat=rbind(sd.intvl1, sd.intvl2)
  sd.intvl1 <- sd.mat[1,]
  sd.intvl2 <- sd.mat[2,]
  vl1 <- sd.intvl1[1]
  vr1 <- sd.intvl1[2]
  vl2 <- sd.intvl2[1]
  vr2 <- sd.intvl2[2]
  
  rhl <- rho.intvl[1]
  rhr <- rho.intvl[2]
  
  a <- par.weight
  C.l <- matrix(c(vl1^2, vl1*vl2*rhl, 
                  vl1*vl2*rhl, vl2^2), 2, 2)
  C.r <- matrix(c(vr1^2, vr1*vr2*rhr, 
                  vr1*vr2*rhr, vr2^2), 2, 2)
  if (i==1){
    return(a*C.l+(1-a)*C.r)
  } 
  if (i==2){
    return((1-a/2)*C.l+(3/2*a-1)*C.r)
  }
}

rho.sig.prod <- function(rho,sig){
  re <- outer(rho,sig)
  as.numeric(re)
}

#min and max
range.set <- function(x1,x2=1-x1, 
                      rho.intvl = c(-0.8, -0.4),
                      sig1.intvl = c(0.5, 1),
                      sig2.intvl = c(1, 2), 
                      print.ind = FALSE){
  if (x1==0) return(x2*sig2.intvl^2)
  if (x1==1) return(x1*sig1.intvl^2)
  r.l <- rho.intvl[1]
  r.r <- rho.intvl[2]
  s1.l <- sig1.intvl[1]
  s1.r <- sig1.intvl[2]
  s2.l <- sig2.intvl[1]
  s2.r <- sig2.intvl[2]
  
  s1.set <- c(s1.l, s1.r)
  s2.set <- c(s2.l, s2.r)
  
  # s1.add <- (-x2/x1)*rho.sig.prod(rho.intvl, sig2.intvl)
  # s1.set <- c(s1.set, s1.add[s1.l < s1.add & s1.add < s1.r])
  # 
  # s2.add<- (-x1/x2)*rho.sig.prod(rho.intvl, sig1.intvl)
  # s2.set <- c(s2.set, s2.add[s2.l < s2.add & s2.add < s2.r])
  
  H.val <- NULL
  
  for (rho in c(r.l, r.r)){
    #for each rho, create the set to be added
    s1.add <- -(x2/x1)*rho*sig2.intvl
    s2.add <- -(x1/x2)*rho*sig1.intvl
    
    #only keep the values that are in the range of sigma
    s1.set <- c(s1.set, s1.add[s1.l < s1.add & s1.add < s1.r])
    s2.set <- c(s2.set, s2.add[s2.l < s2.add & s2.add < s2.r])
    
    if (print.ind){
      for (s1 in s1.set){
        for (s2 in s2.set){
          H.val <- c(H.val, H(sig1 = s1, sig2 = s2, 
                              rho = rho, 
                              x1 = x1, x2 = x2))
          i.end <- length(H.val)
          cat("i = ", i.end, ",",
              "sig1 = ", s1, ",", 
              "sig2 = ", s2, ",",
              "rho = ", rho, ",",
              "H = ", H.val[i.end], "\n")
        }
      }
    } else {
      for (s1 in s1.set){
        for (s2 in s2.set){
          H.val <- c(H.val, H(sig1 = s1, sig2 = s2, 
                              rho = rho, 
                              x1 = x1, x2 = x2))
        }
      }
    }
  }
  
  return(range(H.val))
  
}

cenamb2intvl <- function(cen, amb){
  r <- amb/2
  c(cen-r, cen+r)
}

lowamb2intvl <- function(low, amb){
  c(low, low+amb)
}
#draw EF


x.thres.f <- function(r.set){
  r.thres <- r.set
  x.thres <- 1/(1+r.thres)
  #sort(x.thres)
  x.thres
}

H <- function(sig1, sig2, rho = -0.5, x1=0.8, x2=1-x1){
  x1^2*sig1^2 + 2*x1*x2*sig1*sig2*rho + x2^2*sig2^2
}

H.min <- function(x1, x2=1-x1, sig1, sig2, rhol){
  #double check this function
  r <- x2/x1 
  #t.break 
  if (x1 <= t.break[1]){
    re.mat.c[i,] <- c(sdr1, sdl2)
  } else if (x1 <= t.break[2]) {
    re.mat.c[i,] <- c(r*sdl2/q, sdl2)
  } else if (x1 <= t.break[3]) {
    re.mat.c[i,] <- c(sdl1, sdl2)
  } else if (x1 <= t.break[4]) {
    re.mat.c[i,] <- c(sdl1, sdl1/(r*q))
  } else {
    re.mat.c[i,] <- c(sdl1, sdr2)
  }
}

#also provide H.max function

plot.ambcen.var <- function(sd1, sd2, rho, x.step = 1e-4, 
                            return.ind = FALSE, num.ind = TRUE, 
                            weight.obj = 0.5, 
                            plot.ind = TRUE,
                            plot.minmax.var = FALSE, 
                            weight.obj2 = 0.5){
  #weight.obj: the weight on amb, in objective function,
  #(then the weight on cen is 1- weight.obj)
  
  #weight.obj2: the weight on min, if plot.minmax.var is TRUE
  k <- weight.obj
  d <- x.step
  x.seq <- seq(d, 1-d,d)
  
  sd1.c <- mean(sd1)
  sd2.c <- mean(sd2)

  rhol  <- rho[1]
  rhor <-  rho[2]
  
  val.min.seq <- val.max.seq <- numeric(length(x.seq))
  
  if(num.ind){
    for (i in seq_along(x.seq)){
      x1 <- x.seq[i]
      #min var
      re.min <- optim(c(sd1.c,sd2.c), function(s) H(s[1],s[2], rho = rhol, x1 = x1), 
                      lower = c(sd1[1], sd2[1]), 
                      upper = c(sd1[2], sd2[2]), 
                      method = "L-BFGS-B")
      val.min.seq[i] <- re.min$value
      #max var
      re.max <- optim(c(sd1.c,sd2.c), function(s) -H(s[1],s[2], rho = rhor, x1 = x1), 
                      lower = c(sd1[1], sd2[1]), 
                      upper = c(sd1[2], sd2[2]), 
                      method = "L-BFGS-B")
      val.max.seq[i] <- -re.max$value
    }
    
    #matplot(x.seq, cbind(val.min.seq,val.max.seq), type = "l", lty = c(1,1), main = title.main)
    #center
    var.cen.seq <- (val.min.seq+val.max.seq)/2
    
    #plot(x.seq, var.cen.seq, type = "l",
    #     main = title.main)
    
    
    #amb
    var.amb.seq <- val.max.seq - val.min.seq
    
    # plot(x.seq, var.amb.seq, type = "l",
    #      main = title.main)
    obj.seq <- k*var.amb.seq + (1-k)*var.cen.seq
    x1.opt <- x.seq[which.min(obj.seq)]
    
    if(plot.ind){
      title.main <- paste0("sd1 = ", "[", sd1[1], ",", sd1[2],"], ", "sd2 = ", "[", sd2[1], ",", sd2[2],"]", ", ", 
                           "rho = ", "[", rhol, ",", rhor,"], ", "\n",
                           "obj = ", 
                           k, "*amb+(", 1-k, ")*cen, ", 
                           "x1.opt = ", x1.opt)
      x1 <- x.seq
      dat1 <- data.frame(x=var.cen.seq, y=var.amb.seq)
      p1 <- ggplot(dat1, aes(x=x, y=y)) + geom_point(aes(colour = x1)) + 
        xlab("Var.center") + ylab("Var.amb") + ggtitle(title.main) +
        scale_colour_gradient(low = "#00FFFF", high = "#000099")
      
      print(p1)
      
      if(plot.minmax.var){
        l <- weight.obj2
        obj.seq <- l*val.min.seq + (1-l)*val.max.seq
        x1.opt <- x.seq[which.min(obj.seq)]
        title.main <- paste0("sd1 = ", "[", sd1[1], ",", sd1[2],"], ", "sd2 = ", "[", sd2[1], ",", sd2[2],"]", ", ", 
                             "rho = ", "[", rhol, ",", rhor,"], ", "\n",
                             "obj = ", 
                             l, "*min_var+(", 1-l, ")*max_var, ", 
                             "x1.opt = ", x1.opt)
        x1 <- x.seq
        dat1 <- data.frame(x=val.min.seq, y=val.max.seq)
        p2 <- ggplot(dat1, aes(x=x, y=y)) + geom_point(aes(colour = x1)) + 
          xlab("Var.min") + ylab("Var.max") + ggtitle(title.main) +
          scale_colour_gradient(low = "#00FFFF", high = "#000099")
        
        print(p2)
      }
    }
    
    if(return.ind){
      return(list(x1.seq=x.seq, 
                  var.cen.seq = var.cen.seq, 
                  var.amb.seq=var.amb.seq, 
                  x1.opt = x1.opt))
    }
  } else {
    #use theoretical values
    for (i in seq_along(x.seq)){
      x1 <- x.seq[i]
      #min var
      re.min <- optim(c(sd1.c,sd2.c), function(s) H(s[1],s[2], rho = rhol, x1 = x1), 
                      lower = c(sd1[1], sd2[1]), 
                      upper = c(sd1[2], sd2[2]), 
                      method = "L-BFGS-B")
      val.min.seq[i] <- re.min$value
      
      #max var
      re.max <- optim(c(sd1.c,sd2.c), function(s) -H(s[1],s[2], rho = rhor, x1 = x1), 
                      lower = c(sd1[1], sd2[1]), 
                      upper = c(sd1[2], sd2[2]), 
                      method = "L-BFGS-B")
      val.max.seq[i] <- -re.max$value
    }
    
    #matplot(x.seq, cbind(val.min.seq,val.max.seq), type = "l", lty = c(1,1), main = title.main)
    #center
    var.cen.seq <- (val.min.seq+val.max.seq)/2
    
    #plot(x.seq, var.cen.seq, type = "l",
    #     main = title.main)
    
    
    #amb
    var.amb.seq <- val.max.seq - val.min.seq
    
    # plot(x.seq, var.amb.seq, type = "l",
    #      main = title.main)
    
    x1 <- x.seq
    dat1 <- data.frame(x=var.cen.seq, y=var.amb.seq)
    p1 <- ggplot(dat1, aes(x=x, y=y)) + geom_point(aes(colour = x1)) + 
      xlab("Var.center") + ylab("Var.amb") + ggtitle(title.main) +
      scale_colour_gradient(low = "#00FFFF", high = "#000099")
    
    print(p1)
    
    if(return.ind){
      return(list(x1.seq=x.seq, 
                  var.cen.seq = var.cen.seq, 
                  var.amb.seq=var.amb.seq))
    }
  }
  
}

#create the objective function 
#optimization procedure 
opt.ambcen <- function(sd1, sd2, rho, 
                       weight.obj = 0.5, x.range = c(1e-4, 1-1e-4)){
  #sd1, sd2, rho are intervals 
  #weight.obj: the weight on amb, in objective function,
  #(then the weight on cen is 1- weight.obj)
  #also use the five different regions to compute the min and max variance
  k <- weight.obj
  sd1.c <- mean(sd1)
  sd2.c <- mean(sd2)
  rhol  <- rho[1]
  rhor <-  rho[2]
  obj.f <- function(x){
    x1 <- x
    re.min <- optim(c(sd1.c,sd2.c), function(s) H(s[1],s[2], rho = rhol, x1 = x1), 
                    lower = c(sd1[1], sd2[1]), 
                    upper = c(sd1[2], sd2[2]), 
                    method = "L-BFGS-B")
    var.min <- re.min$value
    re.max <- optim(c(sd1.c,sd2.c), function(s) -H(s[1],s[2], rho = rhor, x1 = x1), 
                    lower = c(sd1[1], sd2[1]), 
                    upper = c(sd1[2], sd2[2]), 
                    method = "L-BFGS-B")
    var.max <- -re.max$value
    var.cen <- (var.min+var.max)/2
    var.amb <- var.max-var.min
    obj.val <- k*var.amb+(1-k)*var.cen
    return(obj.val)
  }
  re <- optimize(obj.f, interval = x.range)
  x1.opt <- re$minimum
  x2.opt <- 1-x1.opt
  return(list(x.opt = c(x1.opt, x2.opt), val.opt = re$objective, 
              weight.obj = weight.obj))
}

