# Simulate data
simulateData <- function(b,a,s, tau2, n, kernel, timemax = 25, fixrange = FALSE){
  increment =1
  if(isTRUE(fixrange)){
    t1_tilde <- matrix(seq(0, timemax, length.out = n), ncol=1)
  }else{
    t1_tilde <- matrix(seq(0, 50, length.out = n), ncol=1)
  }
  t1 <- t1_tilde + runif(n,-increment/4,increment/4)
  # concatenate two time series time together
  t <- rep(t1,2)
  group.index <- rep(c(1,2),each=n)
  group.size <- rle(group.index)
  group.size$values <- s
  s.expand <- inverse.rle(group.size)
  # add time shift to the second time series
  tinput <- t+s.expand
  # create group matrix
  delta = Matrix(1-bdiag(replicate(2,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
  A = a^2*delta+1
  if(kernel == "seperate"){
    # perform exponential kernel(matern1/2)
    tinput <- tinput+runif(2*n,-increment/4,increment/4)
    D <- sqrt(distance(tinput))
    Vk = tau2 * exp(-b * D)
    g = 0.1
    sigma = rbind(cbind(Vk[1:n,1:n] + diag(g, n),matrix(0,n,n)),cbind(matrix(0,n,n),Vk[(n+1):(2*n), (n+1):(2*n)]+diag(g, n)))
  }
  if(kernel == "exp"){
    # perform LExp
    D <- sqrt(plgp::distance(tinput))
    C <- exp(-b*D) ## covariance matrix
    Vk = 1/A * C
    sigma = tau2*Vk
  }
  if(kernel == "rbf"){
    # perform LRBF
    D <- plgp::distance(tinput) / A
    C <- exp(-b*D) ## covariance matrix
    Vk = (1/sqrt(A)*C)
    sigma = tau2*Vk
  }
  if(kernel == "matern"){
    # perform LMat 3/2
    v = 3/2
    dist = sqrt(plgp::distance(tinput))
    dist[dist == 0] <- 1e-8
    part1 = tau2*(b*dist)^v
    part2 = (A^(v+0.5))
    part3 = besselK(b*dist, nu = v)
    sigma = part1 * part3 /part2
  }
  ## generate 100 replicates under a prior with MVN(0,sigma)
  Y <- rmvnorm(100, sigma= sigma %>% as.matrix()) %>% t()
  rownames(Y) <- tinput - s.expand
  return(Y)
}

# Visualize two time series in i replicate
## s is the time lag of the first and second time series, respectatively
plotSimulate <- function(Y, i = 1, s = c(0,2)){
  n = nrow(Y)/2
  cumsum.index <- cumsum(c(n,n))
  t1 = rownames(Y)[1:n] %>% as.numeric()
  matplot(t1+s[1],Y[1:cumsum.index[1],i],pch = 16, col="gray", lty=1, xlab="x", ylab="y",ylim = c(-10,10))
  points(t1+s[2],Y[(cumsum.index[1]+1):(cumsum.index[2]),i],lty=2, col = "blue")
}


## LRBF loglikelihood without vertical shift and jitter for simulation
nl_rbf <- function(par, t, Y, delta, group.size) # group.size is rle of group.index
{
  ngroup <- length(group.size$values)
  b <- softplus(par[1])
  a <- softplus(par[2])
  s <- par[3:(1+ngroup)]
  k <- length(Y)
  group.size$values <- c(0,s)
  s.expand <- inverse.rle(group.size) # input t vector
  tinput <- t+s.expand
  A = (a^2)*delta+1
  D <- plgp::distance(tinput)/A
  C <- exp(-b*D) ## covariance matrix
  Vk = (1/sqrt(A)*C) %>% as.matrix()
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}

  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y / k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a, -ll, "\n")
  return(-ll)
}


## LExp loglikelihood without vertical shift and jitter for simulation
nl_exp <- function(par, t, Y, delta, group.size) # group.size is rle of group.index
{
  ngroup <- length(group.size$values)
  b <- softplus(par[1])
  a <- softplus(par[2])
  s <- par[3:(1+ngroup)]
  k <- length(Y)
  group.size$values <- c(0,s)
  s.expand <- inverse.rle(group.size) # input t vector
  tinput <- t+s.expand
  A = (a^2)*delta+1
  D <- sqrt(plgp::distance(tinput))
  C <- exp(-b*D) ## covariance matrix
  Vk = (1/A*C) %>% as.matrix()
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}

  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y / k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a, -ll, "\n")
  return(-ll)
}


## LMat3/2 loglikelihood without vertical shift and jitter for simulation
nl_matern <- function(par, t, Y, delta, group.size) # group.size is rle of group.index
{
  ngroup <- length(group.size$values)
  b <- softplus(par[1])
  a <- softplus(par[2])
  s <- par[3:(1+ngroup)]
  k <- length(Y)
  group.size$values <- c(0,s)
  s.expand <- inverse.rle(group.size) # input t vector
  tinput <- t+s.expand
  A = (a^2)*delta+1
  dist = sqrt(plgp::distance(tinput))
  dist[dist == 0] <- 1e-8
  v = 3/2
  part1 = tau2*(b*dist)^v
  part2 = (A^(v+0.5))
  part3 = besselK(b*dist, nu = v)
  Vk = part1 * part3 /part2
  Vi <- pinv(Vk %>% as.matrix())
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}

  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y / k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a, -ll, "\n")
  return(-ll)
}

## LRBF loglikelihood with vertical shift and jitter for real application
nl4 <- function(par, t, Y, delta,group.size) # group.size is rle of group.index
{
  ngroup <- length(group.size$values)
  b <- softplus(par[1])
  a <- softplus(par[2])
  s <- par[3:(1+ngroup)]
  mu <- par[(2+ngroup):(2*ngroup)]
  g <- par[2*ngroup+1]
  k <- length(Y)
  group.size$values <- c(0,s)
  s.expand <- inverse.rle(group.size) # input t vector
  tinput <- t+s.expand
  A = (a^2)*delta+1
  D <- plgp::distance(tinput)/A
  C <- exp(-b*D) ## covariance matrix
  Vk = (1/sqrt(A)*C) %>% as.matrix()+diag(g,k)
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}

  group.size$values <- c(0,mu)
  mu.expand <- inverse.rle(group.size)
  Y_shift <- Y - mu.expand
  ll <- - (k/2)*log(t(Y_shift) %*% Vi %*% Y_shift) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a, s, -ll, "\n")
  return(-ll)
}


## GP RBF loglikelihood with for simulation
nlgp <- function(par, t, Y) # group.size is rle of group.index
{
  # b <- exp(par) ##comment it this if not estimate b
  b <- softplus(par)
  k <- length(Y)
  D <- plgp::distance(t)
  C <- exp(-b*D) ## covariance matrix
  Vi <- pinv(C)
  ldetK <- determinant(C, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}
  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y / k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b, mu, -ll, "\n")
  return(-ll)
}

## GP exp kernel (matern1/2) loglikelihood with for simulation
nlgp.matern <- function(par, t, Y) # group.size is rle of group.index
{
  b <- softplus(par)
  v = 1/2
  D <-  sqrt(plgp::distance(t))
  D[D == 0] <- 1e-8
  part1 = (b*D)^v
  part3 = besselK(b*D, nu = v)
  Vk = (part1 * part3) %>% as.matrix()
  k <- length(Y)
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}
  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y /k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b, mu, -ll, "\n")
  return(-ll)
}

## derive parameter estimation
deriveEstimation <-  function(Y, kernel, sl = 0, su = 4){
  k = nrow(Y)
  res = matrix(0,ncol = 5,nrow = ncol(Y))
  n = nrow(Y)/2
  t = rownames(Y)[1:n] %>% as.numeric()
  group.index <- rep(c(1,2),each=n)
  group.size <- rle(group.index)
  delta = Matrix(1-bdiag(replicate(2,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
  if(kernel == "sep.rbf"){
    # perform RBF kernel parameter estimation for SGP
    res = matrix(0,ncol = 3,nrow = ncol(Y))
    t1 = rownames(Y)[1:n] %>% as.numeric()
    t2 = rownames(Y)[(n+1):(2*n)] %>% as.numeric()
    for (i in 1:ncol(Y)) {
      yi = Y[,i]
      y1 = yi[1:n]
      y2 = yi[(n+1):(2*n)]
      out1 <- optim(1, nlgp, method="L-BFGS-B",lower=-10,
                   upper=20,
                   t= t1, Y=y1)
      bhat1 = softplus(out1$par)
      out2 <- optim(1, nlgp, method="L-BFGS-B",lower=-10,
                   upper=20,
                   t= t2, Y=y2)
      bhat2 = softplus(out2$par)

      D1 <-  plgp::distance(t1)
      Chat1 <- exp(-bhat1*D1) ## covariance matrix
      Vihat1 <- pinv(Chat1) %>% as.matrix()
      ldetK1 <- determinant(Chat1, logarithm=TRUE)$modulus
      ll1 <- - (n/2)*log(t(y1) %*% Vihat1 %*% y1 /n) - (1/2)*ldetK1 -n/2

      D2 <-  plgp::distance(t2)
      Chat2 <- exp(-bhat2*D2) ## covariance matrix
      Vihat2 <- pinv(Chat2) %>% as.matrix()
      ldetK2 <- determinant(Chat2, logarithm=TRUE)$modulus
      ll2 <- - (n/2)*log(t(y2) %*% Vihat2 %*% y2 /n) - (1/2)*ldetK2 -n/2
      res[i,]<- c(bhat1,bhat2,ll1+ll2)
    }
  }
  if(kernel == "sep.exp"){
    # perform exponential kernel parameter estimation for SGP
    res = matrix(0,ncol = 3,nrow = ncol(Y))
    t1 = rownames(Y)[1:n] %>% as.numeric()
    t2 = rownames(Y)[(n+1):(2*n)] %>% as.numeric()
    for (i in 1:ncol(Y)) {
      yi = Y[,i]
      y1 = yi[1:n]
      y2 = yi[(n+1):(2*n)]
      out1 <- optim(1, nlgp.matern, method="L-BFGS-B",lower=-10,
                    upper=20,
                    t= t1, Y=y1)
      bhat1 = softplus(out1$par)
      out2 <- optim(1, nlgp.matern, method="L-BFGS-B",lower=-10,
                    upper=20,
                    t= t2, Y=y2)
      bhat2 = softplus(out2$par)

      # MLE estimates tau2hat
      v = 1/2
      D1 <-  sqrt(plgp::distance(t1))
      D1[D1 == 0] <- 1e-8
      part1 = (sqrt(2*v) * bhat1*D1)^v
      part2 = 2^(1-v) / gamma(v)
      part3 = besselK(sqrt(2*v) * bhat1*D1, nu = v)
      Vk1 = (part1 * part2 *part3) %>% as.matrix()

      Vihat1 <- pinv(Vk1) %>% as.matrix()
      ldetK1 <- determinant(Vk1, logarithm=TRUE)$modulus
      ll1 <- - (n/2)*log(t(y1) %*% Vihat1 %*% y1 /n) - (1/2)*ldetK1 -n/2

      D2 <-  sqrt(plgp::distance(t2))
      D2[D2 == 0] <- 1e-8
      part1 = (sqrt(2*v) * bhat2*D2)^v
      part3 = besselK(sqrt(2*v) * bhat2*D2, nu = v)
      Vk2 = (part1 * part2 *part3) %>% as.matrix()
      Vihat2 <- pinv(Vk2) %>% as.matrix()
      ldetK2 <- determinant(Vk2, logarithm=TRUE)$modulus
      ll2 <- - (n/2)*log(t(y2) %*% Vihat2 %*% y2 /n) - (1/2)*ldetK2 -n/2
      res[i,]<- c(bhat1,bhat2,ll1+ll2)
    }
  }
  if(kernel == "exp"){
    # perform LExp parameter estimation
    for (i in 1:ncol(Y)) {
        yi = Y[,i]
        out <- optim(c(1,1,1), nl_exp, method="L-BFGS-B",lower=c(-5,-5,sl),
                     upper=c(20,Inf,su),
                     t= t, Y=yi, delta=delta,group.size = group.size)
        bhat = softplus(out$par[1])
        ahat = softplus(out$par[2])
        shat = out$par[3]
        # MLE estimates tau2hat
        group.size$values <- c(0,shat)
        s.expand <- inverse.rle(group.size) # input t vector
        that <- t+s.expand
        Ahat = (ahat^2)*delta+1
        Dhat <-  sqrt(plgp::distance(that))
        Chat <- exp(-bhat*Dhat) ## covariance matrix
        Vkhat = (1/Ahat*Chat) %>% as.matrix()

        Vihat <- pinv(Vkhat)
        ldetK <- determinant(Vkhat, logarithm=TRUE)$modulus
        tau2hat2 <- drop(t(yi) %*% Vihat %*% yi / length(that))
        ll <- - (k/2)*log(t(yi) %*% Vihat %*% yi / k) - (1/2)*ldetK - k/2

        res[i,]<- c(bhat,ahat,shat,tau2hat2, ll)
      }
  }
  if(kernel == "rbf"){
    # perform LRBF parameter estimation
    t = rownames(Y) %>% as.numeric()
    for (i in 1:ncol(Y)) {
      yi = Y[,i]
      out <- optim(c(0,1,2), nl_rbf, method="L-BFGS-B",lower=c(-10,-10,sl),
                   upper=c(20,Inf,su),
                   t= t, Y=yi, delta=delta,group.size = group.size)
      bhat = softplus(out$par[1])
      ahat = softplus(out$par[2])
      shat = out$par[3]
      # MLE estimates tau2hat
      group.size$values <- c(0,shat)
      s.expand <- inverse.rle(group.size) # input t vector
      that <- t+s.expand
      Ahat = (ahat^2)*delta+1
      Dhat <- plgp::distance(that) / Ahat
      Chat <- exp(-bhat*Dhat) ## covariance matrix
      Vkhat = (1/sqrt(Ahat)*Chat) %>% as.matrix()
      Vihat <- pinv(Vkhat)
      ldetK <- determinant(Vkhat, logarithm=TRUE)$modulus
      tau2hat2 <- drop(t(yi) %*% Vihat %*% yi / length(that))
      ll <- - (k/2)*log(t(yi) %*% Vihat %*% yi / k) - (1/2)*ldetK -k/2
      res[i,]<- c(bhat,ahat,shat,tau2hat2, ll)
    }
  }
  if(kernel == "matern"){
    # perform LMat parameter estimation
    for (i in 1:ncol(Y)) {
      yi = Y[,i]
      out <- optim(c(1,1,1), nl_matern, method="L-BFGS-B",lower=c(-10,-10,sl),
                   upper=c(20,Inf,su),
                   t= t, Y=yi, delta=delta,group.size = group.size)
      bhat = softplus(out$par[1])
      ahat = softplus(out$par[2])
      shat = out$par[3]
      # MLE estimates tau2hat
      group.size$values <- c(0,shat)
      s.expand <- inverse.rle(group.size) # input t vector
      that <- t+s.expand
      Ahat = (ahat^2)*delta+1
      v = 3/2
      Dhat <-  sqrt(plgp::distance(that))
      Dhat[Dhat == 0] <- 1e-8
      part1 = (bhat*Dhat)^v
      part2 = (Ahat^(v+0.5))
      part3 = besselK(bhat*Dhat, nu = v)
      Vkhat = (part1 * part3 /part2) %>% as.matrix()

      Vihat <- pinv(Vkhat)
      ldetK <- determinant(Vkhat, logarithm=TRUE)$modulus
      tau2hat2 <- drop(t(yi) %*% Vihat %*% yi / length(that))
      ll <- - (k/2)*log(t(yi) %*% Vihat %*% yi / k) - (1/2)*ldetK

      res[i,]<- c(bhat,ahat,shat,tau2hat2, ll)
    }
  }
  return(res)
}

# helper function for prediction applying LRBF kernel
## ntest : the number of test points
## delta: group matrix
## group.size : RLE vector indicate time series labels
## params : output from deriveEstimation
interpolation.rbf <- function(t, Y_shift, delta, group.size, params, ntest, plot = TRUE){
  tt <- matrix(seq(min(t), max(t), length.out=ntest), ncol=1)
  # tt2 <- tt+shat # use shat=0 if directly model the predicted value after shift
  ttest <- rbind(tt,tt) # concatenate two time series time together
  group.index.test <- rep(c(1,2),each=ntest)
  group.size.test <- rle(group.index.test)
  # test sets group matrix
  deltatest = Matrix(1-bdiag(replicate(2,matrix(1,ntest,ntest),simplify=FALSE)),sparse = TRUE)

  bhat <- params[1]
  ahat <- params[2]
  shat <- params[3]
  muhat <- params[4]
  # ghat <- params[5]
  tau2hat2 <- params[5]

  ## derive original inverse covariance kernel
  t2hat <- t+shat
  Ahat = (ahat^2)*delta+1
  D <-  plgp::distance(c(t,t2hat))/Ahat
  Chat <- exp(-bhat*D)
  Vk = (1/sqrt(Ahat)*Chat) %>% as.matrix() ## covariance matrix
  Vi <- pinv(Vk) ## psduo-inverse covariance matrix
  Vi <- (1/2)*(Vi+t(Vi))
  group.size$values <- c(0,muhat)
  mu.expand <- inverse.rle(group.size)
  Y_shifthat <- Y_shift - mu.expand

  ## derive test sets inverse covariance kernel
  Ahat = (ahat^2)*deltatest+1
  Dtest <- plgp::distance(ttest)/Ahat
  Ctt <- exp(-bhat*Dtest) ## covariance matrix
  Vtt = (1/sqrt(Ahat)*Ctt) %>% as.matrix() + diag(ghat,2*ntest)

  ## derive off diagonal inverse covariance kernel
  deltaoff = Matrix(1-bdiag(replicate(2,matrix(1,ntest,n),simplify=FALSE)),sparse = TRUE)
  Ahatoff = (ahat^2)*deltaoff+1
  Dtestoff <- plgp::distance(ttest,c(t,t2hat))/Ahatoff
  Ct <- exp(-bhat*Dtestoff) ## covariance matrix
  Vt = (1/sqrt(Ahatoff)*Ct) %>% as.matrix()

  mup2_shift <- Vt %*% Vi %*% (Y_shifthat)
  Sigmap2 <- tau2hat2*(Vtt-Vt%*%Vi%*%t(Vt))
  Sigmap2 <- (1/2)*(Sigmap2+t(Sigmap2))

  ## 100 prediction samples based on derived posterior mean and variance
  YY <- rmvnorm(100, mup2_shift, Sigmap2)
  ## correlation based on 100 prediction samples
  cor_mean <- cor(mup2_shift[1:ntest],mup2_shift[(ntest+1):(2*ntest)])
  ## derive correlation Credible interval
  cor <- sapply(1:100, function(j){
    return(cor(YY[j,1:ntest],YY[j,(ntest+1):(2*ntest)]))
  })

  if(isTRUE(plot)){
    YY <- YY %>% t() %>% as.data.frame()
    YY$group <- c(rep("gene",ntest),rep("peak",ntest))
    YY$t <- ttest
    YY2 <- YY %>% pivot_longer(cols=1:100, names_to = "iter")

    p2 <- ggplot(YY2,aes(t,value, color=iter)) +
      geom_line()+facet_wrap(~group)+
      theme_classic(base_size = 16)+
      # scale_color_igv()+
      theme(legend.position="none",
            panel.grid.minor = element_blank(),
            panel.border = element_rect(fill = 'transparent'))
    print(p2)
  }

  return(list(cor_mean,cor))
}


# helper function for simulation prediction which only need posterior mean
interpolation_sim <- function(t,Y,params,ttest, kernel = "lexp", prediction = TRUE){
  ntest = length(ttest)/2
  # tt2 <- tt+shat # use shat=0 if directly model the predicted value after shift

  ## derive original inverse covariance kernel
  n = length(Y)/2
  t1 <- t[1:n]
  t2 <- t[(n+1):(2*n)]

  if(kernel == "lexp"){
    bhat <- params[1]
    ahat <- params[2]
    shat <- params[3]
    tau2hat2 <- params[4]

    t2hat <- t2+shat
    delta = Matrix(1-bdiag(replicate(2,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
    deltatest = Matrix(1-bdiag(replicate(2,matrix(1,ntest,ntest),simplify=FALSE)),sparse = TRUE)
    group.index.test <- rep(c(0,shat),each=ntest)
    if(isTRUE(prediction)){
      ttesthat <- ttest + group.index.test
    }else{
      ttesthat <- ttest
    }

    Ahat = (ahat^2)*delta+1
    D <-  sqrt(plgp::distance(c(t1,t2hat)))
    Chat <- exp(-bhat*D)
    Vk = (1/Ahat*Chat) %>% as.matrix() ## covariance matrix
    Vi <- pinv(Vk) ## psduo-inverse covariance matrix

    deltaoff = Matrix(1-bdiag(replicate(2,matrix(1,ntest,n),simplify=FALSE)),sparse = TRUE)
    Ahatoff = (ahat^2)*deltaoff+1
    Dtestoff <- sqrt(plgp::distance(ttesthat,c(t1,t2hat)))
    Ct <- exp(-bhat*Dtestoff) ## covariance matrix
    Vt = (1/Ahatoff*Ct) %>% as.matrix()
    mup2_shift <- Vt %*% Vi %*% Y

  }
  if(kernel == "lrbf"){
    bhat <- params[1]
    ahat <- params[2]
    shat <- params[3]
    tau2hat2 <- params[4]

    t2hat <- t2+shat
    delta = Matrix(1-bdiag(replicate(2,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
    deltatest = Matrix(1-bdiag(replicate(2,matrix(1,ntest,ntest),simplify=FALSE)),sparse = TRUE)
    group.index.test <- rep(c(0,shat),each=ntest)
    if(isTRUE(prediction)){
      ttesthat <- ttest + group.index.test
    }else{
      ttesthat <- ttest
    }

    Ahat = (ahat^2)*delta+1
    D <-  plgp::distance(c(t1,t2hat))
    Chat <- exp(-bhat*D/Ahat)
    Vk = (1/sqrt(Ahat)*Chat) %>% as.matrix() ## covariance matrix
    Vi <- pinv(Vk) ## psduo-inverse covariance matrix

    deltaoff = Matrix(1-bdiag(replicate(2,matrix(1,ntest,n),simplify=FALSE)),sparse = TRUE)
    Ahatoff = (ahat^2)*deltaoff+1
    Dtestoff <- plgp::distance(ttesthat,c(t1,t2hat))
    Ct <- exp(-bhat*Dtestoff/Ahatoff) ## covariance matrix
    Vt = (1/sqrt(Ahatoff)*Ct) %>% as.matrix()
    mup2_shift <- Vt %*% Vi %*% Y

  }
  if(kernel == "sep.matern"){
    bhat1 <- params[1]
    bhat2 <- params[2]

    v = 1/2
    D1 <-  sqrt(plgp::distance(t1))
    D1[D1 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat1*D1)^v
    part2 = 2^(1-v) / gamma(v)
    part3 = besselK(sqrt(2*v) * bhat1*D1, nu = v)
    Vk1 = (part1 * part2 *part3) %>% as.matrix()

    Vihat1 <- pinv(Vk1) %>% as.matrix()

    Dtestoff1 <- sqrt(plgp::distance(ttest[1:n],t1))
    Dtestoff1[Dtestoff1 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat1*Dtestoff1)^v
    part3 = besselK(sqrt(2*v) * bhat1*Dtestoff1, nu = v)
    Vkoff1 = (part1 * part2 *part3) %>% as.matrix()
    mup2_shift1 <- Vkoff1 %*% Vihat1 %*% Y[1:n]

    ## Another group
    D2 <-  sqrt(plgp::distance(t2))
    D2[D2 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat2*D2)^v
    part3 = besselK(sqrt(2*v) * bhat2*D2, nu = v)
    Vk2 = (part1 * part2 *part3) %>% as.matrix()
    Vihat2 <- pinv(Vk2) %>% as.matrix()

    Dtestoff2 <- sqrt(plgp::distance(ttest[(n+1):(2*n)],t2))
    Dtestoff2[Dtestoff2 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat2*Dtestoff2)^v
    part3 = besselK(sqrt(2*v) * bhat2*Dtestoff2, nu = v)
    Vkoff2 = (part1 * part2 *part3) %>% as.matrix()
    mup2_shift2 <- Vkoff2 %*% Vihat2 %*% Y[(n+1):(2*n)]

    mup2_shift <- c(mup2_shift1,mup2_shift2)
  }

  return(mup2_shift)
}

