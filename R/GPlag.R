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
  if(kernel == "exp.wa"){
    # perform LExp
    D <- sqrt(plgp::distance(tinput))
    C <- exp(-b*D/sqrt(A)) ## covariance matrix
    Vk = 1/sqrt(A) * C
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
plotSimulate <- function(Y, i = 1, s = c(0,2), pointSize = 1.5, textSize = 1.2){
  n = nrow(Y)/2
  cumsum.index <- cumsum(c(n,n))
  t1 = rownames(Y)[1:n] %>% as.numeric()
  matplot(t1+s[1], Y[1:cumsum.index[1],i], pch = 16, cex = pointSize, col="gray", lty=1, xlab="x", ylab="y", ylim = c(-10,10), cex.axis=textSize, cex.lab=textSize)
  points(t1+s[2], Y[(cumsum.index[1]+1):(cumsum.index[2]),i], lty=2, col = "black", cex = pointSize)
}

plotMultiSimulate <- function(Y, i = 1, s = c(0,2,4), pointSize = 1.5, textSize = 1.2){
  n = nrow(Y)/3
  cumsum.index <- cumsum(c(n,n,n))
  t1 = rownames(Y)[1:n] %>% as.numeric()
  matplot(t1+s[1], Y[1:cumsum.index[1],i], pch = 16, cex = pointSize, col="gray", lty=1, xlab="x", ylab="y", ylim = c(-10,10), cex.axis=textSize, cex.lab=textSize)
  points(t1+s[2], Y[(cumsum.index[1]+1):(cumsum.index[2]),i], lty=2, col = "black", cex = pointSize)
  points(t1+s[3], Y[(cumsum.index[2]+1):(cumsum.index[3]),i], lty=2, col = "blue", cex = pointSize)
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
  epsilon <- 1e-6
  Vk <- Vk + diag(epsilon, nrow(Vk))
  if (any(is.infinite(Vk)) || any(is.na(Vk))) {
    stop("Covariance matrix Vk contains infinite or NaN values.")
  }
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}

  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y / k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a, s, -ll, "\n")
  return(-ll)
}

softplusinv <- function(y) {
  log(exp(y) - 1)
}

softplus <- function(x, cap = 1e6) {
  result <- log(1 + exp(x))
  result[result > cap] <- cap
  return(result)
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
  # Debugging and validation
  # cat("Parameters: b =", b, ", a =", a, ", s =", s, "\n")
  
  A <- (a^2) * delta + 1
  if (any(is.infinite(A)) || any(is.na(A))) {
    stop("A contains infinite or NaN values.")
  }
  
  D <- sqrt(plgp::distance(tinput))
  if (any(is.infinite(D)) || any(is.na(D))) {
    stop("Distance matrix D contains infinite or NaN values.")
  }
  
  C <- exp(-b * D)
  if (any(is.infinite(C)) || any(is.na(C))) {
    stop("Covariance matrix C contains infinite or NaN values.")
  }
  
  Vk <- (1 / A * C) %>% as.matrix()
  if (any(is.infinite(Vk)) || any(is.na(Vk))) {
    stop("Covariance matrix Vk contains infinite or NaN values before regularization.")
  }
  
  epsilon <- 1e-6
  Vk <- Vk + diag(epsilon, nrow(Vk))
  if (any(is.infinite(Vk)) || any(is.na(Vk))) {
    stop("Covariance matrix Vk contains infinite or NaN values after regularization.")
  }
  
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm = TRUE)$modulus
  if (is.infinite(ldetK)) {
    ldetK <- -1e6
  }
  
  ll <- - (k / 2) * log(t(Y) %*% Vi %*% Y / k) - (1 / 2) * ldetK
  return(-ll)
}

## LExp with a in denominator loglikelihood
nl_exp.wa <- function(par, t, Y, delta, group.size) # group.size is rle of group.index
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
  C <- exp(-b*D/sqrt(A)) ## covariance matrix

  Vk = (1/sqrt(A)*C) %>% as.matrix()
  epsilon <- 1e-6
  Vk <- Vk + diag(epsilon, nrow(Vk))
  if (any(is.infinite(Vk)) || any(is.na(Vk))) {
    stop("Covariance matrix Vk contains infinite or NaN values.")
  }
  # cat(sum(is.na(Vk)), "\n", sum(is.infinite(A)), "\n")
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}

  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y / k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a,s, ldetK,-ll, "\n")
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
  epsilon <- 1e-6
  Vk <- Vk + diag(epsilon, nrow(Vk))
  if (any(is.infinite(Vk)) || any(is.na(Vk))) {
    stop("Covariance matrix Vk contains infinite or NaN values.")
  }
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

## LExp loglikelihood with vertical shift and jitter for real application
nl_exp_advanced <- function(par, t, Y, delta, group.size) # group.size is rle of group.index
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
  D <- sqrt(plgp::distance(tinput))
  C <- exp(-b*D) ## covariance matrix
  Vk = (1/A*C) %>% as.matrix() + diag(g, k)
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}

  group.size$values <- c(0,mu)
  mu.expand <- inverse.rle(group.size)
  Y_shift <- Y - mu.expand
  ll <- - (k/2)*log(t(Y_shift) %*% Vi %*% Y_shift) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a, s, mu, g, -ll, "\n")
  return(-ll)
}

nl_matern_advanced <- function(par, t, Y, delta, group.size) {
  ngroup <- length(group.size$values)
  b <- softplus(par[1])
  a <- softplus(par[2])
  s <- par[3:(1+ngroup)]
  mu <- par[(2+ngroup):(2*ngroup)]
  g <- par[2*ngroup+1]
  k <- length(Y)
  
  group.size$values <- c(0, s)
  s.expand <- inverse.rle(group.size)  # input t vector
  tinput <- t + s.expand
  
  A = (a^2) * delta + 1
  dist = sqrt(plgp::distance(tinput))
  dist[dist == 0] <- 1e-8
  
  v = 3/2
  part1 = (b * dist)^v
  part2 = (A^(v + 0.5))
  part3 = besselK(b * dist, nu = v)
  Vk = (part1 * part3 / part2) %>% as.matrix() + diag(g, k)  # add g to the diagonal for stability
  Vi <- pinv(Vk)
  
  ldetK <- determinant(Vk, logarithm = TRUE)$modulus
  if (is.infinite(ldetK)) {
    ldetK <- -1e6
  }
  
  group.size$values <- c(0, mu)
  mu.expand <- inverse.rle(group.size)
  Y_shift <- Y - mu.expand
  
  ll <- - (k / 2) * log(t(Y_shift) %*% Vi %*% Y_shift / k) - (1 / 2) * ldetK
  # counter <<- counter + 1
  # cat(b,a, s, mu, g, -ll, "\n")
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
deriveEstimation <-  function(Y, kernel, b0 = 1, a0=1, s0=1, bu =20, sl = 0, su = 4){
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
        out <- optim(c(b0,a0,s0), nl_exp, method="L-BFGS-B",lower=c(-5,-2,sl),
                     upper=c(bu,Inf,su),
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
  if(kernel == "exp.wa"){
    # perform LExp parameter estimation
    for (i in 1:ncol(Y)) {
      yi = Y[,i]
      out <- optim(c(b0,1,1), nl_exp.wa, method="L-BFGS-B",lower=c(-5,-2,sl),
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
      Chat <- exp(-bhat*Dhat/sqrt(Ahat)) ## covariance matrix
      Vkhat = (1/sqrt(Ahat)*Chat) %>% as.matrix()

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
      out <- optim(c(b0,1,2), nl_rbf, method="L-BFGS-B",lower=c(-5,-2,sl),
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
      out <- optim(c(b0,1,1), nl_matern, method="L-BFGS-B",lower=c(-10,-10,sl),
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

# helper function for prediction
## ntest : the number of test points
## delta: group matrix
## group.size : RLE vector indicate time series labels
## params : output from deriveEstimation
interpolation.kernel <- function(t, Y_shift, delta, group.size, params, ntest, kernel = "lrbf", plot = TRUE){
  tt <- matrix(seq(min(t), max(t), length.out=ntest), ncol=1)
  ttest <- rbind(tt, tt) # concatenate two time series time together
  group.index.test <- rep(c(1, 2), each = ntest)
  group.size.test <- rle(group.index.test)
  deltatest = Matrix(1 - bdiag(replicate(2, matrix(1, ntest, ntest), simplify = FALSE)), sparse = TRUE)

  bhat <- params[1]
  ahat <- params[2]
  shat <- params[3]
  muhat <- params[4]
  ghat <- params[5]
  tau2hat2 <- params[6]

  t2hat <- t + shat
  ## derive original inverse covariance kernel
  Ahat = (ahat^2) * delta + 1
  if (kernel == "lrbf") {
    D <- plgp::distance(c(t, t2hat)) / Ahat
    Chat <- exp(-bhat * D)
    Vk = (1 / sqrt(Ahat) * Chat) %>% as.matrix() ## covariance matrix
  } else if (kernel == "lexp") {
    D <- sqrt(plgp::distance(c(t, t2hat)))
    Chat <- exp(-bhat * D)
    Vk = (1/Ahat*Chat) %>% as.matrix()
  } else if (kernel == "lmatern") {
    v = 1.5  # Matérn ν parameter
    D <- sqrt(plgp::distance(c(t, t2hat)))
    D[D == 0] <- 1e-8  # Avoid division by zero
    part1 = (bhat * D)^v
    part2 = (Ahat^(v + 0.5))
    part3 = besselK(bhat * D, nu = v)
    Vk = (part1 * part3 / part2)%>% as.matrix()
  } else {
      stop("Unsupported kernel type")
  }

  Vi <- pinv(Vk) ## psuedo-inverse covariance matrix
  Vi <- (1 / 2) * (Vi + t(Vi))
  group.size$values <- c(0, muhat)
  mu.expand <- inverse.rle(group.size)
  Y_shifthat <- Y_shift - mu.expand

  ## derive test sets inverse covariance kernel
  Ahat = (ahat^2)*deltatest+1
  if (kernel == "lrbf") {
    Dtest <- plgp::distance(ttest)/Ahat
    Ctt <- exp(-bhat*Dtest) ## covariance matrix
    Vtt = (1/sqrt(Ahat)*Ctt) %>% as.matrix() + diag(ghat,2*ntest)
  } else if (kernel == "lexp") {
    Dtest <- sqrt(plgp::distance(ttest))
    Ctt <- exp(-bhat * Dtest)
    Vtt = (1/Ahat*Ctt) %>% as.matrix() + diag(ghat,2*ntest)
  } else if (kernel == "lmatern") {
    v = 1.5  # Matérn ν parameter
    Dtest <- sqrt(plgp::distance(ttest))
    Dtest[Dtest == 0] <- 1e-8  # Avoid division by zero
    part1 = (bhat * Dtest)^v
    part2 = (Ahat^(v + 0.5))
    part3 = besselK(bhat * Dtest, nu = v)
    Vtt = (part1 * part3 / part2)%>% as.matrix() + diag(ghat,2*ntest)
  } else {
      stop("Unsupported kernel type")
  }

  ## derive off diagonal inverse covariance kernel
  deltaoff = Matrix(1 - bdiag(replicate(2, matrix(1, ntest, length(t)), simplify = FALSE)), sparse = TRUE)
  Ahatoff = (ahat^2) * deltaoff + 1
  if (kernel == "lrbf") {
    Dtestoff <- plgp::distance(ttest,c(t,t2hat))/Ahatoff
    Ct <- exp(-bhat*Dtestoff) ## covariance matrix
    Vt = (1/sqrt(Ahatoff)*Ct) %>% as.matrix()
  } else if (kernel == "lexp") {
    Dtestoff <- sqrt(plgp::distance(ttest,c(t,t2hat)))
    Ct <- exp(-bhat * Dtestoff)
    Vt = (1/Ahatoff*Ct) %>% as.matrix()
  } else if (kernel == "lmatern") {
    v = 1.5  # Matérn ν parameter
    Dtestoff <- sqrt(plgp::distance(ttest,c(t,t2hat)))
    Dtestoff[Dtestoff == 0] <- 1e-8  # Avoid division by zero
    part1 = (bhat * Dtestoff)^v
    part2 = (Ahatoff^(v + 0.5))
    part3 = besselK(bhat * Dtestoff, nu = v)
    Vt = (part1 * part3 / part2)%>% as.matrix() 
  } else {
      stop("Unsupported kernel type")
  }

  mup2_shift <- Vt %*% Vi %*% (Y_shifthat)
  Sigmap2 <- tau2hat2 * (Vtt - Vt %*% Vi %*% t(Vt))
  Sigmap2 <- (1 / 2) * (Sigmap2 + t(Sigmap2))

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
  
  ## return results
  return(list(cor_mean,cor))
}

# helper function for simulation prediction which only need posterior mean
interpolation_sim <- function(t,Y,ytest,params,ttest, kernel = "lexp", prediction = TRUE){
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
    ## derive original inverse covariance kernel
    Ahat = (ahat^2)*delta+1
    D <-  sqrt(plgp::distance(c(t1,t2hat)))
    Chat <- exp(-bhat*D)
    Vk = (1/Ahat*Chat) %>% as.matrix() ## covariance matrix
    Vi <- pinv(Vk) ## psduo-inverse covariance matrix

    ## derive test sets inverse covariance kernel
    Ahat = (ahat^2)*deltatest+1
    Dtest <- sqrt(plgp::distance(ttest))
    Ctt <- exp(-bhat*Dtest) ## covariance matrix
    Vtt = (1/Ahat*Ctt) %>% as.matrix()

    ## derive off diagonal inverse covariance kernel
    deltaoff = Matrix(1-bdiag(replicate(2,matrix(1,ntest,n),simplify=FALSE)),sparse = TRUE)
    Ahatoff = (ahat^2)*deltaoff+1
    Dtestoff <- sqrt(plgp::distance(ttesthat,c(t1,t2hat)))
    Ct <- exp(-bhat*Dtestoff) ## covariance matrix
    Vt = (1/Ahatoff*Ct) %>% as.matrix()

    mup2_shift <- Vt %*% Vi %*% Y
    # prediction variance
    Sigmap2 <- tau2hat2*(Vtt-Vt%*%Vi%*%t(Vt))
    yi= ytest-mup2_shift
    ldetK <- determinant(Sigmap2, logarithm=TRUE)$modulus
    ssq = t(yi) %*% pinv(Sigmap2) %*% yi
    ll <- - (1/2 * ssq) - (1/2)*ldetK - k/2*log(2*pi)

    mup2_shift <- c(mup2_shift,ll)
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

    ## derive test sets inverse covariance kernel
    Ahat = (ahat^2)*deltatest+1
    Dtest <- plgp::distance(ttest)
    Ctt <- exp(-bhat*Dtest/Ahat) ## covariance matrix
    Vtt = (1/sqrt(Ahat)*Ctt) %>% as.matrix()

    deltaoff = Matrix(1-bdiag(replicate(2,matrix(1,ntest,n),simplify=FALSE)),sparse = TRUE)
    Ahatoff = (ahat^2)*deltaoff+1
    Dtestoff <- plgp::distance(ttesthat,c(t1,t2hat))
    Ct <- exp(-bhat*Dtestoff/Ahatoff) ## covariance matrix
    Vt = (1/sqrt(Ahatoff)*Ct) %>% as.matrix()
    mup2_shift <- Vt %*% Vi %*% Y

    # prediction variance
    Sigmap2 <- tau2hat2*(Vtt-Vt%*%Vi%*%t(Vt))
    yi= ytest-mup2_shift
    ldetK <- determinant(Sigmap2, logarithm=TRUE)$modulus
    ssq = t(yi) %*% pinv(Sigmap2) %*% yi
    ll <- - (1/2 * ssq) - (1/2)*ldetK - k/2*log(2*pi)

    mup2_shift <- c(mup2_shift,ll)

  }
  if(kernel == "sep.matern"){
    bhat1 <- params[1]
    bhat2 <- params[2]
    ytest1 <- ytest[1:n]
    ytest2 <- ytest[(n+1):(2*n)]

    v = 1/2
    D1 <-  sqrt(plgp::distance(t1))
    D1[D1 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat1*D1)^v
    part2 = 2^(1-v) / gamma(v)
    part3 = besselK(sqrt(2*v) * bhat1*D1, nu = v)
    Vk1 = (part1 * part2 *part3) %>% as.matrix()

    Vihat1 <- pinv(Vk1) %>% as.matrix()

    ## derive test sets inverse covariance kernel
    Dtest1 <- sqrt(plgp::distance(ttest[1:n]))
    Dtest1[Dtest1 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat1*Dtest1)^v
    part3 = besselK(sqrt(2*v) * bhat1*Dtest1, nu = v)
    Vttest1 = (part1 * part2 *part3) %>% as.matrix()

    Dtestoff1 <- sqrt(plgp::distance(ttest[1:n],t1))
    Dtestoff1[Dtestoff1 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat1*Dtestoff1)^v
    part3 = besselK(sqrt(2*v) * bhat1*Dtestoff1, nu = v)
    Vkoff1 = (part1 * part2 *part3) %>% as.matrix()

    mup2_shift1 <- Vkoff1 %*% Vihat1 %*% Y[1:n]

    # prediction variance
    Sigmap1 <- Vttest1-Vkoff1%*%Vihat1%*%t(Vkoff1)
    yi= ytest1-mup2_shift1
    ldetK1 <- determinant(Sigmap1, logarithm=TRUE)$modulus
    ssq = t(yi) %*% pinv(Sigmap1) %*% yi
    # ll1 <- - (1/2 * ssq) - (1/2)*ldetK - k/2*log(2*pi)

    ll1 <- - (n/2)*log(ssq /n) - (1/2)*ldetK1 -n/2 - n/2*log(2*pi)

    ## Another group
    D2 <-  sqrt(plgp::distance(t2))
    D2[D2 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat2*D2)^v
    part3 = besselK(sqrt(2*v) * bhat2*D2, nu = v)
    Vk2 = (part1 * part2 *part3) %>% as.matrix()
    Vihat2 <- pinv(Vk2) %>% as.matrix()

    ## derive test sets inverse covariance kernel
    Dtest2 <- sqrt(plgp::distance(ttest[(n+1):(2*n)]))
    Dtest2[Dtest2 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat1*Dtest2)^v
    part3 = besselK(sqrt(2*v) * bhat1*Dtest2, nu = v)
    Vttest2 = (part1 * part2 *part3) %>% as.matrix()

    Dtestoff2 <- sqrt(plgp::distance(ttest[(n+1):(2*n)],t2))
    Dtestoff2[Dtestoff2 == 0] <- 1e-8
    part1 = (sqrt(2*v) * bhat2*Dtestoff2)^v
    part3 = besselK(sqrt(2*v) * bhat2*Dtestoff2, nu = v)
    Vkoff2 = (part1 * part2 *part3) %>% as.matrix()
    mup2_shift2 <- Vkoff2 %*% Vihat2 %*% Y[(n+1):(2*n)]

    # prediction variance
    Sigmap2 <- Vttest2-Vkoff2%*%Vihat2%*%t(Vkoff2)
    yi= ytest2-mup2_shift2
    ldetK2 <- determinant(Sigmap2, logarithm=TRUE)$modulus
    ssq = t(yi) %*% pinv(Sigmap2) %*% yi
    # ll1 <- - (1/2 * ssq) - (1/2)*ldetK - k/2*log(2*pi)

    ll2 <- - (n/2)*log(ssq /n) - (1/2)*ldetK2 -n/2 - n/2*log(2*pi)

    mup2_shift <- c(mup2_shift1,mup2_shift2, ll1+ll2)
  }
  if(kernel == "sep.rbf"){
    bhat1 <- params[1]
    bhat2 <- params[2]
    ytest1 <- ytest[1:ntest]
    ytest2 <- ytest[(ntest+1):(2*ntest)]

    D1 <-  plgp::distance(t1)
    Chat1 <- exp(-bhat1*D1) ## covariance matrix
    Vihat1 <- pinv(Chat1) %>% as.matrix()

    ## derive test sets inverse covariance kernel
    Dtest1 <- plgp::distance(ttest[1:ntest])
    Vttest1 = exp(-bhat1*Dtest1)

    Dtestoff1 <- plgp::distance(ttest[1:ntest],t1)
    Vkoff1 = exp(-bhat1*Dtestoff1)

    mup2_shift1 <- Vkoff1 %*% Vihat1 %*% Y[1:n]

    # prediction variance
    Sigmap1 <- Vttest1-Vkoff1%*%Vihat1%*%t(Vkoff1)
    yi= ytest1-mup2_shift1
    ldetK1 <- determinant(Sigmap1, logarithm=TRUE)$modulus
    ssq = t(yi) %*% pinv(Sigmap1) %*% yi
    # ll1 <- - (1/2 * ssq) - (1/2)*ldetK - k/2*log(2*pi)

    ll1 <- - (n/2)*log(ssq /n) - (1/2)*ldetK1 -n/2 - n/2*log(2*pi)

    ## Another group
    D2 <-  plgp::distance(t2)
    Chat2 <- exp(-bhat1*D2) ## covariance matrix
    Vihat2 <- pinv(Chat2) %>% as.matrix()

    ## derive test sets inverse covariance kernel
    Dtest2 <- plgp::distance(ttest[(ntest+1):(2*ntest)])
    Vttest2 = exp(-bhat2*Dtest2)

    Dtestoff2 <- plgp::distance(ttest[(ntest+1):(2*ntest)],t2)
    Vkoff2 = exp(-bhat2*Dtestoff2)

    mup2_shift2 <- Vkoff2 %*% Vihat2 %*% Y[(n+1):(2*n)]

    # prediction variance
    Sigmap2 <- Vttest2-Vkoff2%*%Vihat2%*%t(Vkoff2)
    yi= ytest2-mup2_shift2
    ldetK2 <- determinant(Sigmap2, logarithm=TRUE)$modulus
    ssq = t(yi) %*% pinv(Sigmap2) %*% yi
    # ll1 <- - (1/2 * ssq) - (1/2)*ldetK - k/2*log(2*pi)

    ll2 <- - (n/2)*log(ssq /n) - (1/2)*ldetK2 -n/2 - n/2*log(2*pi)

    mup2_shift <- c(mup2_shift1,mup2_shift2, ll1+ll2)
  }

  return(mup2_shift)
}

simulateMultiData <- function(b,a,s, tau2, n, kernel, timemax = 25, fixrange = FALSE, nrep = 20){
  increment =1
  if(isTRUE(fixrange)){
    t1_tilde <- matrix(seq(0, timemax, length.out = n), ncol=1)
  }else{
    t1_tilde <- matrix(seq(0, 50, length.out = n), ncol=1)
  }
  t1 <- t1_tilde + runif(n,-increment/4,increment/4)
  # concatenate three time series time together
  t <- rep(t1,3)
  group.index <- rep(c(1,2,3),each=n)
  group.size <- rle(group.index)
  group.size$values <- s
  s.expand <- inverse.rle(group.size)
  # add time shift to the second time series
  tinput <- t+s.expand
  # create group matrix
  delta = Matrix(1-bdiag(replicate(3,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
  A = delta+1
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
  ## generate 100 replicates under a prior with MVN(0,sigma)
  Y <- rmvnorm(nrep, sigma= sigma %>% as.matrix()) %>% t()
  rownames(Y) <- tinput - s.expand
  return(Y)
}

nlMultiexp <- function(par, t, Y, delta, group.size) # group.size is rle of group.index
{
  ngroup <- length(group.size$values)
  b <- softplus(par[1])
  a12 <- softplus(par[2])
  a13 <- softplus(par[3])
  a23 <- softplus(par[4])
  s <- par[5:(3+ngroup)]
  k <- length(Y)
  group.size$values <- c(0,s)
  s.expand <- inverse.rle(group.size) # input t vector
  tinput <- t+s.expand
  n = k /ngroup
  # Initialize delta2 as a copy of delta
  delta2 <- delta

  # Apply the transformations
  delta2[1:n, (n+1):(2*n)] <- delta[1:n, (n+1):(2*n)] * a12
  delta2[(n+1):(2*n), 1:n] <- delta[(n+1):(2*n), 1:n] * a12 # Symmetric part

  delta2[1:n, (2*n+1):(3*n)] <- delta[1:n, (2*n+1):(3*n)] * a13
  delta2[(2*n+1):(3*n), 1:n] <- delta[(2*n+1):(3*n), 1:n] * a13 # Symmetric part

  delta2[(n+1):(2*n), (2*n+1):(3*n)] <- delta[(n+1):(2*n), (2*n+1):(3*n)] * a23
  delta2[(2*n+1):(3*n), (n+1):(2*n)] <- delta[(2*n+1):(3*n), (n+1):(2*n)] * a23 # Symmetric part

  A = delta2+1
  D <- sqrt(plgp::distance(tinput))
  C <- exp(-b*D) ## covariance matrix
  Vk = (1/A*C) %>% as.matrix()
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}
  if (a12 + a13 < a23 || a12 + a23 < a13 || a13 + a23 < a12) return(-1e6)
  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y / k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a, -ll, "\n")
  return(-ll)
}

nlMultirbf <- function(par, t, Y, delta, group.size) # group.size is rle of group.index
{
  ngroup <- length(group.size$values)
  b <- softplus(par[1])
  a12 <- softplus(par[2])
  a13 <- softplus(par[3])
  a23 <- softplus(par[4])
  s <- par[5:(3+ngroup)]
  k <- length(Y)
  group.size$values <- c(0,s)
  s.expand <- inverse.rle(group.size) # input t vector
  tinput <- t+s.expand
  n = k /ngroup
  # Initialize delta2 as a copy of delta
  delta2 <- delta

  # Apply the transformations
  delta2[1:n, (n+1):(2*n)] <- delta[1:n, (n+1):(2*n)] * a12
  delta2[(n+1):(2*n), 1:n] <- delta[(n+1):(2*n), 1:n] * a12 # Symmetric part

  delta2[1:n, (2*n+1):(3*n)] <- delta[1:n, (2*n+1):(3*n)] * a13
  delta2[(2*n+1):(3*n), 1:n] <- delta[(2*n+1):(3*n), 1:n] * a13 # Symmetric part

  delta2[(n+1):(2*n), (2*n+1):(3*n)] <- delta[(n+1):(2*n), (2*n+1):(3*n)] * a23
  delta2[(2*n+1):(3*n), (n+1):(2*n)] <- delta[(2*n+1):(3*n), (n+1):(2*n)] * a23 # Symmetric part

  A = delta2+1
  D <- plgp::distance(tinput)/A
  C <- exp(-b*D) ## covariance matrix
  Vk = (1/sqrt(A)*C) %>% as.matrix()
  Vi <- pinv(Vk)
  ldetK <- determinant(Vk, logarithm=TRUE)$modulus
  if(is.infinite(ldetK)){ldetK=-1e6}
  if (a12 + a13 < a23 || a12 + a23 < a13 || a13 + a23 < a12) return(-1e6)
  ll <- - (k/2)*log(t(Y) %*% Vi %*% Y / k) - (1/2)*ldetK
  # counter <<- counter + 1
  # cat(b,a12,a13,a23,s, -ll, "\n")
  return(-ll)
}

deriveMultiEstimation <-  function(Y, kernel, sl = 0, su = 5){
  k = nrow(Y)
  res = matrix(0,ncol = 8,nrow = ncol(Y))
  n = nrow(Y)/3
  t = rownames(Y)[1:n] %>% as.numeric()
  group.index <- rep(c(1,2,3),each=n)
  group.size <- rle(group.index)
  delta = Matrix(1-bdiag(replicate(3,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
  if(kernel == "exp"){
    # perform LExp parameter estimation
    for (i in 1:ncol(Y)) {
      yi = Y[,i]
      out <- optim(c(1,1,1,1,1,3), nlMultiexp, method="L-BFGS-B",lower=c(-2,rep(-2,3),rep(sl,2)),
                   upper=c(10,rep(2,3),rep(su,2)),
                   t= t, Y=yi, delta=delta,group.size = group.size)
      bhat = softplus(out$par[1])
      a12hat = softplus(out$par[2])
      a13hat = softplus(out$par[3])
      a23hat = softplus(out$par[4])
      shat = out$par[5:6]
      # MLE estimates tau2hat
      group.size$values <- c(0,shat)
      s.expand <- inverse.rle(group.size) # input t vector
      that <- t+s.expand
      delta2 <- delta

      # Apply the transformations
      delta2[1:n, (n+1):(2*n)] <- delta[1:n, (n+1):(2*n)] * a12hat
      delta2[(n+1):(2*n), 1:n] <- delta[(n+1):(2*n), 1:n] * a12hat # Symmetric part

      delta2[1:n, (2*n+1):(3*n)] <- delta[1:n, (2*n+1):(3*n)] * a13hat
      delta2[(2*n+1):(3*n), 1:n] <- delta[(2*n+1):(3*n), 1:n] * a13hat # Symmetric part

      delta2[(n+1):(2*n), (2*n+1):(3*n)] <- delta[(n+1):(2*n), (2*n+1):(3*n)] * a23hat
      delta2[(2*n+1):(3*n), (n+1):(2*n)] <- delta[(2*n+1):(3*n), (n+1):(2*n)] * a23hat # Symmetric part

      Ahat = delta2+1
      Dhat <-  sqrt(plgp::distance(that))
      Chat <- exp(-bhat*Dhat) ## covariance matrix
      Vkhat = (1/Ahat*Chat) %>% as.matrix()

      Vihat <- pinv(Vkhat)
      ldetK <- determinant(Vkhat, logarithm=TRUE)$modulus
      tau2hat2 <- drop(t(yi) %*% Vihat %*% yi / length(that))
      ll <- - (k/2)*log(t(yi) %*% Vihat %*% yi / k) - (1/2)*ldetK - k/2

      res[i,]<- c(bhat,a12hat,a13hat,a23hat,shat,tau2hat2, ll)
    }
  }
  if(kernel == "rbf"){
    # perform LRBF parameter estimation
    t = rownames(Y) %>% as.numeric()
    for (i in 1:ncol(Y)) {
      yi = Y[,i]
      out <- optim(c(1,1,1,1,1,3), nlMultirbf, method="L-BFGS-B",lower=c(-2,rep(-2,3),rep(sl,2)),
                   upper=c(10,rep(2,3),rep(su,2)),
                   t= t, Y=yi, delta=delta,group.size = group.size)
      bhat = softplus(out$par[1])
      a12hat = softplus(out$par[2])
      a13hat = softplus(out$par[3])
      a23hat = softplus(out$par[4])
      shat = out$par[5:6]
      # MLE estimates tau2hat
      group.size$values <- c(0,shat)
      s.expand <- inverse.rle(group.size) # input t vector
      that <- t+s.expand
      delta2 <- delta

      # Apply the transformations
      delta2[1:n, (n+1):(2*n)] <- delta[1:n, (n+1):(2*n)] * a12hat
      delta2[(n+1):(2*n), 1:n] <- delta[(n+1):(2*n), 1:n] * a12hat # Symmetric part

      delta2[1:n, (2*n+1):(3*n)] <- delta[1:n, (2*n+1):(3*n)] * a13hat
      delta2[(2*n+1):(3*n), 1:n] <- delta[(2*n+1):(3*n), 1:n] * a13hat # Symmetric part

      delta2[(n+1):(2*n), (2*n+1):(3*n)] <- delta[(n+1):(2*n), (2*n+1):(3*n)] * a23hat
      delta2[(2*n+1):(3*n), (n+1):(2*n)] <- delta[(2*n+1):(3*n), (n+1):(2*n)] * a23hat # Symmetric part

      Ahat = delta2+1
      Dhat <- plgp::distance(that) / Ahat
      Chat <- exp(-bhat*Dhat) ## covariance matrix
      Vkhat = (1/sqrt(Ahat)*Chat) %>% as.matrix()
      Vihat <- pinv(Vkhat)
      ldetK <- determinant(Vkhat, logarithm=TRUE)$modulus
      tau2hat2 <- drop(t(yi) %*% Vihat %*% yi / length(that))
      ll <- - (k/2)*log(t(yi) %*% Vihat %*% yi / k) - (1/2)*ldetK -k/2
      res[i,]<- c(bhat,a12hat,a13hat,a23hat,shat,tau2hat2, ll)
    }
  }
  return(res)
}

# helper function for simulation prediction which only need posterior mean
interpolation_Multisim <- function(t,Y,ytest,params,ttest, kernel = "lexp", prediction = TRUE){
  ntest = length(ttest)/3
  # tt2 <- tt+shat # use shat=0 if directly model the predicted value after shift

  ## derive original inverse covariance kernel
  n = length(Y)/3
  t1 <- t[1:n]
  t2 <- t[(n+1):(3*n)]
  if(kernel == "lexp"){
    bhat <- params[1]
    a12hat <- params[2]
    a13hat <- params[3]
    a23hat <- params[4]
    shat <- params[5:6]
    tau2hat2 <- params[7]

    t2hat <- t2+shat
    delta = Matrix(1-bdiag(replicate(3,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
    deltatest = Matrix(1-bdiag(replicate(3,matrix(1,ntest,ntest),simplify=FALSE)),sparse = TRUE)
    group.index.test <- rep(c(0,shat),each=ntest)
    if(isTRUE(prediction)){
      ttesthat <- ttest + group.index.test
    }else{
      ttesthat <- ttest
    }
    ## derive original inverse covariance kernel
    delta2 <- delta

    # Apply the transformations
    delta2[1:n, (n+1):(2*n)] <- delta[1:n, (n+1):(2*n)] * a12hat
    delta2[(n+1):(2*n), 1:n] <- delta[(n+1):(2*n), 1:n] * a12hat # Symmetric part

    delta2[1:n, (2*n+1):(3*n)] <- delta[1:n, (2*n+1):(3*n)] * a13hat
    delta2[(2*n+1):(3*n), 1:n] <- delta[(2*n+1):(3*n), 1:n] * a13hat # Symmetric part

    delta2[(n+1):(2*n), (2*n+1):(3*n)] <- delta[(n+1):(2*n), (2*n+1):(3*n)] * a23hat
    delta2[(2*n+1):(3*n), (n+1):(2*n)] <- delta[(2*n+1):(3*n), (n+1):(2*n)] * a23hat # Symmetric part

    Ahat = delta2+1
    D <-  sqrt(plgp::distance(c(t1,t2hat)))
    Chat <- exp(-bhat*D)
    Vk = (1/Ahat*Chat) %>% as.matrix() ## covariance matrix
    Vi <- pinv(Vk) ## psduo-inverse covariance matrix

    ## derive test sets inverse covariance kernel
    delta2test <- deltatest

    # Apply the transformations
    delta2test[1:n, (n+1):(2*n)] <- deltatest[1:n, (n+1):(2*n)] * a12hat
    delta2test[(n+1):(2*n), 1:n] <- deltatest[(n+1):(2*n), 1:n] * a12hat # Symmetric part

    delta2test[1:n, (2*n+1):(3*n)] <- deltatest[1:n, (2*n+1):(3*n)] * a13hat
    delta2test[(2*n+1):(3*n), 1:n] <- deltatest[(2*n+1):(3*n), 1:n] * a13hat # Symmetric part

    delta2test[(n+1):(2*n), (2*n+1):(3*n)] <- deltatest[(n+1):(2*n), (2*n+1):(3*n)] * a23hat
    delta2test[(2*n+1):(3*n), (n+1):(2*n)] <- deltatest[(2*n+1):(3*n), (n+1):(2*n)] * a23hat # Symmetric part

    Ahattest = delta2test+1
    Dtest <- sqrt(plgp::distance(ttest))
    Ctt <- exp(-bhat*Dtest) ## covariance matrix
    Vtt = (1/Ahattest*Ctt) %>% as.matrix()

    ## derive off diagonal inverse covariance kernel
    deltaoff = Matrix(1-bdiag(replicate(3,matrix(1,ntest,n),simplify=FALSE)),sparse = TRUE)
    delta2off <- deltaoff

    # Apply the transformations
    delta2off[1:n, (n+1):(2*n)] <- deltaoff[1:n, (n+1):(2*n)] * a12hat
    delta2off[(n+1):(2*n), 1:n] <- deltaoff[(n+1):(2*n), 1:n] * a12hat # Symmetric part

    delta2off[1:n, (2*n+1):(3*n)] <- deltaoff[1:n, (2*n+1):(3*n)] * a13hat
    delta2off[(2*n+1):(3*n), 1:n] <- deltaoff[(2*n+1):(3*n), 1:n] * a13hat # Symmetric part

    delta2off[(n+1):(2*n), (2*n+1):(3*n)] <- deltaoff[(n+1):(2*n), (2*n+1):(3*n)] * a23hat
    delta2off[(2*n+1):(3*n), (n+1):(2*n)] <- deltaoff[(2*n+1):(3*n), (n+1):(2*n)] * a23hat # Symmetric part

    Ahatoff = delta2off+1

    Dtestoff <- sqrt(plgp::distance(ttesthat,c(t1,t2hat)))
    Ct <- exp(-bhat*Dtestoff) ## covariance matrix
    Vt = (1/Ahatoff*Ct) %>% as.matrix()

    mup2_shift <- Vt %*% Vi %*% Y
    # prediction variance
    Sigmap2 <- tau2hat2*(Vtt-Vt%*%Vi%*%t(Vt))
    yi= ytest-mup2_shift
    ldetK <- determinant(Sigmap2, logarithm=TRUE)$modulus
    ssq = t(yi) %*% pinv(Sigmap2) %*% yi
    ll <- - (1/2 * ssq) - (1/2)*ldetK - k/2*log(2*pi)

    mup2_shift <- c(mup2_shift,ll)
  }
  if(kernel == "lrbf"){
    bhat <- params[1]
    a12hat <- params[2]
    a13hat <- params[3]
    a23hat <- params[4]
    shat <- params[5:6]
    tau2hat2 <- params[7]

    t2hat <- t2+shat
    delta = Matrix(1-bdiag(replicate(3,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
    deltatest = Matrix(1-bdiag(replicate(3,matrix(1,ntest,ntest),simplify=FALSE)),sparse = TRUE)
    group.index.test <- rep(c(0,shat),each=ntest)

    if(isTRUE(prediction)){
      ttesthat <- ttest + group.index.test
    }else{
      ttesthat <- ttest
    }

    delta2 <- delta

    # Apply the transformations
    delta2[1:n, (n+1):(2*n)] <- delta[1:n, (n+1):(2*n)] * a12hat
    delta2[(n+1):(2*n), 1:n] <- delta[(n+1):(2*n), 1:n] * a12hat # Symmetric part

    delta2[1:n, (2*n+1):(3*n)] <- delta[1:n, (2*n+1):(3*n)] * a13hat
    delta2[(2*n+1):(3*n), 1:n] <- delta[(2*n+1):(3*n), 1:n] * a13hat # Symmetric part

    delta2[(n+1):(2*n), (2*n+1):(3*n)] <- delta[(n+1):(2*n), (2*n+1):(3*n)] * a23hat
    delta2[(2*n+1):(3*n), (n+1):(2*n)] <- delta[(2*n+1):(3*n), (n+1):(2*n)] * a23hat # Symmetric part

    Ahat = delta2+1
    D <-  plgp::distance(c(t1,t2hat))
    Chat <- exp(-bhat*D/Ahat)
    Vk = (1/sqrt(Ahat)*Chat) %>% as.matrix() ## covariance matrix
    Vi <- pinv(Vk) ## psduo-inverse covariance matrix

    deltaoff = Matrix(1-bdiag(replicate(3,matrix(1,ntest,n),simplify=FALSE)),sparse = TRUE)
    delta2off <- deltaoff

    # Apply the transformations
    delta2off[1:n, (n+1):(2*n)] <- deltaoff[1:n, (n+1):(2*n)] * a12hat
    delta2off[(n+1):(2*n), 1:n] <- deltaoff[(n+1):(2*n), 1:n] * a12hat # Symmetric part

    delta2off[1:n, (2*n+1):(3*n)] <- deltaoff[1:n, (2*n+1):(3*n)] * a13hat
    delta2off[(2*n+1):(3*n), 1:n] <- deltaoff[(2*n+1):(3*n), 1:n] * a13hat # Symmetric part

    delta2off[(n+1):(2*n), (2*n+1):(3*n)] <- deltaoff[(n+1):(2*n), (2*n+1):(3*n)] * a23hat
    delta2off[(2*n+1):(3*n), (n+1):(2*n)] <- deltaoff[(2*n+1):(3*n), (n+1):(2*n)] * a23hat # Symmetric part

    Ahatoff = delta2off+1
    Dtestoff <- plgp::distance(ttesthat,c(t1,t2hat))
    Ct <- exp(-bhat*Dtestoff/Ahatoff) ## covariance matrix
    Vt = (1/sqrt(Ahatoff)*Ct) %>% as.matrix()

    ## derive test sets inverse covariance kernel
    delta2test <- deltatest

    # Apply the transformations
    delta2test[1:n, (n+1):(2*n)] <- deltatest[1:n, (n+1):(2*n)] * a12hat
    delta2test[(n+1):(2*n), 1:n] <- deltatest[(n+1):(2*n), 1:n] * a12hat # Symmetric part

    delta2test[1:n, (2*n+1):(3*n)] <- deltatest[1:n, (2*n+1):(3*n)] * a13hat
    delta2test[(2*n+1):(3*n), 1:n] <- deltatest[(2*n+1):(3*n), 1:n] * a13hat # Symmetric part

    delta2test[(n+1):(2*n), (2*n+1):(3*n)] <- deltatest[(n+1):(2*n), (2*n+1):(3*n)] * a23hat
    delta2test[(2*n+1):(3*n), (n+1):(2*n)] <- deltatest[(2*n+1):(3*n), (n+1):(2*n)] * a23hat # Symmetric part

    Ahattest = delta2test+1
    Dtest <- plgp::distance(ttesthat)
    Ctt <- exp(-bhat*Dtest/Ahattest) ## covariance matrix
    Vtt = (1/sqrt(Ahattest)*Ctt) %>% as.matrix()

    mup2_shift <- Vt %*% Vi %*% Y
    # prediction variance
    Sigmap2 <- tau2hat2*(Vtt-Vt%*%Vi%*%t(Vt))
    yi= ytest-mup2_shift
    ldetK <- determinant(Sigmap2, logarithm=TRUE)$modulus
    ssq = t(yi) %*% pinv(Sigmap2) %*% yi
    ll <- - (1/2 * ssq) - (1/2)*ldetK - k/2*log(2*pi)

    mup2_shift <- c(mup2_shift,ll)
  }

  return(mup2_shift)
}
