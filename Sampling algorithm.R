library(mvtnorm)
Adaptive_MALA <- function(dlogpi, logpi,nits, h, x_curr,delta = 0.7,target_a = 0.574) {
  
  logpi_curr <- logpi(x_curr)
  dlogpi_curr <- dlogpi(x_curr)
  d <- length(x_curr)
  accepted <- 0
  x_store <- matrix(nrow = nits, ncol = d)
  I_d <- diag(rep(1,d))
  V <- I_d
  mu <- x_curr

  for (i in 1:nits) {
    # propose a candidate move
    x_prop <- x_curr + c(h*V%*%dlogpi_curr/2)+c(sqrt(h)*rmvnorm(1,sigma = V))
    logpi_prop <- logpi(x_prop)
    dlogpi_prop <- dlogpi(x_prop)
    
    # accept-reject
    logq_pgivenc <- dmvnorm(x_prop,mean = x_curr+h*c(V%*%dlogpi_curr)/2,sigma = h*V, log=T)
    logq_cgivenp <- dmvnorm(x_curr,mean = x_prop+h*c(V%*%dlogpi_prop)/2,sigma = h*V, log=T)
    loga <- logpi_prop+logq_cgivenp - logpi_curr-logq_pgivenc
    u <- runif(1)
    if ( is.finite(loga) && log(u) < loga) {
      x_curr <- x_prop
      logpi_curr <- logpi_prop
      dlogpi_curr <- dlogpi_prop
      accepted <- accepted + 1
    }
    # store x_curr
    x_store[i, ] <- x_curr
    
    # set the learning rate
    gamma <- 1/(i)^delta
    # update the covariance matrix
    V <- V + gamma * ((x_curr - mu) %*% t(x_curr - mu) - V)
    V <- V + 1e-6 * I_d 
    # update the mean
    mu <- mu + gamma * (x_curr - mu)
    # update the step size
    alpha <-  ifelse(!is.nan(loga),min(1, exp(loga)),target_a)
    h <- h * exp(gamma * (alpha - target_a))
    
    if (floor(i/10000)==i/10000) { cat(i," iterations completed.","\n") }
  }
  
  return(list(x_store = x_store, a_rate = accepted/nits,V = V,mu = mu,step_size = h))
}

Adaptive_RWM <- function(logpi, nits, h, x_curr,delta = 0.7,target_a = 0.23,dlogpi = NULL){
  accepted <- 0
  logpi_curr <- logpi(x_curr)
  d <- length(x_curr)
  x_store <- matrix(0,nrow = nits, ncol = d)
  V = diag(rep(1,d))
  I_d <- V
  mu <- rnorm(d)
  
  for(i in 1:nits){
    
    x_prop <- x_curr + h * as.vector(rmvnorm(1, sigma =  V))
    logpi_prop <- logpi(x_prop)
    loga <- logpi_prop - logpi_curr
    
    if (is.finite(logpi_prop) && log(runif(1)) < loga) {
      x_curr <- x_prop
      logpi_curr <- logpi_prop
      accepted <- accepted + 1
    }
    
    x_store[i, ] <- x_curr
    
    gamma <- 1/(i)^delta
    V <- V + gamma * ((x_curr - mu) %*% t(x_curr - mu) - V)
    V <- V + 1e-6 * I_d 
    mu <- mu + gamma * (x_curr - mu)
    alpha <- min(1, exp(loga))
    h <- h * sqrt(exp(gamma * (alpha - target_a)))
    
    if(i %% 10000 == 0){
      cat("iteration",i)
    }
    
  }
  return(list(x_store = x_store, a_rate = accepted/nits,
              V = V,mu = mu,step_size = h))
}

newadaptiveMCMC <- function(logpi, nits, h, x_curr, delta = 0.7, target_a = 0.23,chain_len = 5){
  library(mvtnorm)
  accepted <- 0
  d <- length(x_curr)
  x_store <- matrix(0, nrow = nits, ncol = d)
  v_store <- sort(rnorm(nits))
  V <- diag(d)
  I_d <- diag(d)
  mu <- x_curr
  
  for(i in 1:nits){
    
    v <- v_store[i]
    logpi_curr <- logpi(x_curr, v)
    
    for (j in 1:chain_len) {
      x_prop <- x_curr + h * as.vector(rmvnorm(1, sigma = V))
      logpi_prop <- logpi(x_prop, v)
      
      if (!is.finite(logpi_prop)) next  # skip the abnormal value
      
      loga <- logpi_prop - logpi_curr
      
      if (log(runif(1)) < loga) {
        x_curr <- x_prop
        logpi_curr <- logpi_prop
        accepted <- accepted + 1
      }
    }
    
    x_store[i,] <- x_curr
    
    gamma <- 1/(i)^delta
    V <- V + gamma * ((x_curr - mu) %*% t(x_curr - mu) - V)
    V <- V + 1e-6 * I_d 
    mu <- mu + gamma * (x_curr - mu)
    alpha <- min(1, exp(loga))
    h <- h * sqrt(exp(gamma * (alpha - target_a)))
    
    if(i %% 10000 == 0){ cat("iteration",i)}
    
  }
  
  return(list(x_store = x_store, v_store = v_store,
              a_rate = accepted / (nits * chain_len),
              V = V, mu = colMeans(x_store),
              step_size  = h))
}

