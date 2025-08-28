# define a function to generate the samples from skew-normal
generate_sn <- function(nits,xi,omega,alpha){
  delta  = alpha / sqrt(1+alpha^2)
  return(xi + omega * (delta* abs(rnorm(nits)) + sqrt(1 - delta^2) *rnorm(nits)))
}

# three different method to fit the parameter
# using moment in high dimension to fit parameter, the final output Omega is a matrix with non-zero 
# element in non-diag position

# not component-wise fit a joint distribution
fit_highdim_sn <- function(samples) {
  if(!is.matrix(samples)){
    stop("The input must be a matrix")
  }
  n <- nrow(samples)
  d <- ncol(samples)
  
  if(d > 1){
    sample_mean <- colMeans(samples)
    sample_cov <- cov(samples)
    third_moments <- sapply(1:d, function(j) 
      mean((samples[,j] - sample_mean[j])^3) / var(samples[,j])^1.5)
  } else {
    sample_mean <- mean(samples)
    sample_cov <- var(samples)
    third_moments <- mean((samples - sample_mean)^3) / sample_cov^1.5
  }
  
  # 初始值
  xi_init <- sample_mean
  Omega_init <- sample_cov[lower.tri(sample_cov, diag = TRUE)]
  alpha_init <- rep(0.1, d)
  start <- c(xi_init, Omega_init, alpha_init)
  
  make_positive_definite <- function(mat, max_attempts = 100, epsilon = 1e-6) {
    attempt <- 0
    while (attempt < max_attempts) {
      # Try Cholesky to check if PD
      result <- tryCatch({
        chol(mat)
        TRUE
      }, error = function(e) FALSE)
      
      if (result) break
      mat <- mat + epsilon * diag(ncol(mat))
      attempt <- attempt + 1
      epsilon <- epsilon * 1.6
    }
    return(mat)
  }
  
  moment_equations <- function(params) {
    pi <- 3.1415926535
    xi <- params[1:d]
    Omega_vec <- params[(d+1):(d + d*(d+1)/2)]
    Omega <- matrix(0, d, d)
    Omega[lower.tri(Omega, diag = TRUE)] <- Omega_vec
    Omega <- Omega + t(Omega) - diag(diag(Omega),ncol = d)  # 保证对称性
    
    # make sure Omgea is positive definite
    Omega <- make_positive_definite(Omega)
    
    alpha <- params[(d + d*(d+1)/2 + 1):(d + d*(d+1)/2 + d)]
    
    omega_diag <- sqrt(diag(Omega))
    omega <- diag(omega_diag, nrow = d)
    Omega_bar <- solve(omega) %*% Omega %*% solve(omega)
    denom <- sqrt(1 + t(alpha) %*% Omega_bar %*% alpha)
    delta <- as.vector(Omega_bar %*% alpha / as.numeric(denom))
    
    # 均值
    mean_theory <- xi + sqrt(2/pi) * omega_diag * delta
    
    # 协方差矩阵
    Sigma_theory <- Omega - (2/pi) * (omega %*% (delta %*% t(delta)) %*% omega)
    
    # 三阶中心矩
    alpha_marginal <- rep(0, d)
    for (j in 1:d) {
      denom_j <- 1 + sum((alpha[-j]^2) * diag(Omega)[-j] / diag(Omega)[j])
      alpha_marginal[j] <- alpha[j] / sqrt(denom_j)
    }
    delta_marginal <- alpha_marginal / sqrt(1 + alpha_marginal^2)
    mu_z_marginal <- sqrt(2/pi) * delta_marginal
    third_theory <- ((4 - pi)/2) * (mu_z_marginal^3) / (1 - mu_z_marginal^2)^1.5
    
    weight <- 1e5
    c(
      (mean_theory - sample_mean) * weight,
      (as.vector(Sigma_theory[lower.tri(Sigma_theory, diag = TRUE)] - 
                   sample_cov[lower.tri(sample_cov, diag = TRUE)])) * weight,
      (third_theory - third_moments) * 100 * weight
    )
  }
  
  res <- nleqslv(start, moment_equations, method = "Broyden", control = list(ftol = 1e-8))
  
  xi <- res$x[1:d]
  Omega_vec <- res$x[(d+1):(d + d*(d+1)/2)]
  Omega <- matrix(0, d, d)
  Omega[lower.tri(Omega, diag = TRUE)] <- Omega_vec
  Omega <- Omega + t(Omega) - diag(diag(Omega),ncol = d)
  Omega <- make_positive_definite(Omega)
  
  alpha <- res$x[(d + d*(d+1)/2 + 1):(d + d*(d+1)/2 + d)]
  omega_diag <- sqrt(diag(Omega))
  omega <- diag(omega_diag, nrow = d)
  Omega_bar <- solve(omega) %*% Omega %*% solve(omega)
  denom <- sqrt(1 + t(alpha) %*% Omega_bar %*% alpha)
  delta <- as.vector(Omega_bar %*% alpha / as.numeric(denom))
  
  list(xi = xi, Omega = Omega, alpha = alpha, delta = delta, gap = res$fvec)
}
# Transformation map with Omega is a matrix 
transformation_highdim_map_sn <- function(xi,Omega,alpha,log_pi,nits,
                                          chain_len = 5,target_a= 0.23){
  
  d <- length(xi) # get the dimension
  omega_diag <- sqrt(diag(Omega)) 
  omega <- diag(omega_diag,nrow = d)  # transfer omega_diag into a matrix omega
  Omega_bar <- solve(omega) %*% Omega %*% solve(omega) # calculate Omega_bar
  denom <- sqrt(1 + t(alpha) %*% Omega_bar %*% alpha) 
  delta <- as.vector(Omega_bar %*% alpha / as.numeric(denom)) #calculate delta
  I_d <- diag(d)
  Delta <- diag(delta,nrow = d) # get the matrix Delta
  
  # Get A,B,C
  A <- omega%*%sqrt(I_d - Delta^2)
  B <- omega%*% Delta
  C <- xi
  
  log_d_tildeu <- function(y,v){
    x <- A %*% y + B %*% rep(abs(v),d) + C
    logf <- log_pi(x) - v^2/2
    logf <- max(logf,-1e30)
    return(logf)
  }
  
  # Sampling from log_d_tildeu
  chain<- newadaptiveMCMC(log_d_tildeu,nits = nits,h = 0.5,x_curr =rep(0,d),
                          target_a = target_a,chain_len = chain_len)
  
  # Get the samples of v and y
  samples_y <- chain$x_store
  samples_v <- chain$v_store
  
  # transport back to get the samples from original target distribution
  samples_x<-  samples_y%*%t(A) + matrix(rep(abs(samples_v),d),ncol = d) %*% t(B)+ 
    matrix(rep(C,nits),nrow = nits,ncol = d,byrow = TRUE)
  
  
  return(list(samples_x = samples_x,samples_y = samples_y,samples_v = samples_v,
              a_rate = chain$a_rate,step_size = chain$step_size))
}

# Transformation map with component-wise method and direct transformation
transformation_componentwise_direct_map_sn <- function(xi,omega,alpha,log_pi,
                                                       nits, target_a = 0.23){
  # Define the Johnson SU CDF and PDF
  d <- length(xi)
  
  logf_x <- function(x) {
    pi = 3.1415926535
    z <- (x - xi) / omega
    term1 <- -z^2/2
    term2 <- pnorm(alpha*z,log.p = TRUE)
    term3 <- log(2) - log(omega) - log(2*pi)/2
    log_pdf <- sum(term1 + term2 + term3)
    return(log_pdf) }
  
  
  # Define inverse map: x = F⁻¹(Φ(y))
  x_from_y <- function(y) {
    x <- rep(0,d)
    for(i in 1:d){
      x[i] <- qsn(pnorm(y[i]),xi = xi[i],omega = omega[i],alpha = alpha[i],solver = "RFB")
    }
    return(x)
  } 
  
  # Define the transformed density p_Y(y)
  log_p_Y <- function(y) {
    x <- x_from_y(y)
    logf <- log_pi(x) - sum(y^2/2) - logf_x(x)
    logf <- max(logf,-1e30)
    if(is.nan(logf)) logf <- -1e30
    return(logf)
  }
  
  # Sampling from log_p_Y
  chain<- Adaptive_RWM(log_p_Y,nits = nits,h = 0.5,x_curr = rep(0,d),target_a = target_a)
  
  # Get the samples of y
  samples_y <- chain$x_store
  
  # transform y back
  samples_x<- matrix(0,ncol = d,nrow = nits)
  for (i in 1:d){
    samples_x[,i] <- qsn(pnorm(samples_y[,i]),xi = xi[i],omega = omega[i],alpha = alpha[i],solver = "RFB")
  }
  return(list(samples_x = samples_x,samples_y = samples_y,
              a_rate = chain$a_rate,step_size = chain$step_size))
}


# using moment method find skew-normal approximation for each dimension
fit_sn_moment <- function(mu,variance,skew){

  sigma <- sqrt(variance)
  pi <- 3.1415926535
  est_delta <- sign(skew) * pmin(sqrt( pi / 2 * abs(skew)^(2/3) / (( abs(skew)^(2/3)  +  (2-pi/2)^(2/3) )) ),0.999999999)
  est_alpha <- est_delta / sqrt(1-est_delta^2)
  est_omega <- sigma / sqrt(max(1-2 * est_delta^2 / pi,1e-300))
  est_xi <- mu - est_omega* est_delta *sqrt(2 / pi)
  return(list(alpha = est_alpha,omega = est_omega,xi = est_xi))
}

# use MLE method to find approximation for each dimension separately
fit_sn_MLE <- function(X){
  fit_oned_ml <- function(samples) {
    
    loglik_skewnormal <- function(par, x) {
      xi <- par[1]
      omega <- par[2]
      alpha <- par[3]
      
      if (omega <= 0) return(1e30)  # 防止负尺度
      log_pdf <- sn::dsn(x, xi, omega, alpha,log = TRUE)
      return(-sum(log_pdf))
    }
    
    mu <- mean(samples)
    sigma <- sqrt(var(samples))
    skew <- skewness(samples)
    
    start <- c(1, 1, 1)
    
    result_ml <- optim(par=start, fn=loglik_skewnormal, x=samples, method="L-BFGS-B", 
                       lower=c(-Inf, 1e-6, -Inf))
    
    xi <- result_ml$par[1]
    omega <- result_ml$par[2]
    alpha <- result_ml$par[3]
    c(xi, omega, alpha)
  }
  if(!is.matrix(X)){
    return("the input must be a matrix")
  }
  r <- as.data.frame(t(apply(X,2,fit_oned_ml)))
  colnames(r) <- c("xi","omega","alpha")
  return(r)
}

# Transformation map with component-wise method and auxiliary variable v
transformation_map_sn <- function(xi,omega,alpha,log_pi,nits,target_a= 0.23,chain_len = 5){
  
  d <- length(xi)
  delta <- alpha / sqrt(1+alpha^2)
  # Get A,B,C
  A <- omega * sqrt(1 - delta^2)
  B <- omega * delta
  C <- xi
  
  log_piy <- function(y,v){
    x <- A * y + B * rep(abs(v),d) + C
    logf <- log_pi(x) - v^2/2
    logf <- max(logf,-1e30)
    return(logf)
    }
  
  # Sampling from log_piy
  chain<- newadaptiveMCMC(log_piy,nits = nits,h = 0.5,x_curr =rep(0,d),
                          target_a = target_a,chain_len = chain_len)
  
  # Get the samples of v and y
  samples_y <- chain$x_store
  samples_v <- chain$v_store

  # transport back to get the samples from original target distribution
  samples_x<-  samples_y* matrix(rep(A,nits),nrow = nits,ncol = d,byrow = TRUE) + 
    matrix(rep(abs(samples_v),d),ncol = d) * matrix(rep(B,nits),nrow = nits,ncol = d,byrow = TRUE)+ 
    matrix(rep(C,nits),nrow = nits,ncol = d,byrow = TRUE)

  return(list(samples_x = samples_x,samples_y = samples_y,samples_v = samples_v,
              a_rate = chain$a_rate,step_size = chain$step_size))
}


# Sampling using transformation map
Sampling_trsanformation_sn <- function(samples,log_pi,nits,approximation ="MLE" ,
                                        map = "indirect",target_a = 0.44){
  if(approximation == "highdim"&map == "direct"){
    return("Please change approximation method ot mapping method")}
  if(approximation == "highdim"){
    sn_approx <- fit_highdim_sn (samples)
    xi <- sn_approx$ xi
    Omega <- sn_approx$ Omega
    alpha <- sn_approx$ alpha
    result <- transformation_highdim_map_sn(xi,Omega,alpha,log_pi,nits,target_a)
    return(list(xi = xi,Omega = Omega,alpha = alpha,samples_x = result$samples_x,
                samples_y = result$samples_y, a_rate = result$a_rate, step_size = result$step_size))
  }
  if(approximation =="MLE"){
    sn_approx <- fit_sn_MLE (samples)
  }
  if(approximation == "Moment"){
    X <- samples
    mu <- apply(X,2,mean)
    variance <- apply(X,2,var)
    skew <- apply(X,2,skewness)
    sn_approx <- fit_sn_moment (mu,variance,skew)
  }
  xi <- sn_approx$ xi
  omega <- sn_approx$ omega
  alpha <- sn_approx$ alpha
  if(map == "direcct"){
    result <- transformation_componentwise_direct_map_sn (xi,omega,alpha,log_pi,nits,target_a)
  }
  if(map == "indirect"){
    result <- transformation_map_sn (xi,omega,alpha,log_pi,nits,target_a)
  }
  return(list(xi = xi,omega = omega,alpha = alpha,samples_x = result$samples_x,
              samples_y = result$samples_y, a_rate = result$a_rate, step_size = result$step_size))
}

# Sampling with online MLE
Transformation_sn_MLE <- function(log_pi,nits,x_curr,Sampling_algorithm = Adaptive_RWM, d_logpi=NULL,
                                            y_curr = rnorm(length(x_curr)),chain_len = 5,h = 0.9,
                                            t0 = 2000,delta = 0.9,c = 200, t1 = 3000,
                                            target_a = 0.44){
  
  # get the dimension
  d <- length(x_curr)
  # Get the pre samples 
  pre_sampling <- Sampling_algorithm(logpi = log_pi,nits = t0, h = h,
                                     x_curr = x_curr,dlogpi = d_logpi )
  pre_samples <- as.matrix(pre_sampling$x_store)
  #effectiveSize(pre_samples)
  
  # get the parameters of approximation distribution
  para <- fit_sn_MLE (pre_samples)
  xi_est <- para$xi; omega_est <- para$omega;alpha_est <- para$alpha
  #print(para)
  # create matrix to store samples
  y_store <- matrix(nrow = nits, ncol = d)
  x_store <- matrix(nrow = nits, ncol = d)
  v_store <- sort(rnorm(nits))
  
  # define the z_curr
  y_curr <- rep(0,d)
  accepted <- 0
  
  # define the identical matrix
  I_d <- diag(rep(1,d))
  # initialise V and mu
  V <- I_d
  mu <- y_curr
  
  for( i in 1:nits){
    
    delta_est <- alpha_est / sqrt(1+alpha_est^2)
    # Get A,B,C
    A <- omega_est * sqrt(1 - delta_est^2)
    B <- omega_est * delta_est
    C <- xi_est
    
    # Define the transformed density p_Y(z)
    log_piy <- function(y,v){
      x <- A * y + B * rep(abs(v),d) + C
      logf <- log_pi(x) - v^2/2
      logf <- max(logf,-1e30)
      return(logf)
    }
    
    v <- v_store[i]
    
    if(i == 1){
      logpi_curr <- log_piy(y_curr, v)
    }
    
    # sample chain_len y from the y|v_i
    for (j in 1:chain_len) {
      y_prop <- y_curr + h * as.vector(rmvnorm(1, sigma = V))
      logpi_prop <- log_piy(y_prop, v)
      
      if (!is.finite(logpi_prop)) next  # skip the abnormal value
      
      loga <- logpi_prop - logpi_curr
      
      if (log(runif(1)) < loga) {
        y_curr <- y_prop
        logpi_curr <- logpi_prop
        accepted <- accepted + 1
      }
    }
    
    # store y and x
    y_store[i, ] <- y_curr
    x_curr <- A * y_curr + B * rep(abs(v),d) + C
    x_store[i,] <- x_curr
    
    # set the learning rate
    gamma <- min(0.7,c* 1/(i)^delta) # avoid the V explosion
    # update the covariance matrix
    V <- V + gamma * ((y_curr - mu) %*% t(y_curr - mu) - V)
    V <- V + 1e-7 * I_d 
    # update the mean
    mu <- mu + gamma * (y_curr - mu)
    # update the step size
    alpha <- ifelse(!is.nan(loga),min(1, exp(loga)),target_a)
    h <- h * sqrt(exp(gamma * (alpha - target_a)))
    
    #print(sample_cov)
    
    if(i %% t1 == 0){
      para <- sgd_mle(rbind(as.matrix(x_store[1:i,]),pre_samples),
                      initial_params = list(xi = xi_est,omega = omega_est,alpha = alpha_est),
                      gradient = sn_log_gradient,
                      loglik_single = log_dsn,
                      force_change = function(params){params$omega = pmax(params$omega,0.01);
                      params$xi = pmin(pmax(params$xi,-1e5),1e5);
                      params$alpha = pmin(pmax(params$alpha,-1e10),1e10);
                      return(params)}
      )$parameters
      xi_est <- para$xi; omega_est <- para$omega;alpha_est <- para$alpha
      #print(as.data.frame(para))
    }
    if (floor(i/10000)==i/10000) { cat(i," iterations")}
  }
  return(list(xi = xi_est,omega = omega_est,alpha = alpha_est,
              samples_x = x_store,samples_y = y_store))
}

# Sampling with online Moment Methods
Transformation_sn_mom <- function(log_pi,nits,x_curr,Sampling_algorithm = Adaptive_RWM, d_logpi=NULL,
                                         y_curr = rnorm(length(x_curr)),chain_len = 5,h = 0.9,
                                         t0 = 2000,delta = 0.7,c = 1, t1 = 100,
                                         target_a = 0.44){
  
  # get the dimension
  d <- length(x_curr)
  # Get the pre samples 
  pre_sampling <- Sampling_algorithm(logpi = log_pi,nits = t0, h = h,
                                     x_curr = x_curr,dlogpi = d_logpi )
  pre_samples <- as.matrix(pre_sampling$x_store)
  #effectiveSize(pre_samples)
  
  # construct matrix to store moments
  mu_store <- matrix(nrow = nits, ncol = d)
  var_store <- matrix(nrow = nits, ncol = d)
  skew_store <- matrix(nrow = nits, ncol = d)
  kurt_store <- matrix(nrow = nits, ncol = d)
  
  
  # calculate the sample moments
  sample_mean <- apply(pre_samples,2,mean)
  sample_cov <- as.matrix(cov(pre_samples))
  sample_skew <- apply(pre_samples,2,skewness)
  sample_kurt <- apply(pre_samples,2,kurtosis)
  
  # get the parameters of approximation distribution
  para <- fit_sn_moment(sample_mean,diag(sample_cov),sample_skew)
  xi_est <- para$xi; omega_est = para$omega ; alpha_est <- para$alpha
  #print(as.data.frame( para))
  # create matrix to store samples
  y_store <- matrix(nrow = nits, ncol = d)
  x_store <- matrix(nrow = nits, ncol = d)
  v_store <- sort(rnorm(nits))
  
  # define the z_curr
  y_curr <- rep(0,d)
  accepted <- 0
  
  # define the identical matrix
  I_d <- diag(rep(1,d))
  # initialise V and mu
  V <- I_d
  mu <- y_curr
  
  for( i in 1:nits){
    
    delta_est <- alpha_est / sqrt(1+alpha_est^2)
    # Get A,B,C
    A <- omega_est * sqrt(1 - delta_est^2)
    B <- omega_est * delta_est
    C <- xi_est
    
    # Define the transformed density p_Y(z)
    log_piy <- function(y,v){
      x <- A * y + B * rep(abs(v),d) + C
      logf <- log_pi(x) - v^2/2
      logf <- max(logf,-1e30)
      return(logf)
    }
    
    v <- v_store[i]
    
    if(i == 1){
      logpi_curr <- log_piy(y_curr, v)
    }
    
    # sample chain_len y from the y|v_i
    for (j in 1:chain_len) {
      y_prop <- y_curr + h * as.vector(rmvnorm(1, sigma = V))
      logpi_prop <- log_piy(y_prop, v)
      
      if (!is.finite(logpi_prop)) next  # skip the abnormal value
      
      loga <- logpi_prop - logpi_curr
      
      if (log(runif(1)) < loga) {
        y_curr <- y_prop
        logpi_curr <- logpi_prop
        accepted <- accepted + 1
      }
    }
    
    # store y and x
    y_store[i, ] <- y_curr
    x_curr <- A * y_curr + B * rep(abs(v),d) + C
    x_store[i,] <- x_curr
    
    # set the learning rate
    gamma <- min(0.7,c* 1/(i)^delta) # avoid the V explosion
    # update the covariance matrix
    V <- V + gamma * ((y_curr - mu) %*% t(y_curr - mu) - V)
    V <- V + 1e-7 * I_d 
    
    # update the mean
    mu <- mu + gamma * (y_curr - mu)
    # update the step size
    alpha <- ifelse(!is.nan(loga),min(1, exp(loga)),target_a)
    h <- h * sqrt(exp(gamma * (alpha - target_a)))
    
    #update the moments
    mom <- update_moments(x_curr,t0 + i -1,sample_mean,sample_cov,sample_skew)
    sample_mean <- mom$mean
    sample_cov <- as.matrix(mom$variance)
    sample_skew <- mom$skew
    sample_kurt <- mom$kurt
    mu_store[i,] <- sample_mean
    var_store[i,] <- diag(sample_cov)
    skew_store[i,] <- sample_skew
    kurt_store[i,] <- sample_kurt
    
    #print(sample_cov)
    
    if(i%% t1 == 0){
      para <- fit_sn_moment(sample_mean,diag(sample_cov),sample_skew)
      xi_est <- para$xi; omgea_est = para$omega ; alpha_est <- para$alpha
      # print(as.data.frame( para))
      # cat("the y is ",y_curr, " the x is ",x_curr)

    }
    if (floor(i/10000)==i/10000) { cat(i," iterations")}
  }
  
  return(list(xi = xi_est,omega = omega_est,alpha = alpha_est,
              samples_x = x_store,samples_y = y_store, mu_store = mu_store,
              var_store = var_store,skew_store = skew_store,kurt_store = kurt_store))
}

