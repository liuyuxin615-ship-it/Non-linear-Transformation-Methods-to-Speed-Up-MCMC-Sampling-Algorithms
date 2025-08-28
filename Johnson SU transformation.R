library(moments)
# define a funciton to generate samples from johnsonSU distribution
generate_john <- function(nits,xi,lambda,gamma,delta){
  return(lambda * sinh((rnorm(nits) - gamma) / delta) + xi)
}

# Define the fitting function in one dimension using MLE
fit_oned_johnsonsu_MLE <- function(samples){
  # Define a negative log-likelihood function of jognson su distribution 
  loglik_johnson_su <- function(par, x) {
    gamma <- par[1]
    delta <- par[2]
    xi    <- par[3]
    lambda<- par[4]
    
    pi <- 3.1415926535
    if (delta <= 0 || lambda <= 0) return(Inf)
    y <- (x - xi) / lambda
    z <- gamma + delta * asinh(y)
    
    log_pdf <- log(delta) - log(lambda) - 0.5 * log(2*pi) - 0.5 * z^2 + log(1 / sqrt(1 + y^2))
    log_pdf <- pmax(log_pdf, -1e30)
    return(-sum(log_pdf))  # pay attention it is negative 
  }
  start <- c(0, 1, 0, 1)
  
  # calculate the maximum likelihood estimator
  result_ml <- optim(par=start, fn=loglik_johnson_su, x=samples, method="L-BFGS-B", 
                     lower=c(-Inf, 1e-6, -Inf, 1e-6))
  # extrat it from the result
  gamma <- result_ml$par [1]
  delta <- result_ml$par [2]
  xi    <- result_ml$par [3]
  lambda<- result_ml$par [4]
  
  # return the parameters
  c(xi, lambda, gamma, delta)
}

# Find the approximation in high dimension using component-wise method
fit_johnsonsu_MLE <- function(X) {
  if(!is.matrix(X)){
    return("the input must be a matrix")
  }
  # Find the Johnson approximation in each dimension
  r <- as.data.frame(t(apply(X,2,fit_oned_johnsonsu_MLE)))
  colnames(r) <- c("xi","lambda","gamma","delta")
  return(r)
}

# for one dimension, the input is the first four moment of samples
fit_oned_johnsonsu_moment <- function(moment_vec) {
  
  sample_mean = moment_vec[1];sample_var = moment_vec[2]
  sample_skew = moment_vec[3];sample_kurt = moment_vec[4]
  #print(moment_vec)
  # 初始值（log 参数确保lambda和delta为正）
  xi_init <- sample_mean
  log_lambda_init <- log(sqrt(sample_var))
  gamma_init <- 1
  log_delta_init <- log(1)
  
  start <- c(xi_init, log_lambda_init, gamma_init, log_delta_init)
  
  moment_equations <- function(params) {
    
    xi <- params[1]
    lambda <- exp(params[2])  # ensure > 0
    gamma <- params[3]
    delta <- exp(params[4])   # ensure > 0
    clip <- function(x, bound=1e20) {
      pmax(pmin(x, bound), -bound)
    }
    # 理论矩
    mean_theory <- xi - lambda * exp(delta^-2 / 2) * sinh(gamma / delta)
    var_theory <- (lambda^2 / 2) * (exp(delta^-2) - 1) *
      (exp(delta^-2) * cosh(2 * gamma / delta) + 1)
    
    var_theory <- max(var_theory, 1e-12)
    
    skew_theory <- -(lambda^3 * exp(-delta^-2/2) * (exp(delta^-2) - 1)^2 *
                       ((exp(delta^-2) * (exp(delta^-2) + 2)) * sinh(3 * gamma / delta) +
                          3 * sinh(gamma / delta))) / 4*(var_theory)^1.5
    
    term <- exp(pmin(delta^-2, 50))
    K1 <- term^2 * (term^4 + 2 * term^3 + 3 * term^2 - 3) * cosh(4 * gamma / delta)
    K2 <- 4 * term^2 * (term + 2) * cosh(3 * gamma / delta)
    K3 <- 3 * (2 * term + 1)

    
    kurt_theory <- clip((lambda^4 * (term - 1)^2 * (K1 + K2 + K3)) / pmin(8 * (var_theory)^2,1e300) + 3)
   
    
    
    c(
      clip(mean_theory - sample_mean),
      clip(var_theory - sample_var),
      clip(skew_theory - sample_skew),
      clip(kurt_theory - sample_kurt)
    )
  }
  
  res <- nleqslv(start, moment_equations, method = "Broyden",
                 control = list(maxit = 500, ftol = 1e-8))
  
  xi <- res$x[1]
  lambda <- exp(res$x[2])
  gamma <- res$x[3]
  delta <- exp(res$x[4])
  
  c(xi, lambda, gamma, delta)
}





# required the X is the matrix of moment in each dimension, each row represent one dimension
fit_johnsonsu_moment <- function(X) {
  if(!is.matrix(X)){
    return("the input must be a matrix")
  }
  r <- as.data.frame(t(apply(X,1,fit_oned_johnsonsu_moment)))
  colnames(r) <- c("xi","lambda","gamma","delta")
  return(r)
}


# Define the transformation map
transformation_map_johnsonsu <- function(xi,lambda,gamma,delta,nits,
                                                 log_pi,x_curr = rep(0,length(xi)), target_a= 0.44){
  d <- length(xi)
  # Define inverse map: x = F⁻¹(Φ(y))
  
  x_from_y <- function(y) sinh((y - gamma)/delta)*lambda + xi 
  
  # Define the transformed density p_Y(y)
  log_p_Y <- function(y) {
    x <- x_from_y(y)
    logf <- log_pi(x) + sum(log(1 + ((x - xi) / lambda)^2 ) /2)
    logf <- max(logf,-1e30)
    if(is.nan(logf)) logf <- -1e30
    return(logf)
  }
  
  # Sampling from log_p_Y
  chain<- Adaptive_RWM(log_p_Y,nits = nits,h = 0.5,x_curr = x_curr,target_a = target_a)
  
  # Get the samples of y
  samples_y <- chain$x_store
  
  # transform y back
  samples_x<- t(apply(samples_y,1,x_from_y)) 
  if(dim(samples_x)[1]<dim(samples_x)[2]) samples_x = t(samples_x)
  
  return(list(samples_x = samples_x,samples_y = samples_y,
              a_rate = chain$a_rate,step_size = chain$step_size))
}
# Define the final sampling algorithm
Sampling_trsanformation_johnsonsu <- function(samples,log_pi,nits,method = "MLE",
                                              x_curr = rep(0,ncol(samples)), target_a = 0.23){
  
  if(method == "MLE") johnsonsu_approx <- fit_johnsonsu_MLE(samples)
  if(method == "Moment") {
    sample_mean <- apply(samples,2,mean)
    sample_var <- apply(samples,2,var)
    sample_skew <- apply(samples,2,skewness)
    sample_kurt <- apply(samples,2,kurtosis)
    # get the parameters of approximation distribution
    moment_matrix <- cbind(sample_mean,sample_var,sample_skew,sample_kurt)
    johnsonsu_approx <- fit_johnsonsu_moment(moment_matrix)
    }
  xi <- johnsonsu_approx$xi
  lambda <- johnsonsu_approx$lambda
  gamma <- johnsonsu_approx$gamma
  delta <- johnsonsu_approx$delta
  result <- transformation_map_johnsonsu(xi,lambda,gamma,delta,nits,log_pi,x_curr,target_a)
  return(list(xi = xi,lambda = lambda,gamma = gamma,delta = delta,
              samples_x = result$samples_x,
              samples_y = result$samples_y,
              a_rate = result$a_rate,step_size = result$step_size))
}

# The online MLE for sampling
Transformation_johnsonsu_MLE <- function(log_pi,nits,x_curr,Sampling_algorithm = Adaptive_RWM, d_logpi=NULL,
                                         y_curr = rnorm(length(x_curr)),h = 0.9,
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
  para <- fit_johnsonsu_MLE(pre_samples)
  xi_est <- para$xi; lambda_est <- para$lambda ; gamma_est <- para$gamma; delta_est = para$delta
  
  # create matrix to store samples
  y_store <- matrix(nrow = nits, ncol = d)
  x_store <- matrix(nrow = nits, ncol = d)
  
  # define the y_curr
  y_curr <- rep(0,d)
  accepted <- 0
  
  # define the identical matrix
  I_d <- diag(rep(1,d))
  # initialise V and mu
  V <- I_d
  mu <- y_curr
  
  for( i in 1:nits){
    
    x_from_y <- function(y) sinh((y - gamma_est)/delta_est)*lambda_est + xi_est 
    
    # Define the transformed density p_Y(y)
    log_piy <- function(y) {
      x <- x_from_y(y)
      logf <- log_pi(x) + sum(log(1 + ((x - xi_est) / lambda_est)^2 ) /2)
      logf <- max(logf,-1e30)
      if(is.nan(logf)) logf <- -1e30
      return(logf)}
    
    if(i == 1){
      logpi_curr <- log_piy(y_curr)
    }
    
    # generate y_propose
    y_prop <- y_curr + h * as.vector(rmvnorm(1, sigma =  V))
    logpi_prop <- log_piy(y_prop) # calculate the density
    loga <- logpi_prop - logpi_curr
    
    u <- runif(1)
    if ( is.finite(loga) && log(u) < loga) {
      y_curr <- y_prop
      logpi_curr <- logpi_prop
      accepted <- accepted + 1
    }
    
    # store y and x
    y_store[i, ] <- y_curr
    x_curr <- x_from_y(y_curr)
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
      para <- sgd_mle(rbind(pre_samples,as.matrix(x_store[1:i,])),
                      initial_params = list (lambda = lambda_est,delta = delta_est,
                                             xi = xi_est,gamma = gamma_est))$parameters
      xi_est <- para$xi; lambda_est <- para$lambda ; gamma_est <- para$gamma; delta_est = para$delta
    }
    if (floor(i/10000)==i/10000) { cat(i," iterations")}
  }
  return(list(xi = xi_est,lambda = lambda_est,gamma = gamma_est,delta = delta_est,
              samples_x = x_store,samples_y = y_store))
}

# The online Mom for sampling
Transformation_johnsonsu_mom <- function(log_pi,nits,x_curr,Sampling_algorithm = Adaptive_RWM, d_logpi=NULL,
                                         y_curr = rnorm(length(x_curr)),h = 0.9,
                                         t0 = 2000,delta = 0.7,c = 1, t1 = 5000,
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
  moment_matrix <- cbind(sample_mean,diag(sample_cov),sample_skew,sample_kurt)
  para <- fit_johnsonsu_moment(moment_matrix)
  xi_est <- para$xi; lambda_est <- para$lambda ; gamma_est <- para$gamma; delta_est = para$delta
  
  # create matrix to store samples
  y_store <- matrix(nrow = nits, ncol = d)
  x_store <- matrix(nrow = nits, ncol = d)
  
  # define the y_curr
  y_curr <- rep(0,d)
  accepted <- 0
  
  # define the identical matrix
  I_d <- diag(rep(1,d))
  # initialise V and mu
  V <- I_d
  mu <- y_curr
  
  for( i in 1:nits){
    
    x_from_y <- function(y) sinh((y - gamma_est)/delta_est)*lambda_est + xi_est 
    
    # Define the transformed density p_Y(y)
    log_piy <- function(y) {
      x <- x_from_y(y)
      logf <- log_pi(x) + sum(log(1 + ((x - xi_est) / lambda_est)^2 ) /2)
      logf <- max(logf,-1e30)
      if(is.nan(logf)) logf <- -1e30
      return(logf)}
    
    if(i == 1){
      logpi_curr <- log_piy(y_curr)
    }
    
    # generate y_propose
    y_prop <- y_curr + h * as.vector(rmvnorm(1, sigma =  V))
    logpi_prop <- log_piy(y_prop) # calculate the density
    loga <- logpi_prop - logpi_curr
    
    u <- runif(1)
    if ( is.finite(loga) && log(u) < loga) {
      y_curr <- y_prop
      logpi_curr <- logpi_prop
      accepted <- accepted + 1
    }
    
    # store y and x
    y_store[i, ] <- y_curr
    x_curr <- x_from_y(y_curr)
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
      # get the new approximation
      moment_matrix <- cbind(sample_mean,diag(sample_cov),sample_skew,sample_kurt)
      #print(moment_matrix)
      para <- fit_johnsonsu_moment(moment_matrix)
      xi_est <- para$xi; lambda_est <- para$lambda ; gamma_est <- para$gamma; delta_est = para$delta
    }
    if (floor(i/10000)==i/10000) { cat(i," iterations")}
  }
  
  return(list(xi = xi_est,lambda = lambda_est,gamma = gamma_est,delta = delta_est,
              samples_x = x_store,samples_y = y_store, mu_store = mu_store,
              var_store = var_store,skew_store = skew_store,kurt_store = kurt_store))
}


