#define the function to genertae sample form sinh-arcsinh
generate_sinh <- function(nits,xi,epsilon,eta,delta){
  return(pmin(xi + eta * sinh((asinh(rnorm(nits)) +epsilon)/delta),1e300))
}

# Sinh-arcsinh moment method to estimation
modified_bessel_K <- function(nu, z) {
  return(besselK(z, nu))
}

# 计算P_q值的函数
compute_Pq <- function(q) {
  # 根据图片中的最终公式
  # P_q = (e^(1/4) / (8π)^(1/2)) * {K_{(q+1)/2}(1/4) + K_{(q-1)/2}(1/4)}
  
  coeff <- exp(1/4) / sqrt(8 * pi)
  
  # 计算修正贝塞尔函数
  K_term1 <- modified_bessel_K((q + 1)/2, 1/4)
  K_term2 <- modified_bessel_K((q - 1)/2, 1/4)
  
  return(coeff * (K_term1 + K_term2))
}

# using moment method to estimate the parameter
Sinh_para_est_mom <- function(x_bar, m2, m3, m4) {
  # 输入：
  # mu1: 一阶中心矩（均值，通常为0）
  # mu2: 二阶中心矩（方差）
  # mu3: 三阶中心矩
  # mu4: 四阶中心矩
  
  # 预计算P值（这些在实际应用中可能需要数值积分）
  # 根据公式，我们需要P_{1/δ}, P_{2/δ}, P_{3/δ}, P_{4/δ}
  
  # 定义目标函数来求解参数
  objective_function <- function(params) {
    epsilon <- params[1]
    delta <- params[2]
    xi <- params[3]
    eta <- params[4]
    
    if (delta <= 0 ||eta<=0 ) return(1e30)  # delta必须为正
    
    # 计算P值
    P_1_delta <- compute_Pq(1/delta)
    P_2_delta <- compute_Pq(2/delta)
    P_3_delta <- compute_Pq(3/delta)
    P_4_delta <- compute_Pq(4/delta)
    
    # 根据图片中的矩公式计算理论矩
    # E(X_{ε,δ}) = sinh(ε/δ) * P_{1/δ}
    mu1 <- sinh(epsilon/delta) * P_1_delta
    
    # E(X²_{ε,δ}) = (1/2){cosh(2ε/δ) * P_{2/δ} - 1}
    mu2 <- 0.5 * (cosh(2*epsilon/delta) * P_2_delta - 1)
    
    # E(X³_{ε,δ}) = (1/4){sinh(3ε/δ) * P_{3/δ} - 3*sinh(ε/δ) * P_{1/δ}}
    mu3 <- 0.25 * (sinh(3*epsilon/delta) * P_3_delta - 
                     3*sinh(epsilon/delta) * P_1_delta)
    
    # E(X⁴_{ε,δ}) = (1/8){cosh(4ε/δ) * P_{4/δ} - 4*cosh(2ε/δ) * P_{2/δ} + 3}
    mu4 <- 0.125 * (cosh(4*epsilon/delta) * P_4_delta - 
                      4*cosh(2*epsilon/delta) * P_2_delta + 3)
    
    # compute the central moments and after scale
    theoretical_x_bar <- eta * mu1 + xi
    theoretical_m2 <- (mu2 - mu1^2) * eta^2
    theoretical_m3 <- (mu3 - 3*mu1*mu2 + 2*mu1^3) *eta^3
    theoretical_m4 <- (mu4 - 4*mu1*mu3 + 6*mu1^2*mu2 - 3*mu1^4) * eta^4
    
    # 计算误差平方和
    error <- (theoretical_x_bar - x_bar)^2 + 
      (theoretical_m2 - m2)^2 + 
      (theoretical_m3 - m3)^2 + 
      (theoretical_m4 - m4)^2
    
    return(pmin(error,1e200))
  }

  # 使用优化算法求解
  # 初始猜测值
  initial_guess <- c(0, 1,1,1)  # epsilon = 0, delta = 1
  
  # 使用Nelder-Mead方法优化
  result <- optim(initial_guess, objective_function, 
                  method = "Nelder-Mead",
                  control = list(maxit = 10000))
  
  if (result$convergence != 0) {
    #warning("优化可能未收敛")
  }
  
  epsilon_hat <- result$par[1]
  delta_hat <- result$par[2]
  xi_hat <- result$par[3]
  eta_hat <- result$par[4]
  
  # 计算拟合优度
  fitted_error <- result$value
  
  return(list(
    epsilon = epsilon_hat,
    delta = delta_hat,
    xi = xi_hat,
    eta = eta_hat,
    convergence = result$convergence,
    fitted_error = fitted_error,
    optimization_result = result
  ))
}

# using MLE to estimate the parameter
Sinh_para_est_MLE <- function(samples){
  # Define a negative log-likelihood function of sinh_arcsinh distribution 
  loglik_sinh_arcsinh<- function(par, x) {
    epsilon <- par[1];delta <- par[2];xi<- par[3];eta <- par[4]
    pi <- base::pi
    y <- (x - xi) / eta
    z <- pmin(pmax(sinh(delta * asinh(y) - epsilon),-1e20),1e20)
    term1 <- log(delta) - log(eta)
    term2 <- 0.5 * (log(1 + z^2) - log(2*pi*(1+y^2)))
    term3 <- -0.5 * z^2
    log_pdf <- sum(term1 + term2 + term3)
    log_pdf <- pmin(pmax(log_pdf, -1e300),1e300)
    return(-(log_pdf))  # pay attention it is negative 
  }
  start <- c(1, 1, mean(samples), sd(samples))
  
  # calculate the maximum likelihood estimator
  result_ml <- optim(par=start, fn=loglik_sinh_arcsinh, x=samples, method="L-BFGS-B", 
                     lower=c(-Inf, 1e-6, -Inf, 1e-6))
  # extrat it from the result
  epsilon<- result_ml$par [1]
  delta <- result_ml$par [2]
  xi    <- result_ml$par [3]
  eta<- result_ml$par [4]
  
  # return the parameters
  c(xi, epsilon, eta, delta)
}


# Extend the parameter estimation into high dimension
fit_highdim_sinh_arcsinh <- function(X,method = "MLE"){
  
  if(!is.matrix(X)){
    return("the input must be a matrix")
  }
  # Find the Johnson approximation in each dimension
  if(method == "MLE"){
    r <- as.data.frame(t(apply(X,2,Sinh_para_est_MLE)))
    colnames(r) <- c("xi","epsilon","eta","delta")
  }
  if(method == "Moment"){
    x_bar <- apply(X,2,mean)
    m2 <- pmin(apply(X,2,FUN = function(x) mean((x-mean(x))^2)),1e200)
    m3 <- apply(X,2,FUN = function(x) mean((x-mean(x))^3))
    m4 <- pmin(apply(X,2,FUN = function(x) mean((x-mean(x))^4)),1e200)
    d <- dim(X)[2]
    r <- as.data.frame(matrix(0,ncol = 4,nrow = d))
    colnames(r) <- c("xi","epsilon","eta","delta")
    for (i in 1:d){
      result <- Sinh_para_est_mom(x_bar[i],m2[i],m3[i],m4[i])
      r[i,"xi"] <- result$xi
      r[i,"epsilon"] <- result$epsilon
      r[i,"eta"] <- result$eta
      r[i,"delta"] <- result$delta
    }
  }
  return(r)
}

# Define the transformation map
transformation_map_sinh_arcsinh <- function(xi,epsilon,eta,delta,nits,log_pi,
                                                    x_curr = rep(0,length(xi)), target_a= 0.23){
  # Get the dimension
  d <- length(xi)
  
  # Define inverse map: x = F⁻¹(Φ(z))
  x_from_z <- function(z) xi + eta * sinh((asinh(z)+epsilon)/delta)
  y_from_x <- function(x) (x - xi) / eta
  
  # Define the transformed density p_Y(z)
  log_p_z <- function(z) {
    x <- x_from_z(z)
    y <- y_from_x(x)
    
    logf <- log_pi(x) - sum(log(1 + z^2) / 2 + log(1 + y^2) / 2)
    logf <- max(logf,-1e30)
    if(is.nan(logf)) logf <- -1e30
    return(logf)
  }
  
  # Sampling from log_o_Y
  chain<- Adaptive_RWM(log_p_z,nits = nits,x_curr = x_curr,target_a = target_a,h = 0.5)
  
  # Get the samples of y
  samples_z <- chain$x_store
  
  # transform y back
  samples_x<- t(apply(samples_z,1,x_from_z)) 
  if(dim(samples_x)[1]<dim(samples_x)[2]) samples_x = t(samples_x)
  
  return(list(samples_x = samples_x,samples_y = samples_z,a_rate = chain$a_rate,
              step_size = chain$step_size))
}
# define the sampling algorithm
Sampling_trsanformation_sinh_arcsinh<- function(samples,log_pi,nits,method = "Moment",
                                                x_curr = rep(0,ncol(samples)),target_a = 0.23){
  sinh_arcsin_approx <- fit_highdim_sinh_arcsinh(samples,method)
  xi <- sinh_arcsin_approx$xi
  eta <- sinh_arcsin_approx$eta
  epsilon <- sinh_arcsin_approx$epsilon
  delta <- sinh_arcsin_approx$delta
  result <- transformation_map_sinh_arcsinh(xi,epsilon,eta,delta,nits,log_pi,x_curr ,target_a)
  return(list(xi = xi,eta = eta,epsilon = epsilon,delta = delta,
              samples_x = result$samples_x,
              samples_y = result$samples_y,
              a_rate = result$a_rate, step_size = result$step_size))
}

# Define the online MLE sampling
Transformation_sinh_arcsinh_MLE <- function(log_pi,nits,x_curr,Sampling_algorithm = Adaptive_RWM, d_logpi=NULL,
                                         y_curr = rnorm(length(x_curr)),h = 0.9,
                                         t0 = 2000,delta = 0.9,c = 200, t1 = 3000,
                                         target_a = 0.44){
 
  # get the dimension
  d <- length(x_curr)
  # Get the pre samples 
  pre_sampling <- Sampling_algorithm(logpi = log_pi,nits = t0, h = h,
                                     x_curr = x_curr,dlogpi = d_logpi )
  pre_samples <- as.matrix(pre_sampling$x_store)
  # effectiveSize(pre_samples)
  
  # get the parameters of approximation distribution
  para <- fit_highdim_sinh_arcsinh(pre_samples,"Moment")
  xi_est <- para$xi; epsilon_est <- para$epsilon ; eta_est <- para$eta; delta_est = para$delta
  print(para)
  # create matrix to store samples
  z_store <- matrix(nrow = nits, ncol = d)
  x_store <- matrix(nrow = nits, ncol = d)
  
  # define the z_curr
  z_curr <- rep(0,d)
  accepted <- 0
  
  # define the identical matrix
  I_d <- diag(rep(1,d))
  # initialise V and mu
  V <- I_d
  mu <- z_curr
  
  for( i in 1:nits){
    
    # Define inverse map: x = F⁻¹(Φ(z))
    x_from_z <- function(z) xi_est + eta_est * sinh((asinh(z)+epsilon_est)/delta_est)
    y_from_x <- function(x) (x - xi_est) / eta_est
    
    # Define the transformed density p_Y(z)
    log_p_z <- function(z) {
      x <- x_from_z(z)
      y <- y_from_x(x)
      
      logf <- log_pi(x) - sum(log(1 + z^2) / 2 + log(1 + y^2) / 2)
      logf <- max(logf,-1e30)
      if(is.nan(logf)) logf <- -1e30
      return(logf)
    }
    
    if(i == 1){
      logpi_curr <- log_p_z(z_curr)
    }
    
    # generate z_propose
    z_prop <- z_curr + h * as.vector(rmvnorm(1, sigma =  V))
    logpi_prop <- log_p_z(z_prop) # calculate the density
    loga <- logpi_prop - logpi_curr
    
    u <- runif(1)
    if ( is.finite(loga) && log(u) < loga) {
      z_curr <- z_prop
      logpi_curr <- logpi_prop
      accepted <- accepted + 1
    }
    
    # store y and x
    z_store[i, ] <- z_curr
    x_curr <- x_from_z(z_curr)
    x_store[i,] <- x_curr
    
    # set the learning rate
    gamma <- min(0.7,c* 1/(i)^delta) # avoid the V explosion
    # update the covariance matrix
    V <- V + gamma * ((z_curr - mu) %*% t(z_curr - mu) - V)
    V <- V + 1e-7 * I_d 
    # update the mean
    mu <- mu + gamma * (z_curr - mu)
    # update the step size
    alpha <- ifelse(!is.nan(loga),min(1, exp(loga)),target_a)
    h <- h * sqrt(exp(gamma * (alpha - target_a)))
    
    #print(sample_cov)
    
    if(i %% t1 == 0){
      para <- sgd_mle(rbind(pre_samples,as.matrix(x_store[1:i,])),
                      initial_params = list(xi = xi_est,epsilon = epsilon_est,eta = eta_est,delta = delta_est),
                      gradient = sinh_log_gradient,
                      loglik_single = log_dsinh_arcsinh,
                      force_change = function(params){params$eta = pmin(pmax(params$eta,0.001),1e15);
                      params$delta = pmin(pmax(params$delta,0.001),1e15);
                      params$epsilon = pmin(pmax(params$epsilon,-1e15),1e15);
                      params$xi = pmin(pmax(params$xi,-1e15),1e15)
                      return(params)})$parameters
      # grad <- s$grad
      # print("the gradient is:")
      # print(grad)
      # print("the parameter is")
      # print(para)
      xi_est <- para$xi; epsilon_est <- para$epsilon ; eta_est <- para$eta; delta_est = para$delta
    }
    
    if (floor(i/10000)==i/10000) { cat(i," iterations")}
  }
  
  return(list(xi = xi_est,epsilon = epsilon_est,eta= eta_est,delta = delta_est,
              samples_x = x_store,samples_y = z_store))
}

# define a function to get the estimation from moment
fit_sinh_arcsinh_moment <- function(moment_matrix){
  #moment_matrix is a matrix that the first column is mean, second is variance, third is skew, fourth is kurtosis.
  central_mom <- moment_matrix
  central_mom[,3] <- central_mom[,3]*central_mom[,2]^1.5
  central_mom[,4] <- central_mom[,4]*central_mom[,2]^2
  d <- length(moment_matrix[,1])
  r <- as.data.frame(matrix(0,ncol = 4,nrow = d))
  colnames(r) <- c("xi","epsilon","eta","delta")
  for(i in 1:d){
    result <- Sinh_para_est_mom(central_mom[i,1],central_mom[i,2],central_mom[i,3],central_mom[i,4])
    r[i,"xi"] <- result$xi
    r[i,"epsilon"] <- result$epsilon
    r[i,"eta"] <- result$eta
    r[i,"delta"] <- result$delta
  }
  return(r)
}

# The online Mom for sampling
Transformation_sinh_arcsinh_mom <- function(log_pi,nits,x_curr,Sampling_algorithm = Adaptive_RWM, d_logpi=NULL,
                                         y_curr = rnorm(length(x_curr)),h = 0.9,
                                         t0 = 2000,delta = 0.7,c = 1, t1 = 500,
                                         target_a = 0.44){
  
  # get the dimension
  d <- length(x_curr)
  # Get the pre samples 
  pre_sampling <- Sampling_algorithm(logpi = log_pi,nits = t0, h = h,
                                     x_curr = x_curr,dlogpi = d_logpi )
  pre_samples <- as.matrix(pre_sampling$x_store)
  effectiveSize(pre_samples)
  
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
  para <- fit_sinh_arcsinh_moment(moment_matrix)
  xi_est <- para$xi; epsilon_est <- para$epsilon ; eta_est <- para$eta; delta_est = para$delta
  
  # construct matrix to store the samples
  z_store <- matrix(nrow = nits, ncol = d)
  x_store <- matrix(nrow = nits, ncol = d)
  
  # define the z_curr
  z_curr <- rep(0,d)
  accepted <- 0
  
  # define the identical matrix
  I_d <- diag(rep(1,d))
  # initialise V and mu
  V <- I_d
  mu <- z_curr
  
  for( i in 1:nits){
    
    # Define inverse map: x = F⁻¹(Φ(z))
    x_from_z <- function(z) xi_est + eta_est * sinh((asinh(z)+epsilon_est)/delta_est)
    y_from_x <- function(x) (x - xi_est) / eta_est
    
    # Define the transformed density p_Y(z)
    log_p_z <- function(z) {
      x <- x_from_z(z)
      y <- y_from_x(x)
      
      logf <- log_pi(x) - sum(log(1 + z^2) / 2 + log(1 + y^2) / 2)
      logf <- max(logf,-1e30)
      if(is.nan(logf)) logf <- -1e30
      return(logf)
    }
    
    if(i == 1){
      logpi_curr <- log_p_z(z_curr)
    }
    
    # generate z_propose
    z_prop <- z_curr + h * as.vector(rmvnorm(1, sigma =  V))
    logpi_prop <- log_p_z(z_prop) # calculate the density
    loga <- logpi_prop - logpi_curr
    
    u <- runif(1)
    if ( is.finite(loga) && log(u) < loga) {
      z_curr <- z_prop
      logpi_curr <- logpi_prop
      accepted <- accepted + 1
    }
    
    # store y and x
    z_store[i, ] <- z_curr
    x_curr <- x_from_z(z_curr)
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
      para <- fit_sinh_arcsinh_moment(moment_matrix)
      xi_est <- para$xi; epsilon_est <- para$epsilon ; eta_est <- para$eta; delta_est = para$delta
    }
    if (floor(i/10000)==i/10000) { cat(i," iterations")}
  }
  
  return(list(xi = xi_est,epsilon = epsilon_est,eta = eta_est,delta = delta_est,
              samples_x = x_store,samples_y = z_store, mu_store = mu_store,
              var_store = var_store,skew_store = skew_store,kurt_store = kurt_store))
}



