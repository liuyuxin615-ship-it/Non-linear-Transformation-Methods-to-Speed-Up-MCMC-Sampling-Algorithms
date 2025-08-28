
# define the function to draw the contour plot
court_plot <- function(logpi,samples,var1,var2,title,length.out = 200){
  
  x1_min = min(samples[,1]); x1_max = max(samples[,1]);x1_len = x1_max - x1_min
  x2_min = min(samples[,2]); x2_max = max(samples[,2]);x2_len = x2_max - x2_min
  x1 <- seq((x1_min-0.1*x1_len), (x1_max + 0.1*x1_len),length.out = length.out)
  x2 <- seq((x2_min-0.1*x2_len), (x2_max + 0.1*x2_len),length.out = length.out)
  g <- expand.grid(x1 = x1, x2 = x2)
  g$logpost <- as.numeric(apply(g[, c("x1", "x2")], 1, logpi))
  g$logpost_shift <- g$logpost - max(g$logpost)
  g$post <- exp(g$logpost_shift)

  p <- ggplot(g, aes(x = x1, y =x2, z = post)) +
    geom_contour(color = "red", bins = 12) +
    labs(title = title,x = var1,y = var2) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    )
  return(p)
}


realcase_visualization <- function(logpi,x_curr = c(0,0),var1= "beta0",var2="beta1",
                                   nits1 = 50000, nits2 = 5000,h = 0.9){


  # get hte dimension
  d <- length(x_curr)
  # sampling from the original distribution
  sampling_ori <- Adaptive_RWM(logpi,nits1,h,x_curr)
  samples_ori <- sampling_ori$x_store
  p_ori <- court_plot(logpi,samples_ori,var1,var2,title = "Original distribution")
  
  # Get the linear transformation density(X^* = \Sigma^-0.5*X)
  A <- cov(samples_ori) ; eig <- eigen(A) ; V <- eig$vectors ; D <- diag(eig$values)
  # A^(-0.5) = V * D^(-0.5) * V^T
  A_neg_half <- V %*% diag(1 / sqrt(eig$values)) %*% t(V)
  A_half <- solve(A_neg_half)
  logpi_linear_trans <- function(x){logpi(A_half%*% x)}
  sampling_linear_trans <- Adaptive_RWM(logpi_linear_trans,nits1,h,x_curr)
  samples_linear_trans <- sampling_linear_trans$x_store
  p_linear_trans <- court_plot(logpi_linear_trans,samples_linear_trans,
                               var1,var2,title = "Distribution after linear transformation")
  
  # Get the density after the transformation map
  # With Johnson-SU distribution
  #using MLE to find the approximation parameters
  params <- fit_johnsonsu_MLE(samples_ori)
  # get the approximated parameters
  xi <- params$xi
  lambda <- params$lambda
  gamma <- params$gamma
  delta <- params$delta
  
  # define the transformation map of Johnson distribution
  x_from_y <- function(y) sinh((y - gamma)/delta)*lambda + xi 
  # Define the transformed density p_Y(y) in Johnson distribution
  log_p_Y <- function(y) {
    x <- x_from_y(y)
    logf <- logpi(x) + sum(log(1 + ((x - xi) / lambda)^2 ) /2)
    logf <- max(logf,-1e30)
    if(is.nan(logf)) logf <- -1e30
    return(logf)
  }
  # Sampling in the transformed space
  john_sampling <- Adaptive_RWM(log_p_Y,nits1,h ,x_curr=rep(0,d))
  john_samples <- john_sampling$x_store
  A <- cov(john_samples) ; eig <- eigen(A) ; V <- eig$vectors ; D <- diag(eig$values)
  # A = V * D * V^T
  A_half <- V %*%  diag(sqrt(eig$values)) %*% t(V)
  logp_p_Y_linear <- function(x){log_p_Y (A_half%*% x)}
  john_samples <- Adaptive_RWM(logp_p_Y_linear,nits2,h ,x_curr=rep(0,d))$x_store
  p_john <- court_plot(logp_p_Y_linear,john_samples,var1,var2,
                       title = "Distribution after Johnson SU transformation")
  
  
  # Skew-normal transformation
  # using MLE to find the SN approximation distribution
  sn_approx <- fit_sn_MLE (samples_ori)
  # get the each parameter
  xi <- sn_approx$ xi
  omega <- sn_approx$ omega
  alpha <- sn_approx$ alpha
  delta <- alpha / sqrt(1+alpha^2)

  # Get A,B,C
  A <- omega * sqrt(1 - delta^2)
  B <- omega * delta
  C <- xi
  # get the density after transformation
  log_piy <- function(y,v){
    x <- A * y + B * rep(abs(v),d) + C
    # print(x)
    logf <- logpi(x) - v^2/2
    logf <- max(logf,-1e30)
    return(logf)
  }
  sn_sampling <- newadaptiveMCMC(log_piy,nits = nits1,h,x_curr =rep(0,d))
  sn_samples <- sn_sampling$x_store
  A_sn <- cov(sn_samples) ; eig <- eigen(A_sn) ; V <- eig$vectors ; D <- diag(eig$values)
  # A = V * D * V^T
  A_negative_half <- V %*%  diag(1/sqrt(eig$values)) %*% t(V)
  p_sn <-  ggplot()+
    geom_density_2d(data = as.data.frame(sn_samples%*%A_negative_half),aes(x = V1, y = V2), color = "red") +
    labs(title = "Distribution after Skew-normal transformation", x = var1, y = var2) +
    theme_minimal() +  # 保留文本大小
    theme(
      panel.background = element_rect(fill = "white", color = NA),   
      plot.background  = element_rect(fill = "white", color = NA),   
      panel.grid       = element_blank(),                           
      panel.border     = element_blank(),                         
      axis.line        = element_blank()                            
    )
  
  #using Sinh-arcsinh approximation
  # get the approximatio  parameter
  sinh_arcsin_approx <- fit_highdim_sinh_arcsinh(samples_ori,"MLE")
  xi <- sinh_arcsin_approx$xi
  eta <- sinh_arcsin_approx$eta
  epsilon <- sinh_arcsin_approx$epsilon
  delta <- sinh_arcsin_approx$delta
  
  # get the transformation map
  x_from_z <- function(z) xi + eta * sinh((asinh(z)+epsilon)/delta)
  y_from_x <- function(x) (x - xi) / eta
  
  # Define the transformed density p_Y(z)
  log_p_z <- function(z) {
    x <- x_from_z(z)
    y <- y_from_x(x)
    
    logf <- logpi(x) - sum(log(1 + z^2) / 2 + log(1 + y^2) / 2)
    logf <- max(logf,-1e30)
    if(is.nan(logf)) logf <- -1e30
    return(logf)
  }
  
  # Sampling from log_o_Y
  sinh_sampling <- Adaptive_RWM(log_p_z,nits = nits1,x_curr = rep(0,d),h = h)
  sinh_samples <- sinh_sampling$x_store
  A <- cov(sinh_samples) ; eig <- eigen(A) ; V <- eig$vectors ; D <- diag(eig$values)
  # A = V * D * V^T
  A_half <- V %*%  diag(sqrt(eig$values)) %*% t(V)
  log_p_z_linear <- function(z){log_p_z (A_half%*% z)}
  sinh_samples <- Adaptive_RWM(log_p_z_linear,nits2,h ,x_curr=rep(0,d))$x_store

  # the countour plot
  p_sinh <- court_plot(log_p_z_linear ,sinh_samples,var1,var2,
                       title = "Distribution after Sinh-arcsinh transformation")
  
  return(list(p_ori = p_ori,p_linear_trans = p_linear_trans,
              p_john = p_john,p_sn = p_sn,p_sinh = p_sinh))
  
}


# the posterior of logistic regression
set.seed(4231)
# True parameter values
beta1 <- -5 ; beta2 <- 3
# Number of observations
n_logistic <- 1000
# Generate covariates
x_logistic <- rnorm(n_logistic, mean = -2, sd = 1)
# Create success probabilities
p_logistic <- exp(beta1 + beta2*x_logistic)/(1+exp(beta1 + beta2*x_logistic))
# Generate y values
y_logistic <- rbinom(n_logistic, size = 1, prob = p_logistic)
mean(y_logistic) # when the mean is around 0.001 the posterior is pretty skew
# Create posterior
logpi_logistic <- function(beta) {
  logprior <- dnorm(beta[1], mean = 0, sd = 10, log = T) + dnorm(beta[2], mean = 0, sd = 10, log = T)
  p <- 1/(1+exp( - beta[1] - beta[2]*x_logistic))
  loglik <- sum(dbinom(y_logistic, size = 1, prob = p, log = T))
  logpi <- pmin(loglik + logprior,1e30)
  return(pmax(logpi,-1e300))
}

d_logpi_logistic <- function(beta){
  beta1 <- beta[1]; beta2 <- beta[2]
  p_term <- exp( - beta[1] - beta[2]*x_logistic)
  p <- 1/(1+p_term)
  term <- pmax(pmin(y_logistic / p -(1 - y_logistic)/(1-p),1e300),-1e300)
  dp_dbeta1 <- -p_term / (1+p_term)^2
  dp_dbeta2 <- -x_logistic * p_term / (1+p_term)^2
  grad1 <- -beta1 / 100 + sum(term * dp_dbeta1)
  grad2 <- -beta2 / 100 + sum(term * dp_dbeta2)
  return(pmax(c(grad1,grad2),-1e300))
}

#Neal's funnel distribution
log_funnel <- function(para) {
  x <- para[1]; y = para[2]
  - 0.5 * y^2/5 - 0.5 * exp(-y) * x^2
}
d_log_funnel <- function(para){
  x <- para[1]; y = para[2]
  grad1 <- -exp(-y)*x
  grad2 <- -y / 5 + 0.5 * exp(-y) * x^2
  return(c(grad1,grad2))
}


#Rosenbrock distribution
log_rosenbrock <- function(para) {
  x <- para[1]; y = para[2]
  return( (- ((y - x^2)^2)*100- ((1-x)^2) )/20)
}
d_log_rosenbrock <- function(para){
  x <- para[1];y = para[2]
  grad1 <- 2*x*(y - x^2) - x
  grad2 <- -(y - x^2)
  return(c(grad1,grad2))
}




# posterior for high dimension posterior of logistic regression
# the posterior of logistic regression
set.seed(4231)
# dimension 
d = 20
# True parameter values
beta <- c(-5,rep(1,d-1))
# Number of observations
n_logistic_highdim <- 1000
# Generate covariates
x_logistic_highdim <- cbind(rep(1,n_logistic_highdim), 
                            rmvnorm(n_logistic_highdim, mean =rep(-0.5,(d-1)),sigma = diag((d-1))))
# Create success probabilities
p_logistic_highdim <- exp(x_logistic_highdim %*% beta)/(1+exp(x_logistic_highdim %*% beta))
# Generate y values
y_logistic_highdim <- rbinom(n_logistic_highdim, size = 1, prob = p_logistic_highdim)
mean(y_logistic_highdim) # when the mean is around 0.001 the posterior is pretty skew
# Create posterior
logpi_logistic_highdim <- function(beta) {
  logprior <- dmvnorm(beta,mean = rep(0,d),sigma = 100* diag(d),log = T)
  p <- 1/(1+exp( - x_logistic_highdim %*% beta))
  loglik <- sum(dbinom(y_logistic_highdim, size = 1, prob = p, log = T))
  logpi <- pmin(loglik + logprior,1e30)
  return(pmax(logpi,-1e300))
}

d_logpi_logistic_highdim <- function(beta){
  term1 <- exp(- x_logistic_highdim %*% beta)
  term2 <- exp(x_logistic_highdim %*% beta)
  grad <- -beta/100 + t(x_logistic_highdim)%*% (y_logistic_highdim*(term1/(1+ term1))) -
    t(x_logistic_highdim)%*%((1-y_logistic_highdim)*(term2/(1+term2)))
  return(grad)
}


# Neal's funel distribution in high diemnsion
log_funnel_highdim <- function(para) {
  d <- length(para)
  x <- para[1]; y = para[2:d]
  - 0.5 * x^2/5 - sum(0.5 * exp(-x) * y^2)
}
d_log_funnel_highdim <- function(para){
  d <- length(para)
  x <- para[1]; y = para[2:d]
  grad1 <- -x / 5 + sum(0.5 * exp(-x) * y^2)
  grad2 <- -exp(-x)*y
  return(c(grad1,grad2))
}

