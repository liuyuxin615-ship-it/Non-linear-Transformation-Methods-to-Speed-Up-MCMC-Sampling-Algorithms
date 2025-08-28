# all parameters input in gradient function is a list.
# Johnson SU distribution
library(coda)
log_dJohnson <- function(x,lambda,delta,xi,gamma){
  z <- (x - xi) / lambda
  term1 <- log(delta) - log(lambda) - log(2 * pi)/2
  term2 <- -0.5 * log(1 + z^2)
  term3 <- -0.5 * (gamma + delta * asinh(z))^2
  log_pdf <- sum(term1 + term2 + term3)
  #print(term1+term2+term3)
  return(max(log_pdf,-1e30))
}
# define the gradient for x of logpi
d_logJohnson <- function(x,lambda,delta,xi,gamma){
  z <- (x - xi) / lambda
  result <- - 1 / lambda * (z / (1 + z^2) + 1 / sqrt(1 + z^2) * delta *(gamma + delta * asinh(z)))
  return(result)
}
# define the gradient for each parameter for Johnson distribution
johnson_log_gradient <- function(x,params){
  # params = c(lambda,delta,xi,gamma)
  x = as.matrix(x);n_col= ncol(x)
  lambda = params$lambda; delta = params$delta; xi = params$xi; gamma = params$gamma
  grad <- list()
  z <- matrix((apply(x,1,function(x) (x - xi) / lambda)),ncol =n_col,byrow = TRUE)
  grad$lambda <- colMeans(matrix(apply(z,1,function(z)-(-1/lambda + z^2 / lambda / (1 + z^2) + (gamma + delta*asinh(z))*
                                                     delta * z / lambda / sqrt(1 + z^2))),ncol = n_col,byrow = TRUE))
  grad$delta <-colMeans(matrix(apply(z,1,function(z)-(1 / delta -asinh(z) * (gamma + delta * asinh(z)))),ncol = n_col,byrow = TRUE))
  grad$xi <- colMeans(matrix(apply(z,1,function(z)-(z / lambda / (1 + z^2) + 
                                                 (gamma + delta * asinh(z)) * delta / (lambda * sqrt(1 + z^2)))),ncol = n_col,byrow = TRUE))
  grad$gamma <- colMeans(matrix(apply(z,1,function(z) gamma + delta * asinh(z)),ncol = n_col,byrow = TRUE))
  return(grad)
  # grad = (d(-log_pix) / d(lambda),d(-log_pix) / d(delta),
  #         d(-log_pix) / d(xi),d(-log_pix) / d(gamma))
}


# Skew-normal distribution
log_dsn <- function(x,xi,omega,alpha){
  pi <- 3.1415926
  z <- (x - xi)/omega
  term1 <- -z^2/2
  term2 <- log(pnorm(alpha * z))
  term3 <- log(2) - log(omega) - log(2*pi)/2
  logf <- sum(term1 + term2 + term3)
  return(pmax(logf,-1e30))
}
d_logsn <- function(x,xi,omega,alpha){
  z <- (x - xi) / omega
  term1 <- -z / omega
  term2 <- dnorm(alpha * z) / pnorm(alpha * z) * alpha / omega
  result <- term1 + term2
  return(result)
}
sn_log_gradient <- function(x,params){
  # params = list(xi,omega,alpha)
  x = as.matrix(x);n_row = nrow(x);n_col = ncol (x)
  xi = params$xi;omega <- params$omega;alpha <- params$alpha
  z <- matrix(apply(x,1,function(x) (x - xi) / omega),ncol = n_col,byrow = TRUE)
  
  grad <- list()
  grad$xi <- colMeans(matrix(apply(z,1,function(z) -(z / omega- exp(dnorm(alpha * z, log = TRUE) - pnorm(alpha * z, log.p = TRUE))
                                                * alpha / omega)),ncol = n_col,byrow = TRUE))
  grad$omega <- colMeans(matrix(apply(z,1,function(z) -(-1/omega + z^2 / omega + 
                                                     exp(dnorm(alpha * z, log = TRUE) - pnorm(alpha * z, log.p = TRUE)) *
                                                     (-alpha * z / omega))),ncol = n_col,byrow = TRUE))
  grad$alpha <- colMeans(matrix(apply(z,1,function(z)-(exp(dnorm(alpha * z, log = TRUE) - pnorm(alpha * z, log.p = TRUE)) * z)),ncol = n_col,byrow = TRUE))
  return(grad)
}


# Sinh-arcsinh distribution
log_dsinh_arcsinh <- function(x,xi,epsilon,eta,delta){
  pi <- 3.1415926535
  y <- (x - xi) / eta
  z <- pmax(pmin(sinh(delta * asinh(y) - epsilon),1e10),-1e10)
  term1 <- log(delta) - log(eta)
  term2 <- 0.5 * (log(1 + z^2) - log(2*pi*(1+y^2)))
  term3 <- -0.5 * z^2
  log_pdf <- sum(term1 + term2 + term3)
  return(max(log_pdf,-1e30))
}

d_logsinh <- function(x,xi,epsilon,eta,delta){
  y <- (x - xi) / eta
  z <- pmax(pmin(sinh(delta * asinh(y) - epsilon),1e10),-1e10)
  dz_dx <- 1 / eta * cosh(asinh(y) * delta - epsilon) * delta  / sqrt(1+y^2)
  dy_dx <- 1 / eta
  result <- - z / (1 + z^2) * dz_dx - y / (1 + y^2) * dy_dx - z* dz_dx
  return(result)
}

sinh_log_gradient <- function(x,params){
  
  xi <- params$xi; epsilon <- params$epsilon; eta <- params$eta; delta <- params$delta
  x <- as.matrix(x); n_row <- nrow(x); n_col <- ncol(x)
  if(n_col == 1){
    y <- as.matrix(apply(x,1,function(x) (x - xi) / eta))
    z <- as.matrix(apply(y,1,function(y) pmax(pmin(sinh(delta * asinh(y) - epsilon),1e10),-1e10)))
    term <- as.matrix(apply(y,1,function(y) cosh(asinh(y) * delta - epsilon)))
    term_z <- as.matrix(apply(z,1,function(z) z / (1 + z^2)))
    term_y <- as.matrix(apply(y,1,function(y) y / (1 + y^2)))
  }else{
    y <- t(apply(x,1,function(x) (x - xi) / eta))
    z <- t(apply(y,1,function(y) pmax(pmin(sinh(delta * asinh(y) - epsilon),1e10),-1e10)))
    term <- t(apply(y,1,function(y) cosh(asinh(y) * delta - epsilon)))
    term_z <- t(apply(z,1,function(z) z / (1 + z^2)))
    term_y <- t(apply(y,1,function(y) y / (1 + y^2)))
  }
  grad_delta = grad_eta = grad_epsilon = grad_xi =  matrix(0,ncol = n_col,nrow = n_row)
  grad <- list()
  for(i in 1:n_row){
    dz_dy <- term[i,] * delta / sqrt(1 + y[i,]^2)
    grad_delta[i,] <- -(1 / delta + term_z[i,] * term[i,] * asinh(y[i,]) - z[i,] * term[i,] * asinh(y[i,]))
    grad_eta[i,] <- -(-1 / eta + term_z[i,] * dz_dy * (-(x[i,]- xi) / eta^2) 
                      - term_y[i,] * (-(x[i,]- xi) / eta^2) - z[i,] * dz_dy *(-(x[i,]- xi) / eta^2))
    grad_epsilon[i,] <- -(-term_z[i,] * term[i,] + z[i,] * term[i,])
    grad_xi[i,] <- -(term_z[i,] * dz_dy * (-1/eta) - term_y[i,] * (-1/eta) - z[i,]* dz_dy * (-1 / eta))
  }
  grad$xi <- colMeans(grad_xi)
  grad$epsilon <- colMeans(grad_epsilon)
  grad$eta <- colMeans(grad_eta)
  grad$delta <- colMeans(grad_delta)
  return(grad)
} 





