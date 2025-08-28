library(ggplot2)

# SGD-MLE Algorithm 
# the default parameter is for Johnson distribution
sgd_mle <- function(data, 
                    initial_params = list(lambda = 1, delta = 1, xi = 1, gamma = 1),
                    gradient = johnson_log_gradient,
                    loglik_single = log_dJohnson,
                    force_change = function(params) {
                      params$lambda <- pmax(params$lambda, 0.001)  # lambda > 0
                      params$delta  <- pmax(params$delta, 0.001)   # delta > 0
                      return(params)
                    },
                    gamma = 0.9,
                    delta = 0.7,
                    tolerance = 1e-5,
                    batch_size = 100,
                    Max_iter = 20000,
                    verbose = TRUE) {
 
  data <- as.matrix(data)
  n <- nrow(data)
  n_params <- length(initial_params)
  params <- initial_params
  param_names <- names(params)
  stop_sign = FALSE
  
  # estimation_error <- c * (log(n) / n)^alpha
  # optimal_iterations <- 1 / estimation_error
  # max_epochs <- ceiling(optimal_iterations / n)
  max_epochs <- ceiling(Max_iter / n)
  n_iter <- floor(n/batch_size)
  
  loglik_history <- numeric(max_epochs)
  v <- list()
  for (name in param_names){v[[name]] = 0}
  parameter_store <- list()

  for (epoch in 1:max_epochs) {
    
    shuffled_indices <- sample(1:n)
    epoch_loglik <- 0
    loglik_store <- numeric(n_iter)

    for (i in 1:n_iter) {
      
      index <- shuffled_indices[((i-1)*batch_size+1):(i*batch_size)]
      x_i <- as.matrix(data[index, ])
      
      grad <- as.list(gradient(x_i, params))
      
      for (name in param_names) {
        v[[name]] <- 0.8 * v[[name]] + gamma /(i + (epoch - 1) * n_iter)^delta * grad[[name]]
        params[[name]] <- params[[name]] - v[[name]]
        if (!is.finite(max(params[[name]])) || !is.finite(min(params[[name]])))  {
          if(i > 1) params[[name]] <- parameter_store[[i + n_iter *(epoch-1)-1]][[name]]
          if(i == 1 & epoch == 1) params[[name]] <- initial_params[[name]]
      }
      }
      params <- force_change(params)
      parameter_store[[i + n_iter *(epoch-1)]] <- params
      

      for(s in (1:batch_size)){
        epoch_loglik <- epoch_loglik - do.call(loglik_single, c(list(x_i[s,]), params))
      }
      
      loglik_store[i] <- epoch_loglik / (i*batch_size)
      
      if (!is.finite(epoch_loglik )|| epoch_loglik > 1e200 ) {
        stop_sign = TRUE
        #print(epoch_loglik)
        if(i == 1 & epoch == 1) params = initial_params 
        if(epoch >1 ) params =  parameter_store[[i + n_iter *(epoch-1)-1]]
        break }
      # if (i > 1) {
      #   denom <- max(abs(loglik_store[i-1]), .Machine$double.eps)
      #   improvement <- abs(loglik_store[i] - loglik_store[i-1]) / denom
      #   #improvement <- max(grad[[1]],grad[[2]],grad[[3]],grad[[4]])
      #   if (!is.finite(epoch_loglik)|| epoch_loglik > 1e25 || improvement < tolerance) {
      #     print(i)
      #     break
      #   }
      # }
    }
    
    if(stop_sign){break}
    loglik_history[epoch] <- epoch_loglik / (n_iter * batch_size)
    
    if (epoch > 1) {
      denom <- max(abs(loglik_history[epoch-1]), .Machine$double.eps)
      improvement <- abs(loglik_history[epoch-1] - loglik_history[epoch]) / denom
      if (!is.finite(epoch_loglik)|| epoch_loglik > 1e25 ||improvement < tolerance) {
        if (verbose) cat(sprintf("Converged at epoch %d (improvement: %.2e)\n", epoch, improvement))
        break
      }
    }
    # cat("epoch",epoch)
  }
  result <- list(
    parameters = as.data.frame(params),
    grad = as.data.frame(t(unlist(grad))),
    final_loglik = loglik_history[epoch],
    convergence_history = list(loglik = loglik_history[1:epoch]),
    n_epochs = epoch,
    n_samples = n
  )
  return(result)
}

