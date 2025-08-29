library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)

# in one dimension case compare the algorithm's efficiency on sinh-arcsinh distribution across different skewness
compare_skewness_algorithms <- function(epsilon_values, delta_sinh=1,xi_sinh = 1,eta_sinh =1,
                                        nits = 50000, pre_nits = 5000) {
  # delta is the parameter to control the tail behaviour
  # Other parameters

  h <- 0.9; x_curr <- 0
  c_mom <- 500; delta_mom <- 0.5; t1_mom <- 3000
  c_mle <- 200; delta_mle <- 0.9; t1_mle <- 5000
  
  # Define algorithm names
  algorithm_names <- c("epsilon_sinh",
                       "Sinh_mom_ESS", "Sinh_mom_time", "Sinh_mom_ESS_per_time",
                       "Sinh_MLE_ESS", "Sinh_MLE_time", "Sinh_MLE_ESS_per_time", 
                       "Sinh_online_MLE_ESS", "Sinh_online_MLE_time", "Sinh_online_MLE_ESS_per_time",
                       "Sinh_online_mom_ESS", "Sinh_online_mom_time", "Sinh_online_mom_ESS_per_time",
                       "John_MLE_ESS", "John_MLE_time", "John_MLE_ESS_per_time",
                       "John_mom_ESS", "John_mom_time", "John_mom_ESS_per_time",
                       "John_online_MLE_ESS", "John_online_MLE_time", "John_online_MLE_ESS_per_time",
                       "John_online_mom_ESS", "John_online_mom_time", "John_online_mom_ESS_per_time",
                       "SN_MLE_ESS", "SN_MLE_time", "SN_MLE_ESS_per_time",
                       "SN_mom_ESS", "SN_mom_time", "SN_mom_ESS_per_time",
                       "SN_online_MLE_ESS", "SN_online_MLE_time", "SN_online_MLE_ESS_per_time",
                       "SN_online_mom_ESS", "SN_online_mom_time", "SN_online_mom_ESS_per_time",
                       "Adaptive_MALA_ESS", "Adaptive_MALA_time", "Adaptive_MALA_ESS_per_time",
                       "Adaptive_RWM_ESS", "Adaptive_RWM_time", "Adaptive_RWM_ESS_per_time",
                       "Perfect_Map_ESS", "Perfect_Map_time", "Perfect_Map_ESS_per_time")
  
  # Initialize comparison dataframe
  compare_eff <- matrix(0, nrow = length(epsilon_values), ncol = length(algorithm_names))
  colnames(compare_eff) <- algorithm_names
  compare_eff <- data.frame(compare_eff)
  compare_eff[, "epsilon_sinh"] <- epsilon_values
  
  
  # Loop through each epsilon value
  for (i in (1:length(epsilon_values))) {

    epsilon_sinh <- epsilon_values[i]
    cat("Current epsilon_sinh =", epsilon_sinh, "\n")
    
    # Define target distribution functions
    d_logpi <- function(x) d_logsinh(x, xi_sinh, epsilon_sinh, eta_sinh, delta_sinh)
    log_pi <- function(x) log_dsinh_arcsinh(x, xi_sinh, epsilon_sinh, eta_sinh, delta_sinh)
    
    # Generate pre-samples for transformation methods
    start_pre <- Sys.time()
    pre_sampling <- Adaptive_RWM(log_pi, pre_nits,h, x_curr,target_a = 0.44)
    end_pre <- Sys.time()
    pre_samples <- as.matrix(pre_sampling$x_store[floor(pre_nits/3):pre_nits])
    pre_time <- as.numeric(end_pre - start_pre, units = "secs")
    
    # 1. Sinh-arcsinh moment method
    start_time <- Sys.time()
    trans_sinh_mom <- Sampling_trsanformation_sinh_arcsinh(pre_samples,log_pi, nits, method = "Moment")
    sinh_mom_samples <- trans_sinh_mom$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "Sinh_mom_ESS"] <- round(effectiveSize(mcmc(sinh_mom_samples)), 0)
    compare_eff[i, "Sinh_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Sinh-arcsinh moment method")
    
    # 2. Sinh-arcsinh MLE method
    start_time <- Sys.time()
    trans_sinh_mle <- Sampling_trsanformation_sinh_arcsinh(pre_samples,log_pi, nits, method = "MLE")
    sinh_mle_samples <- trans_sinh_mle$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "Sinh_MLE_ESS"] <- round(effectiveSize(mcmc(sinh_mle_samples)), 0)
    compare_eff[i, "Sinh_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Sinh-arcsinh MLE method")
    
    # 3. Sinh-arcsinh online MLE
    start_time <- Sys.time()
    trans_sinh_online_mle <- Transformation_sinh_arcsinh_MLE(log_pi, nits, x_curr)
    sinh_online_mle_samples <- trans_sinh_online_mle$samples_y[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "Sinh_online_MLE_ESS"] <-round(effectiveSize(mcmc(sinh_online_mle_samples[is.finite(sinh_online_mle_samples)])), 0)
    compare_eff[i, "Sinh_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Sinh-arcsinh online MLE method")
    
    
    # 4. Sinh-arcsinh online moment
    start_time <- Sys.time()
    trans_sinh_online_mom <- Transformation_sinh_arcsinh_mom(log_pi, nits, x_curr)
    sinh_online_mom_samples <- trans_sinh_online_mom$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "Sinh_online_mom_ESS"] <- round(effectiveSize(mcmc(sinh_online_mom_samples)), 0)
    compare_eff[i, "Sinh_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Sinh-arcsinh online Moment method")
    
    # 5. Johnson-SU MLE
    start_time <- Sys.time()
    trans_john_mle <- Sampling_trsanformation_johnsonsu(pre_samples,log_pi, nits, method = "MLE")
    john_mle_samples <- trans_john_mle$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "John_MLE_ESS"] <- round(effectiveSize(mcmc(john_mle_samples)), 0)
    compare_eff[i, "John_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Johnson SU MLE method")
    
    
    # 6. Johnson-SU moment
    start_time <- Sys.time()
    trans_john_mom <- Sampling_trsanformation_johnsonsu(pre_samples,log_pi, nits, method = "Moment")
    john_mom_samples <- trans_john_mom$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "John_mom_ESS"] <- round(effectiveSize(mcmc(john_mom_samples)), 0)
    compare_eff[i, "John_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Johnson SU Moment method")
    
    # 7. Johnson-SU online MLE
    start_time <- Sys.time()
    trans_john_online_mle <- Transformation_johnsonsu_MLE(log_pi, nits, x_curr)
    john_online_mle_samples <- trans_john_online_mle$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "John_online_MLE_ESS"] <-round(effectiveSize(mcmc(john_online_mle_samples[is.finite(john_online_mle_samples)])), 0)
    compare_eff[i, "John_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Johnson SU online MLE method")
    
    # 8. Johnson-SU online moment
    start_time <- Sys.time()
    trans_john_online_mom <- Transformation_johnsonsu_mom(log_pi, nits, x_curr)
    john_online_mom_samples <- trans_john_online_mom$samples_y[floor(nits/2):nits]
    end_time <- Sys.time()
    john_online_mom_samples <- rbind(john_online_mom_samples[is.finite(john_online_mom_samples)],1)
    ESS <- ifelse(
      max(john_online_mom_samples) == min(john_online_mom_samples),0,
      round(effectiveSize(mcmc(john_online_mom_samples)), 0)
    )
    compare_eff[i, "John_online_mom_ESS"] <- ESS
    compare_eff[i, "John_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Johnson SU online moment method")
   
    # 9. Skew-Normal MLE
    start_time <- Sys.time()
    trans_sn_mle <- Sampling_trsanformation_sn(pre_samples,log_pi, nits, approximation = "MLE")
    sn_mle_samples <- trans_sn_mle$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "SN_MLE_ESS"] <- round(effectiveSize(mcmc(sn_mle_samples)), 0)
    compare_eff[i, "SN_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished skew-normal MLE method")
    
    # 10. Skew-Normal moment
    start_time <- Sys.time()
    trans_sn_mom <- Sampling_trsanformation_sn(pre_samples,log_pi, nits, approximation = "Moment")
    sn_mom_samples <- trans_sn_mom$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "SN_mom_ESS"] <- round(effectiveSize(mcmc(sn_mom_samples)), 0)
    compare_eff[i, "SN_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished skew-normal Moment method")
    
    # 11. Skew-Normal online MLE
    start_time <- Sys.time()
    trans_sn_online_mle <- Transformation_sn_MLE(log_pi, nits, x_curr)
    sn_online_mle_samples <- trans_sn_online_mle$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "SN_online_MLE_ESS"] <-round(effectiveSize(mcmc(sn_online_mle_samples[is.finite(sn_online_mle_samples)])), 0)
    compare_eff[i, "SN_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished skew-normal online MLE method")
    
    # 12. Skew-Normal online moment
    start_time <- Sys.time()
    trans_sn_online_mom <- Transformation_sn_mom(log_pi, nits, x_curr)
    sn_online_mom_samples <- trans_sn_online_mom$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "SN_online_mom_ESS"] <- round(effectiveSize(mcmc(sn_online_mom_samples[is.finite(sn_online_mom_samples)])), 0)
    compare_eff[i, "SN_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished skew-normal online moment method")
    
    # 13. Adaptive MALA
    start_time <- Sys.time()
    mala_samples <- Adaptive_MALA(d_logpi, log_pi, nits, h, x_curr)
    mala_final_samples <- mala_samples$x_store[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "Adaptive_MALA_ESS"] <- round(effectiveSize(mcmc(mala_final_samples)), 0)
    compare_eff[i, "Adaptive_MALA_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished adaptive MALA")
    
    # 14. Adaptive RWM
    start_time <- Sys.time()
    rwm_samples <- Adaptive_RWM(log_pi, nits,h, x_curr)
    rwm_final_samples <- rwm_samples$x_store[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "Adaptive_RWM_ESS"] <- round(effectiveSize(mcmc(rwm_final_samples)), 0)
    compare_eff[i, "Adaptive_RWM_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished adaptive RWM")
    
    # 15. Perfect Map (transformation map)
    start_time <- Sys.time()
    map_samples <- transformation_map_sinh_arcsinh(xi_sinh, epsilon_sinh, eta_sinh, delta_sinh, nits, log_pi)
    map_final_samples <- map_samples$samples_x[floor(nits/2):nits]
    end_time <- Sys.time()
    compare_eff[i, "Perfect_Map_ESS"] <- round(effectiveSize(mcmc(map_final_samples)), 0)
    compare_eff[i, "Perfect_Map_time"] <- as.numeric(end_time - start_time, units = "secs")
    
    print(paste("Completed epsilon =", epsilon_sinh))
  }

 
  # Calculate ESS per time for all methods
  methods <- c("Sinh_mom", "Sinh_MLE", "Sinh_online_MLE", "Sinh_online_mom",
               "John_MLE", "John_mom", "John_online_MLE", "John_online_mom", 
               "SN_MLE", "SN_mom", "SN_online_MLE", "SN_online_mom",
               "Adaptive_MALA", "Adaptive_RWM", "Perfect_Map")
  
  for (method in methods) {
    ess_col <- paste0(method, "_ESS")
    time_col <- paste0(method, "_time")
    ess_per_time_col <- paste0(method, "_ESS_per_time")
    compare_eff[, ess_per_time_col] <- round(compare_eff[, ess_col] / compare_eff[, time_col], 2)
  }
  #compare_eff
  # Prepare data for plotting
  ess_per_time_cols <- paste0(methods, "_ESS_per_time")
  plot_data <- compare_eff %>%
    select(epsilon_sinh, all_of(ess_per_time_cols)) %>%
    pivot_longer(cols = -epsilon_sinh, names_to = "method", values_to = "ESS_per_time")
  
  # Clean method names for legend
  plot_data <- plot_data %>%
    mutate(method = str_remove(method, "_ESS_per_time")) %>%
    mutate(method = recode(method,
                           "Sinh_mom" = "Sinh_batch_mom",
                           "Sinh_MLE" = "Sinh_batch_MLE", 
                           "Sinh_online_MLE" = "Sinh_incremental_MLE",
                           "Sinh_online_mom" = "Sinh_incremental_mom",
                           "John_MLE" = "John_batch_MLE",
                           "John_mom" = "John_batch_mom",
                           "John_online_MLE" = "John_incremental_MLE", 
                           "John_online_mom" = "John_incremental_mom",
                           "SN_MLE" = "SN_batch_MLE",
                           "SN_mom" = "SN_batch_mom",
                           "SN_online_MLE" = "SN_incremental_MLE",
                           "SN_online_mom" = "SN_incremental_mom", 
                           "Adaptive_MALA" = "Adaptive MALA",
                           "Adaptive_RWM" = "Adaptive RWM",
                           "Perfect_Map" = "Perfect Map"))
  
  p1 <- ggplot(plot_data, aes(x = epsilon_sinh, y = ESS_per_time, color = method, shape = method)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(linewidth = 1, alpha = 0.7) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Algorithm Efficiency Comparison: ESS per Second vs Epsilon\n (normal-like tail)",
      x = "Epsilon (Skewness Parameter)",
      y = "ESS per Second",
      color = "Sampling Algorithm",
      shape = "Sampling Algorithm"
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    scale_x_continuous(breaks = epsilon_values) +
    scale_shape_manual(values = rep(c(16, 17, 15, 18, 19, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2), length.out = 15)) +
    guides(color = guide_legend(ncol = 3),
           shape = guide_legend(ncol = 3))
  return(list(compare_eff = compare_eff, p = p1))
}

# at some fixed skewness level compare the algorithm's efficiency on sinh-arcsinh across different dimensions
compare_dimension_algorithms <- function(dimension_values,compare = min, epsilon_sinh = 1, delta_sinh = 1,
                                         xi_sinh = 1, eta_sinh =1, nits = 50000, pre_nits = 7000){
  h <- 0.9
  c_mom <- 500; delta_mom <- 0.5; t1_mom <- 3000
  c_mle <- 200; delta_mle <- 0.9; t1_mle <- 5000
  
  # Define algorithm names
  algorithm_names <- c("dimension",
                       "Sinh_mom_ESS", "Sinh_mom_time", "Sinh_mom_ESS_per_time",
                       "Sinh_MLE_ESS", "Sinh_MLE_time", "Sinh_MLE_ESS_per_time", 
                       "Sinh_online_MLE_ESS", "Sinh_online_MLE_time", "Sinh_online_MLE_ESS_per_time",
                       "Sinh_online_mom_ESS", "Sinh_online_mom_time", "Sinh_online_mom_ESS_per_time",
                       "John_MLE_ESS", "John_MLE_time", "John_MLE_ESS_per_time",
                       "John_mom_ESS", "John_mom_time", "John_mom_ESS_per_time",
                       "John_online_MLE_ESS", "John_online_MLE_time", "John_online_MLE_ESS_per_time",
                       "John_online_mom_ESS", "John_online_mom_time", "John_online_mom_ESS_per_time",
                       "SN_MLE_ESS", "SN_MLE_time", "SN_MLE_ESS_per_time",
                       "SN_mom_ESS", "SN_mom_time", "SN_mom_ESS_per_time",
                       "SN_online_MLE_ESS", "SN_online_MLE_time", "SN_online_MLE_ESS_per_time",
                       "SN_online_mom_ESS", "SN_online_mom_time", "SN_online_mom_ESS_per_time",
                       "Adaptive_MALA_ESS", "Adaptive_MALA_time", "Adaptive_MALA_ESS_per_time",
                       "Adaptive_RWM_ESS", "Adaptive_RWM_time", "Adaptive_RWM_ESS_per_time",
                       "Perfect_Map_ESS", "Perfect_Map_time", "Perfect_Map_ESS_per_time")
  
  # Initialize comparison dataframe
  compare_eff <- matrix(0, nrow = length(dimension_values), ncol = length(algorithm_names))
  colnames(compare_eff) <- algorithm_names
  compare_eff <- data.frame(compare_eff)
  compare_eff[, "dimension"] <- dimension_values
 
  # Loop through each epsilon value
  for (i in (1:length(dimension_values))) {
    
    d = dimension_values[i]
    x_curr = rep(0,d)
    xi_values = rep(xi_sinh,d); delta_values = rep(delta_sinh,d); 
    epsilon_values = rep(epsilon_sinh,d); eta_values = rep(eta_sinh,d)
    cat("Current dimension =", d, "\n")
    
    # Define target distribution functions
    d_logpi <- function(x) d_logsinh(x, xi_values, epsilon_values, eta_values, delta_values)
    log_pi <- function(x) log_dsinh_arcsinh(x, xi_values, epsilon_values, eta_values, delta_values)
    
    # Generate pre-samples for transformation methods
    start_pre <- Sys.time()
    pre_sampling <- Adaptive_RWM(log_pi, pre_nits,h, x_curr,target_a = 0.23)
    end_pre <- Sys.time()
    pre_samples <- as.matrix(pre_sampling$x_store[floor(pre_nits/3):pre_nits,])
    pre_time <- as.numeric(end_pre - start_pre, units = "secs")
    
    # 1. Sinh-arcsinh moment method
    start_time <- Sys.time()
    trans_sinh_mom <- Sampling_trsanformation_sinh_arcsinh(pre_samples,log_pi, nits, method = "Moment")
    sinh_mom_samples <- as.matrix(trans_sinh_mom$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "Sinh_mom_ESS"] <- compare(round(effectiveSize(mcmc(sinh_mom_samples)), 0))
    compare_eff[i, "Sinh_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Sinh-arcsinh moment method")
    
    # 2. Sinh-arcsinh MLE method
    start_time <- Sys.time()
    trans_sinh_mle <- Sampling_trsanformation_sinh_arcsinh(pre_samples,log_pi, nits, method = "MLE")
    sinh_mle_samples <- as.matrix(trans_sinh_mle$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "Sinh_MLE_ESS"] <- compare(round(effectiveSize(mcmc(sinh_mle_samples)), 0))
    compare_eff[i, "Sinh_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Sinh-arcsinh MLE method")
    
    # 3. Sinh-arcsinh online MLE
    start_time <- Sys.time()
    trans_sinh_online_mle <- Transformation_sinh_arcsinh_MLE(log_pi, nits, x_curr)
    sinh_online_mle_samples <- as.matrix(trans_sinh_online_mle$samples_y[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "Sinh_online_MLE_ESS"] <-compare(round(effectiveSize(mcmc(sinh_online_mle_samples[is.finite(sinh_online_mle_samples)])), 0))
    compare_eff[i, "Sinh_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Sinh-arcsinh online MLE method")
    
    
    # 4. Sinh-arcsinh online moment
    start_time <- Sys.time()
    trans_sinh_online_mom <- Transformation_sinh_arcsinh_mom(log_pi, nits, x_curr)
    sinh_online_mom_samples <- as.matrix(trans_sinh_online_mom$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "Sinh_online_mom_ESS"] <- compare(round(effectiveSize(mcmc(sinh_online_mom_samples)), 0))
    compare_eff[i, "Sinh_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Sinh-arcsinh online Moment method")
    
    # 5. Johnson-SU MLE
    start_time <- Sys.time()
    trans_john_mle <- Sampling_trsanformation_johnsonsu(pre_samples,log_pi, nits, method = "MLE")
    john_mle_samples <- as.matrix(trans_john_mle$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "John_MLE_ESS"] <- compare(round(effectiveSize(mcmc(john_mle_samples)), 0))
    compare_eff[i, "John_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Johnson SU MLE method")
    
    
    # 6. Johnson-SU moment
    start_time <- Sys.time()
    trans_john_mom <- Sampling_trsanformation_johnsonsu(pre_samples,log_pi, nits, method = "Moment")
    john_mom_samples <- as.matrix(trans_john_mom$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "John_mom_ESS"] <- compare(round(effectiveSize(mcmc(john_mom_samples)), 0))
    compare_eff[i, "John_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Johnson SU Moment method")
    
    # 7. Johnson-SU online MLE
    start_time <- Sys.time()
    trans_john_online_mle <- Transformation_johnsonsu_MLE(log_pi, nits, x_curr)
    john_online_mle_samples <- as.matrix(trans_john_online_mle$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "John_online_MLE_ESS"] <-compare(round(effectiveSize(mcmc(john_online_mle_samples[is.finite(john_online_mle_samples)])), 0))
    compare_eff[i, "John_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Johnson SU online MLE method")
    
    # 8. Johnson-SU online moment
    start_time <- Sys.time()
    trans_john_online_mom <- Transformation_johnsonsu_mom(log_pi, nits, x_curr)
    john_online_mom_samples <- as.matrix(trans_john_online_mom$samples_y[floor(nits/2):nits,])
    end_time <- Sys.time()
    john_online_mom_samples <- john_online_mom_samples[is.finite(john_online_mom_samples)]
    if(length(john_online_mom_samples==0)){
      ESS = 0
    }else{
      ESS <- apply(john_online_mom_samples,2,function(samples) ifelse(
        max(samples) == min(samples),rep(0),
        round(effectiveSize(mcmc(samples)), 0)
      ))
      }
    compare_eff[i, "John_online_mom_ESS"] <- compare(ESS)
    compare_eff[i, "John_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Johnson SU online moment method")
    
    # 9. Skew-Normal MLE
    start_time <- Sys.time()
    trans_sn_mle <- Sampling_trsanformation_sn(pre_samples,log_pi, nits, approximation = "MLE")
    sn_mle_samples <- as.matrix(trans_sn_mle$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "SN_MLE_ESS"] <- compare(round(effectiveSize(mcmc(sn_mle_samples)), 0))
    compare_eff[i, "SN_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished skew-normal MLE method")
    
    # 10. Skew-Normal moment
    start_time <- Sys.time()
    trans_sn_mom <- Sampling_trsanformation_sn(pre_samples,log_pi, nits, approximation = "Moment")
    sn_mom_samples <- as.matrix(trans_sn_mom$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "SN_mom_ESS"] <- compare(round(effectiveSize(mcmc(sn_mom_samples)), 0))
    compare_eff[i, "SN_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished skew-normal Moment method")
    
    # 11. Skew-Normal online MLE
    start_time <- Sys.time()
    trans_sn_online_mle <- Transformation_sn_MLE(log_pi, nits, x_curr)
    sn_online_mle_samples <- as.matrix(trans_sn_online_mle$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "SN_online_MLE_ESS"] <-compare(round(effectiveSize(mcmc(sn_online_mle_samples[is.finite(sn_online_mle_samples)])), 0))
    compare_eff[i, "SN_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished skew-normal online MLE method")
    
    # 12. Skew-Normal online moment
    start_time <- Sys.time()
    trans_sn_online_mom <- Transformation_sn_mom(log_pi, nits, x_curr)
    sn_online_mom_samples <- as.matrix(trans_sn_online_mom$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "SN_online_mom_ESS"] <- compare(round(effectiveSize(mcmc(sn_online_mom_samples[is.finite(sn_online_mom_samples)])), 0))
    compare_eff[i, "SN_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished skew-normal online moment method")
    
    # 13. Adaptive MALA
    start_time <- Sys.time()
    mala_samples <- Adaptive_MALA(d_logpi, log_pi, nits, h, x_curr)
    mala_final_samples <- as.matrix(mala_samples$x_store[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "Adaptive_MALA_ESS"] <- compare(round(effectiveSize(mcmc(mala_final_samples)), 0))
    compare_eff[i, "Adaptive_MALA_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished adaptive MALA")
    
    # 14. Adaptive RWM
    start_time <- Sys.time()
    rwm_samples <- Adaptive_RWM(log_pi, nits,h, x_curr)
    rwm_final_samples <- as.matrix(rwm_samples$x_store[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "Adaptive_RWM_ESS"] <- compare(round(effectiveSize(mcmc(rwm_final_samples)), 0))
    compare_eff[i, "Adaptive_RWM_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished adaptive RWM")
    
    # 15. Perfect Map (transformation map)
    start_time <- Sys.time()
    map_samples <- transformation_map_sinh_arcsinh(xi_values, epsilon_values, eta_values, delta_values, nits, log_pi)
    map_final_samples <- as.matrix(map_samples$samples_x[floor(nits/2):nits,])
    end_time <- Sys.time()
    compare_eff[i, "Perfect_Map_ESS"] <- compare(round(effectiveSize(mcmc(map_final_samples)), 0))
    compare_eff[i, "Perfect_Map_time"] <- as.numeric(end_time - start_time, units = "secs")
    
    print(paste("Completed dimension =", d))
  }
 
  # Calculate ESS per time for all methods
  methods <- c("Sinh_mom", "Sinh_MLE", "Sinh_online_MLE", "Sinh_online_mom",
               "John_MLE", "John_mom", "John_online_MLE", "John_online_mom", 
               "SN_MLE", "SN_mom", "SN_online_MLE", "SN_online_mom",
               "Adaptive_MALA", "Adaptive_RWM", "Perfect_Map")
  
  for (method in methods) {
    ess_col <- paste0(method, "_ESS")
    time_col <- paste0(method, "_time")
    ess_per_time_col <- paste0(method, "_ESS_per_time")
    compare_eff[, ess_per_time_col] <- round(compare_eff[, ess_col] / compare_eff[, time_col], 2)
  }
 
  # Prepare data for plotting
  ess_per_time_cols <- paste0(methods, "_ESS_per_time")
  plot_data <- compare_eff %>%
    select(dimension, all_of(ess_per_time_cols)) %>%
    pivot_longer(cols = -dimension, names_to = "method", values_to = "ESS_per_time")
  
  # Clean method names for legend
  plot_data <- plot_data %>%
    mutate(method = str_remove(method, "_ESS_per_time")) %>%
    mutate(method = recode(method,
                           "Sinh_mom" = "Sinh_batch_mom",
                           "Sinh_MLE" = "Sinh_batch_MLE", 
                           "Sinh_online_MLE" = "Sinh_incremental_MLE",
                           "Sinh_online_mom" = "Sinh_incremental_mom",
                           "John_MLE" = "John_batch_MLE",
                           "John_mom" = "John_batch_mom",
                           "John_online_MLE" = "John_incremental_MLE", 
                           "John_online_mom" = "John_incremental_mom",
                           "SN_MLE" = "SN_batch_MLE",
                           "SN_mom" = "SN_batch_mom",
                           "SN_online_MLE" = "SN_incremental_MLE",
                           "SN_online_mom" = "SN_incremental_mom", 
                           "Adaptive_MALA" = "Adaptive MALA",
                           "Adaptive_RWM" = "Adaptive RWM",
                           "Perfect_Map" = "Perfect Map"))

  p1 <- ggplot(plot_data, aes(x = dimension, y = ESS_per_time, color = method, shape = method)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(linewidth = 1, alpha = 0.7) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Algorithm Efficiency Comparison: \n Min ESS per Second vs Dimension(High Skew)",
      x = "Dimension",
      y = "ESS per Second (samples/sec)",
      color = "Sampling Algorithm",
      shape = "Sampling Algorithm"
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    scale_x_continuous(breaks = dimension_values) +
    scale_shape_manual(values = rep(c(16, 17, 15, 18, 19, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2), length.out = 15)) +
    guides(color = guide_legend(ncol = 3),
           shape = guide_legend(ncol = 3))
  return(list(compare_eff = compare_eff, p = p1))
}


# compare the different sample algorithm's efficiency on a spcific distribution
compare_efficiecny <- function(log_pi,d_logpi,dimension  = 2, n_rep = 3,nits = 50000, pre_nits = 5000){
  # delta is the parameter to control the tail behaviour
  # Other parameters
  
  h <- 0.9; x_curr <- rep(0,dimension)
  c_mom <- 500; delta_mom <- 0.5; t1_mom <- 3000
  c_mle <- 200; delta_mle <- 0.9; t1_mle <- 5000

  iteration <- (1:n_rep)
  # Define algorithm names
  algorithm_names <- c("iteration",
                       "Sinh_mom_ESS", "Sinh_mom_time", "Sinh_mom_ESS_per_time",
                       "Sinh_MLE_ESS", "Sinh_MLE_time", "Sinh_MLE_ESS_per_time", 
                       "Sinh_online_MLE_ESS", "Sinh_online_MLE_time", "Sinh_online_MLE_ESS_per_time",
                       "Sinh_online_mom_ESS", "Sinh_online_mom_time", "Sinh_online_mom_ESS_per_time",
                       "John_MLE_ESS", "John_MLE_time", "John_MLE_ESS_per_time",
                       "John_mom_ESS", "John_mom_time", "John_mom_ESS_per_time",
                       "John_online_MLE_ESS", "John_online_MLE_time", "John_online_MLE_ESS_per_time",
                       "John_online_mom_ESS", "John_online_mom_time", "John_online_mom_ESS_per_time",
                       "SN_MLE_ESS", "SN_MLE_time", "SN_MLE_ESS_per_time",
                       "SN_mom_ESS", "SN_mom_time", "SN_mom_ESS_per_time",
                       "SN_online_MLE_ESS", "SN_online_MLE_time", "SN_online_MLE_ESS_per_time",
                       "SN_online_mom_ESS", "SN_online_mom_time", "SN_online_mom_ESS_per_time",
                       "Adaptive_MALA_ESS", "Adaptive_MALA_time", "Adaptive_MALA_ESS_per_time",
                       "Adaptive_RWM_ESS", "Adaptive_RWM_time", "Adaptive_RWM_ESS_per_time")
  
  # Initialize comparison dataframe
  compare_eff <- matrix(0, nrow = length(iteration), ncol = length(algorithm_names))
  colnames(compare_eff) <- algorithm_names
  compare_eff <- data.frame(compare_eff)
  compare_eff[, "iteration"] <- iteration
  
  # Loop through each epsilon value
  for (i in (1:n_rep)) {
    
    cat("Current iteration", i, "\n")
    # Generate pre-samples for transformation methods
    start_pre <- Sys.time()
    pre_sampling <- Adaptive_RWM(log_pi, pre_nits,h, x_curr,target_a = 0.44)
    end_pre <- Sys.time()
    pre_samples <- as.matrix(pre_sampling$x_store[floor(pre_nits/3):pre_nits,])
    pre_time <- as.numeric(end_pre - start_pre, units = "secs")
    
    # 1. Sinh-arcsinh moment method
    start_time <- Sys.time()
    #trans_sinh_mom <- Sampling_trsanformation_sinh_arcsinh(pre_samples,log_pi, nits, method = "Moment")
    #sinh_mom_samples <- trans_sinh_mom$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "Sinh_mom_ESS"] <- min(round(effectiveSize(mcmc(sinh_mom_samples)), 0))
    compare_eff[i, "Sinh_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Sinh-arcsinh moment method")
    
    # 2. Sinh-arcsinh MLE method
    start_time <- Sys.time()
    trans_sinh_mle <- Sampling_trsanformation_sinh_arcsinh(pre_samples,log_pi, nits, method = "MLE")
    sinh_mle_samples <- trans_sinh_mle$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "Sinh_MLE_ESS"] <- min(round(effectiveSize(mcmc(sinh_mle_samples)), 0))
    compare_eff[i, "Sinh_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Sinh-arcsinh MLE method")
    
    # 3. Sinh-arcsinh online MLE
    start_time <- Sys.time()
    # trans_sinh_online_mle <- Transformation_sinh_arcsinh_MLE(log_pi, nits, x_curr)
    # sinh_online_mle_samples <- trans_sinh_online_mle$samples_y[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "Sinh_online_MLE_ESS"] <-0
      #min(round(effectiveSize(mcmc(sinh_online_mle_samples[is.finite(sinh_online_mle_samples)])), 0))
    compare_eff[i, "Sinh_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Sinh-arcsinh online MLE method")
    
    
    # 4. Sinh-arcsinh online moment
    start_time <- Sys.time()
    # trans_sinh_online_mom <- Transformation_sinh_arcsinh_mom(log_pi, nits, x_curr)
    # sinh_online_mom_samples <- trans_sinh_online_mom$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "Sinh_online_mom_ESS"] <- min(round(effectiveSize(mcmc(sinh_online_mom_samples)), 0))
    compare_eff[i, "Sinh_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Sinh-arcsinh online Moment method")
    
    # 5. Johnson-SU MLE
    start_time <- Sys.time()
    trans_john_mle <- Sampling_trsanformation_johnsonsu(pre_samples,log_pi, nits, method = "MLE")
    john_mle_samples <- trans_john_mle$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "John_MLE_ESS"] <- min(round(effectiveSize(mcmc(john_mle_samples)), 0))
    compare_eff[i, "John_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Johnson SU MLE method")
    
    
    # 6. Johnson-SU moment
    start_time <- Sys.time()
    # trans_john_mom <- Sampling_trsanformation_johnsonsu(pre_samples,log_pi, nits, method = "Moment")
    # john_mom_samples <- trans_john_mom$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "John_mom_ESS"] <- min(round(effectiveSize(mcmc(john_mom_samples)), 0))
    compare_eff[i, "John_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished Johnson SU Moment method")
    
    # 7. Johnson-SU online MLE
    start_time <- Sys.time()
    # trans_john_online_mle <- Transformation_johnsonsu_MLE(log_pi, nits, x_curr)
    # john_online_mle_samples <- trans_john_online_mle$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "John_online_MLE_ESS"] <-min(round(effectiveSize(mcmc(john_online_mle_samples[is.finite(john_online_mle_samples)])), 0))
    compare_eff[i, "John_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Johnson SU online MLE method")
    
    # 8. Johnson-SU online moment
    start_time <- Sys.time()
    # trans_john_online_mom <- Transformation_johnsonsu_mom(log_pi, nits, x_curr)
    # john_online_mom_samples <- trans_john_online_mom$samples_y[floor(nits/2):nits,]
    end_time <- Sys.time()
    john_online_mom_samples <- john_online_mom_samples[is.finite(john_online_mom_samples)]
    ESS <- ifelse(
      max(john_online_mom_samples) == min(john_online_mom_samples),0,
      round(effectiveSize(mcmc(john_online_mom_samples)), 0)
    )
    compare_eff[i, "John_online_mom_ESS"] <- min(ESS)
    compare_eff[i, "John_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished Johnson SU online moment method")
    
   
    # 9. Skew-Normal MLE
    start_time <- Sys.time()
    trans_sn_mle <- Sampling_trsanformation_sn(pre_samples,log_pi, nits, approximation = "MLE")
    sn_mle_samples <- trans_sn_mle$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "SN_MLE_ESS"] <- min(round(effectiveSize(mcmc(sn_mle_samples)), 0))
    compare_eff[i, "SN_MLE_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished skew-normal MLE method")
    
    # 10. Skew-Normal moment
    start_time <- Sys.time()
    trans_sn_mom <- Sampling_trsanformation_sn(pre_samples,log_pi, nits, approximation = "Moment")
    sn_mom_samples <- trans_sn_mom$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "SN_mom_ESS"] <- min(round(effectiveSize(mcmc(sn_mom_samples)), 0))
    compare_eff[i, "SN_mom_time"] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished skew-normal Moment method")
    
    # 11. Skew-Normal online MLE
    start_time <- Sys.time()
    # trans_sn_online_mle <- Transformation_sn_MLE(log_pi, nits, x_curr)
    # sn_online_mle_samples <- trans_sn_online_mle$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "SN_online_MLE_ESS"] <-min(round(effectiveSize(mcmc(sn_online_mle_samples[is.finite(sn_online_mle_samples)])), 0))
    compare_eff[i, "SN_online_MLE_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished skew-normal online MLE method")
    
    # 12. Skew-Normal online moment
    start_time <- Sys.time()
    # trans_sn_online_mom <- Transformation_sn_mom(log_pi, nits, x_curr)
    # sn_online_mom_samples <- trans_sn_online_mom$samples_x[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "SN_online_mom_ESS"] <- min(round(effectiveSize(mcmc(sn_online_mom_samples[is.finite(sn_online_mom_samples)])), 0))
    compare_eff[i, "SN_online_mom_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished skew-normal online moment method")
    
    # 13. Adaptive MALA
    start_time <- Sys.time()
    mala_samples <- Adaptive_MALA(d_logpi, log_pi, nits, h, x_curr)
    mala_final_samples <- mala_samples$x_store[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "Adaptive_MALA_ESS"] <- min(round(effectiveSize(mcmc(mala_final_samples)), 0))
    compare_eff[i, "Adaptive_MALA_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished adaptive MALA")
  
    
    
    # 14. Adaptive RWM
    start_time <- Sys.time()
    rwm_samples <- Adaptive_RWM(log_pi, nits,h, x_curr)
    rwm_final_samples <- rwm_samples$x_store[floor(nits/2):nits,]
    end_time <- Sys.time()
    compare_eff[i, "Adaptive_RWM_ESS"] <- min(round(effectiveSize(mcmc(rwm_final_samples)), 0))
    compare_eff[i, "Adaptive_RWM_time"] <- as.numeric(end_time - start_time, units = "secs")
    print("finished adaptive RWM")
  }
  methods <- c("Sinh_mom", "Sinh_MLE", "Sinh_online_MLE", "Sinh_online_mom",
               "John_MLE", "John_mom", "John_online_MLE", "John_online_mom", 
               "SN_MLE", "SN_mom", "SN_online_MLE", "SN_online_mom",
               "Adaptive_MALA", "Adaptive_RWM")
  n_method <- length(methods)
  eff_matrix <- as.data.frame(matrix(0,ncol = 14,nrow = n_rep))
  for (i in (1:n_method)) {
    method <- methods[i]
    ess_col <- paste0(method, "_ESS")
    time_col <- paste0(method, "_time")
    ess_per_time_col <- paste0(method, "_ESS_per_time")
    colnames(eff_matrix)[i] <- ess_per_time_col
    eff_matrix[, ess_per_time_col] <- round(compare_eff[, ess_col] / compare_eff[, time_col], 2)
  }

  return(list(eff_matrix = eff_matrix))
}





library(pracma)  

#define a function to calculate the total variation distance 
total_variation_continuous <- function(logdens_p, logdens_q, lower=-10, upper=10, n=10000) {
  x <- seq(lower, upper, length.out = n)
  f <- exp(sapply(x,logdens_p))
  g <- exp(sapply(x,logdens_q))
  return(0.5 * trapz(x, abs(f - g)))
}

# compare the total variation distance of surrogate distribution gotten by different algorithm with the sinh-arcsinh target distribution
compare_TVdis_algorithms <- function(epsilon_values,delta_sinh=1, xi_sinh = 1, eta_sinh =1,
                                     pre_nits = 5000,nits = 50000, t0 = 2000) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(transport)
  
  # 方法名字（3个分布 × 4个拟合方式 = 12）
  method_names <- c("Sinh_batch_Mom", "Sinh_batch_MLE", "Sinh_incremental_Mom", "Sinh_incremental_MLE",
                    "John_batch_Mom", "John_batch_MLE", "John_incremental_Mom", "John_incremental_MLE",
                    "SN_batch_Mom", "SN_batch_MLE", "SN_incremental_Mom", "SN_incremental_MLE")
  
  # 保存结果
  compare_res <- matrix(0, nrow = length(epsilon_values), ncol = 1 + 2*length(method_names))
  colnames(compare_res) <- c("epsilon", paste0(method_names, "_Wdist"), paste0(method_names, "_time"))
  compare_res <- as.data.frame(compare_res)
  compare_res$epsilon <- epsilon_values
  
  # 循环 epsilon
  for (i in (2:length(epsilon_values))) {
    
    epsilon_sinh <- epsilon_values[i]
    cat("Current epsilon =", epsilon_sinh, "\n")
    
    target_samples = xi_sinh + eta_sinh * sinh((asinh(rnorm(nits)) + epsilon_sinh) / delta_sinh)
    min_x <- min(target_samples);max_x <- max(target_samples)
    log_pi <- function(x) log_dsinh_arcsinh(x, xi_sinh, epsilon_sinh, eta_sinh, delta_sinh)
    
    start_pre <- Sys.time()
    pre_sampling <- Adaptive_RWM(log_pi, pre_nits, h, x_curr, target_a = 0.44)
    end_pre <- Sys.time()
    pre_samples <- as.matrix(pre_sampling$x_store[floor(pre_nits/3):pre_nits])
    pre_time <- as.numeric(end_pre - start_pre, units = "secs")
    
    #1. Sinh_batch_Mom
    method ="Sinh_batch_Mom"
    start_time <- Sys.time()
    para <- fit_highdim_sinh_arcsinh(pre_samples, method = "Moment")
    end_time = Sys.time()
    logd_sin_mom <- function(x) log_dsinh_arcsinh(x,para$xi,para$epsilon,para$eta,para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_sin_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sinh_batch_mom")
    
    #2. Sinh_batch_MLE
    method ="Sinh_batch_MLE"
    start_time <- Sys.time()
    para <- fit_highdim_sinh_arcsinh(pre_samples, method = "MLE")
    end_time = Sys.time()
    logd_sin_MLE<- function(x) log_dsinh_arcsinh(x,para$xi,para$epsilon,para$eta,para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_sin_MLE,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sinh_batch_mle")
    
    #3. Sinh_incremental_Mom
    method ="Sinh_incremental_Mom"
    start_time <- Sys.time()
    para <- Transformation_sinh_arcsinh_mom(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_sin_online_mom<- function(x) log_dsinh_arcsinh(x,para$xi,para$epsilon,para$eta,para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_sin_online_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished sinh_incremental_mom")
    
    #4. Sinh_incremental_MLE
    method ="Sinh_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_sinh_arcsinh_MLE(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_sin_online_MLE<- function(x) log_dsinh_arcsinh(x,para$xi,para$epsilon,para$eta,para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_sin_online_MLE,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished sinh_incremental_mle")
    
    #5. Johnson su_batch_Mom
    method ="John_batch_Mom"
    start_time <- Sys.time()
    sample_mean <- apply(pre_samples, 2, mean)
    sample_var <- apply(pre_samples, 2, var)
    sample_skew <- apply(pre_samples, 2, skewness)
    sample_kurt <- apply(pre_samples, 2, kurtosis)
    moment_matrix <- cbind(sample_mean, sample_var, sample_skew, sample_kurt)
    para <- fit_johnsonsu_moment(moment_matrix)
    end_time = Sys.time()
    logd_john_mom <- function(x) log_dJohnson(x,para$lambda,para$delta,para$xi,para$gamma)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_john_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished john_batch_mom")
    
    #6. Johnson su_batch_MLE
    method ="John_batch_MLE"
    start_time <- Sys.time()
    para <- fit_johnsonsu_MLE(pre_samples)
    end_time = Sys.time()
    logd_john_mle <- function(x) log_dJohnson(x,para$lambda,para$delta,para$xi,para$gamma)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_john_mle,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished john_batch_mle")
    
    #7. Johnson su_incremental_Mom
    method ="John_incremental_Mom"
    start_time <- Sys.time()
    para <- Transformation_johnsonsu_mom(log_pi,(nits -t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_john_online_mom <- function(x) log_dJohnson(x,para$lambda,para$delta,para$xi,para$gamma)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_john_online_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished john_incremental_mom")
    
    #8. Johnson su_incremental_MLE
    method ="John_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_johnsonsu_MLE(log_pi,(nits -t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_john_online_mle <- function(x) log_dJohnson(x,para$lambda,para$delta,para$xi,para$gamma)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_john_online_mle,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished john_incremental_mle")
    
    #9. SN_batch_Mom
    method ="SN_batch_Mom"
    start_time <- Sys.time()
    mu <- apply(pre_samples, 2, mean)
    variance <- apply(pre_samples, 2, var)
    skew <- apply(pre_samples, 2, skewness)
    para <- fit_sn_moment(mu, variance, skew)
    end_time = Sys.time()
    logd_sn_mom <- function(x) log_dsn(x,para$xi,para$omega,para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_sn_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sn_batch_mom")
    
    #10. SN_batch_MLE
    method ="SN_batch_MLE"
    start_time <- Sys.time()
    para <- fit_sn_MLE(pre_samples)
    end_time = Sys.time()
    logd_sn_mle <- function(x) log_dsn(x,para$xi,para$omega,para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_sn_mle,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sn_batch_mle")
    
    #11. SN_incremental_Mom
    method ="SN_incremental_Mom"
    start_time <- Sys.time()
    para <- Transformation_sn_mom(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_sn_online_mom <- function(x) log_dsn(x,para$xi,para$omega,para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_sn_online_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished sn_incremental_mom")
    
    #12. SN_incremental_MLE
    method ="SN_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_sn_MLE(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_sn_online_mle <- function(x) log_dsn(x,para$xi,para$omega,para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(total_variation_continuous(log_pi, logd_sn_online_mle,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished sn_incremental_mle")
  }
  
  
  # ---- 绘图 ----
  plot_data <- compare_res %>%
    select(epsilon, ends_with("_Wdist")) %>%
    pivot_longer(cols = -epsilon, names_to = "method", values_to = "Wdist") %>%
    mutate(method = gsub("_Wdist", "", method))
  
  p <- ggplot(plot_data, aes(x = epsilon, y = Wdist, color = method, shape = method)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(linewidth = 1, alpha = 0.7) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Total variance Distance vs Epsilon",
      x = "Epsilon (Skewness Parameter)",
      y = "Wasserstein Distance",
      color = "Method",
      shape = "Method"
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    scale_x_continuous(breaks = epsilon_values) +
    scale_shape_manual(values = rep(c(16, 17, 15, 18, 19, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2), length.out = 15)) +
    guides(color = guide_legend(ncol = 3),
           shape = guide_legend(ncol = 3))
  
  return(list(compare_res = compare_res, p = p))
}

# compare the total Wasserstein distance of surrogate distribution gotten by different algorithm with the sinh-arcsinh target distribution
compare_wasserstein_algorithms <- function(epsilon_values,dividemax = TRUE,
                                           f = function(x) return(x),
                                           delta_sinh=1, xi_sinh = 1, eta_sinh =1,
                                           pre_nits = 5000, nits = 50000,t0 = 2000) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(transport)
  
  # 方法名字（3个分布 × 4个拟合方式 = 12）
  method_names <- c("Sinh_batch_Mom", "Sinh_batch_MLE", "Sinh_incremental_Mom", "Sinh_incremental_MLE",
                    "John_batch_Mom", "John_batch_MLE", "John_incremental_Mom", "John_incremental_MLE",
                    "SN_batch_Mom", "SN_batch_MLE", "SN_incremental_Mom", "SN_incremental_MLE")
  
  # 保存结果
  compare_res <- matrix(0, nrow = length(epsilon_values), ncol = 1 + 2*length(method_names))
  colnames(compare_res) <- c("epsilon", paste0(method_names, "_Wdist"), paste0(method_names, "_time"))
  compare_res <- as.data.frame(compare_res)
  compare_res$epsilon <- epsilon_values
  
  # 循环 epsilon
  for (i in (1:length(epsilon_values))) {
    epsilon_sinh <- epsilon_values[i]
    cat("Current epsilon =", epsilon_sinh, "\n")
    
    target_samples = xi_sinh + eta_sinh * sinh((asinh(rnorm(nits)) + epsilon_sinh) / delta_sinh)
    
    d_logpi <- function(x) d_logsinh(x, xi_sinh, epsilon_sinh, eta_sinh, delta_sinh)
    log_pi <- function(x) log_dsinh_arcsinh(x, xi_sinh, epsilon_sinh, eta_sinh, delta_sinh)
    
    start_pre <- Sys.time()
    pre_sampling <- Adaptive_RWM(log_pi, pre_nits, h, x_curr, target_a = 0.44)
    end_pre <- Sys.time()
    pre_samples <- as.matrix(pre_sampling$x_store[floor(pre_nits/3):pre_nits])
    pre_time <- as.numeric(end_pre - start_pre, units = "secs")
    
    #1. Sinh_batch_Mom
    method ="Sinh_batch_Mom"
    start_time <- Sys.time()
    para <- fit_highdim_sinh_arcsinh(pre_samples, method = "Moment")
    end_time = Sys.time()
    sinh_mom_samples <- generate_sinh(nits, para$xi, para$epsilon, para$eta, para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(sinh_mom_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sinh_batch_mom")
    
    #2. Sinh_batch_MLE
    method ="Sinh_batch_MLE"
    start_time <- Sys.time()
    para <- fit_highdim_sinh_arcsinh(pre_samples, method = "MLE")
    end_time = Sys.time()
    sinh_mle_samples <- generate_sinh(nits, para$xi, para$epsilon, para$eta, para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(sinh_mle_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sinh_batch_mle")
    
    #3. Sinh_incremental_Mom
    method ="Sinh_incremental_Mom"
    start_time <- Sys.time()
    para <- Transformation_sinh_arcsinh_mom(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    sinh_online_mom_samples <- generate_sinh(nits, para$xi, para$epsilon, para$eta, para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(sinh_online_mom_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished sinh_incremental_mom")
    
    #4. Sinh_incremental_MLE
    method ="Sinh_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_sinh_arcsinh_MLE(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    sinh_online_mle_samples <- generate_sinh(nits, para$xi, para$epsilon, para$eta, para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(sinh_online_mle_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished sinh_incremental_mle")
    
    #5. Johnson su_batch_Mom
    method ="John_batch_Mom"
    start_time <- Sys.time()
    sample_mean <- apply(pre_samples, 2, mean)
    sample_var <- apply(pre_samples, 2, var)
    sample_skew <- apply(pre_samples, 2, skewness)
    sample_kurt <- apply(pre_samples, 2, kurtosis)
    moment_matrix <- cbind(sample_mean, sample_var, sample_skew, sample_kurt)
    para <- fit_johnsonsu_moment(moment_matrix)
    end_time = Sys.time()
    john_mom_samples <- generate_john(nits, para$xi, para$lambda, para$gamma, para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(john_mom_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished john_batch_mom")
    
    #6. Johnson su_batch_MLE
    method ="John_batch_MLE"
    start_time <- Sys.time()
    para <- fit_johnsonsu_MLE(pre_samples)
    end_time = Sys.time()
    john_mle_samples <- generate_john(nits, para$xi, para$lambda, para$gamma, para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(john_mle_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished john_batch_mle")
    
    #7. Johnson su_incremental_Mom
    method ="John_incremental_Mom"
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(john_online_mom_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished john_incremental_mom")
    
    #8. Johnson su_incremental_MLE
    method ="John_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_johnsonsu_MLE(log_pi,(nits -t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    john_online_mle_samples <- generate_john(nits, para$xi, para$lambda, para$gamma, para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(john_online_mle_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished john_incremental_mle")
    
    #9. SN_batch_Mom
    method ="SN_batch_Mom"
    start_time <- Sys.time()
    mu <- apply(pre_samples, 2, mean)
    variance <- apply(pre_samples, 2, var)
    skew <- apply(pre_samples, 2, skewness)
    para <- fit_sn_moment(mu, variance, skew)
    end_time = Sys.time()
    sn_mom_samples <- generate_sn(nits, para$xi, para$omega, para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(sn_mom_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sn_batch_mom")
    
    #10. SN_batch_MLE
    method ="SN_batch_MLE"
    start_time <- Sys.time()
    para <- fit_sn_MLE(pre_samples)
    end_time = Sys.time()
    sn_mle_samples <- generate_sn(nits, para$xi, para$omega, para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(sn_mle_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sn_batch_mle")
    
    #11. SN_incremental_Mom
    method ="SN_incremental_Mom"
    start_time <- Sys.time()
    para <- Transformation_sn_mom(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    sn_online_mom_samples <- generate_sn(nits, para$xi, para$omega, para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(sn_online_mom_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished sn_incremental_mom")
    
    #12. SN_incremental_MLE
    method ="SN_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_sn_MLE(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    sn_online_mle_samples <- generate_sn(nits, para$xi, para$omega, para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(
      wasserstein1d(f(target_samples), f(sn_online_mle_samples), p = 2), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished sn_incremental_mle")
  }
  
  max_wdist <- c()
  for(i in 1:length(epsilon_values)){
    max_d <- max(compare_res[i,paste0(method_names, "_Wdist")]%>%
                   select(-Sinh_incremental_MLE_Wdist, -John_incremental_Mom_Wdist))
    max_wdist[i] <- max_d
  }
  
  if(dividemax == TRUE){
    for(method in method_names){
      compare_res[,paste0(method, "_Wdist")] = compare_res[,paste0(method, "_Wdist")] /max_wdist
    }
  }
  
  # ---- 绘图 ----
  plot_data <- compare_res %>%
    select(epsilon, ends_with("_Wdist"), -Sinh_incremental_MLE_Wdist, -John_incremental_Mom_Wdist) %>%
    pivot_longer(cols = -epsilon, names_to = "method", values_to = "Wdist") %>%
    mutate(method = gsub("_Wdist", "", method))
  
  p <- ggplot(plot_data, aes(x = epsilon, y = Wdist, color = method, shape = method)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(linewidth = 1, alpha = 0.7) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Wasserstein Distance vs Epsilon",
      x = "Epsilon (Skewness Parameter)",
      y = "Wasserstein Distance",
      color = "Method",
      shape = "Method"
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    scale_x_continuous(breaks = epsilon_values) +
    scale_shape_manual(values = rep(c(16, 17, 15, 18, 19, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2), length.out = 15)) +
    guides(color = guide_legend(ncol = 3),
           shape = guide_legend(ncol = 3))
  
  return(list(compare_res = compare_res, p = p,max_wdist = max_wdist))
}

# Define a function to calculate the KL divergence between two distributions
KL_continuous <- function(logdens_p, logdens_q, lower=-10, upper=10, n=10000, eps = 1e-12) {
  x <- seq(lower, upper, length.out = n)
  
  # 对数密度转化为密度
  p <- exp(sapply(x, logdens_p))
  q <- exp(sapply(x, logdens_q))
  
  # 归一化
  dx <- x[2] - x[1]
  p <- p / sum(p * dx)
  q <- q / sum(q * dx)
  
  # 避免 log(0)
  q <- pmax(q, eps)
  
  # KL(P || Q) = ∫ p(x) log(p(x)/q(x)) dx
  kl <- trapz(x, p * log(p / q))
  
  return(kl)
}

# compare the total KL divergence of surrogate distribution gotten by different algorithm with the sinh-arcsinh target distribution
compare_KL_algorithms <- function(epsilon_values,delta_sinh=1, xi_sinh = 1, eta_sinh =1,
                                     nits = 5000, t0 = 2000) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(transport)
  
  # 方法名字（3个分布 × 4个拟合方式 = 12）
  method_names <- c("Sinh_batch_Mom", "Sinh_batch_MLE", "Sinh_incremental_Mom", "Sinh_incremental_MLE",
                    "John_batch_Mom", "John_batch_MLE", "John_incremental_Mom", "John_incremental_MLE",
                    "SN_batch_Mom", "SN_batch_MLE", "SN_incremental_Mom", "SN_incremental_MLE")
  
  # 保存结果
  compare_res <- matrix(0, nrow = length(epsilon_values), ncol = 1 + 2*length(method_names))
  colnames(compare_res) <- c("epsilon", paste0(method_names, "_Wdist"), paste0(method_names, "_time"))
  compare_res <- as.data.frame(compare_res)
  compare_res$epsilon <- epsilon_values
  
  # 循环 epsilon
  for (i in (1:length(epsilon_values))) {
    
    epsilon_sinh <- epsilon_values[i]
    cat("Current epsilon =", epsilon_sinh, "\n")
    
    target_samples = xi_sinh + eta_sinh * sinh((asinh(rnorm(nits)) + epsilon_sinh) / delta_sinh)
    min_x <- min(target_samples); max_x <- max(target_samples)
    log_pi <- function(x) log_dsinh_arcsinh(x, xi_sinh, epsilon_sinh, eta_sinh, delta_sinh)
    
    start_pre <- Sys.time()
    pre_sampling <- Adaptive_RWM(log_pi, nits, h, x_curr, target_a = 0.44)
    end_pre <- Sys.time()
    pre_samples <- as.matrix(pre_sampling$x_store[floor(pre_nits/3):pre_nits])
    pre_time <- as.numeric(end_pre - start_pre, units = "secs")
    
    #1. Sinh_batch_Mom
    method ="Sinh_batch_Mom"
    start_time <- Sys.time()
    para <- fit_highdim_sinh_arcsinh(pre_samples, method = "Moment")
    end_time = Sys.time()
    logd_sin_mom <- function(x) log_dsinh_arcsinh(x,para$xi,para$epsilon,para$eta,para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_sin_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sinh_batch_mom")
    
    #2. Sinh_batch_MLE
    method ="Sinh_batch_MLE"
    start_time <- Sys.time()
    para <- fit_highdim_sinh_arcsinh(pre_samples, method = "MLE")
    end_time = Sys.time()
    logd_sin_MLE<- function(x) log_dsinh_arcsinh(x,para$xi,para$epsilon,para$eta,para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_sin_MLE,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sinh_batch_mle")
    
    #3. Sinh_incremental_Mom
    method ="Sinh_incremental_Mom"
    start_time <- Sys.time()
    para <- Transformation_sinh_arcsinh_mom(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_sin_online_mom<- function(x) log_dsinh_arcsinh(x,para$xi,para$epsilon,para$eta,para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_sin_online_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished sinh_incremental_mom")
    
    #4. Sinh_incremental_MLE
    method ="Sinh_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_sinh_arcsinh_MLE(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_sin_online_MLE<- function(x) log_dsinh_arcsinh(x,para$xi,para$epsilon,para$eta,para$delta)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_sin_online_MLE,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished sinh_incremental_mle")
    
    #5. Johnson su_batch_Mom
    method ="John_batch_Mom"
    start_time <- Sys.time()
    sample_mean <- apply(pre_samples, 2, mean)
    sample_var <- apply(pre_samples, 2, var)
    sample_skew <- apply(pre_samples, 2, skewness)
    sample_kurt <- apply(pre_samples, 2, kurtosis)
    moment_matrix <- cbind(sample_mean, sample_var, sample_skew, sample_kurt)
    para <- fit_johnsonsu_moment(moment_matrix)
    end_time = Sys.time()
    logd_john_mom <- function(x) log_dJohnson(x,para$lambda,para$delta,para$xi,para$gamma)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_john_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished john_batch_mom")
    
    #6. Johnson su_batch_MLE
    method ="John_batch_MLE"
    start_time <- Sys.time()
    para <- fit_johnsonsu_MLE(pre_samples)
    end_time = Sys.time()
    logd_john_mle <- function(x) log_dJohnson(x,para$lambda,para$delta,para$xi,para$gamma)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_john_mle,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished john_batch_mle")
    
    #7. Johnson su_incremental_Mom
    method ="John_incremental_Mom"
    start_time <- Sys.time()
    para <- Transformation_johnsonsu_mom(log_pi,(nits -t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_john_online_mom <- function(x) log_dJohnson(x,para$lambda,para$delta,para$xi,para$gamma)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_john_online_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished john_incremental_mom")
    
    #8. Johnson su_incremental_MLE
    method ="John_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_johnsonsu_MLE(log_pi,(nits -t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_john_online_mle <- function(x) log_dJohnson(x,para$lambda,para$delta,para$xi,para$gamma)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_john_online_mle,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished john_incremental_mle")
    
    #9. SN_batch_Mom
    method ="SN_batch_Mom"
    start_time <- Sys.time()
    mu <- apply(pre_samples, 2, mean)
    variance <- apply(pre_samples, 2, var)
    skew <- apply(pre_samples, 2, skewness)
    para <- fit_sn_moment(mu, variance, skew)
    end_time = Sys.time()
    logd_sn_mom <- function(x) log_dsn(x,para$xi,para$omega,para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_sn_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sn_batch_mom")
    
    #10. SN_batch_MLE
    method ="SN_batch_MLE"
    start_time <- Sys.time()
    para <- fit_sn_MLE(pre_samples)
    end_time = Sys.time()
    logd_sn_mle <- function(x) log_dsn(x,para$xi,para$omega,para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_sn_mle,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time + pre_time, units = "secs")
    print("finished sn_batch_mle")
    
    #11. SN_incremental_Mom
    method ="SN_incremental_Mom"
    start_time <- Sys.time()
    para <- Transformation_sn_mom(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_sn_online_mom <- function(x) log_dsn(x,para$xi,para$omega,para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_sn_online_mom,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time, units = "secs")
    print("finished sn_incremental_mom")
    
    #12. SN_incremental_MLE
    method ="SN_incremental_MLE"
    start_time <- Sys.time()
    para <- Transformation_sn_MLE(log_pi, (nits - t0), x_curr = 0, t0 = t0)
    end_time = Sys.time()
    logd_sn_online_mle <- function(x) log_dsn(x,para$xi,para$omega,para$alpha)
    compare_res[i, paste0(method, "_Wdist")] <- round(KL_continuous(log_pi, logd_sn_online_mle,lower = min_x,upper = max_x), 3)
    compare_res[i, paste0(method, "_time")] <- as.numeric(end_time - start_time , units = "secs")
    print("finished sn_incremental_mle")
  }
  
  
  # ---- 绘图 ----
  plot_data <- compare_res %>%
    select(epsilon, ends_with("_Wdist")) %>%
    pivot_longer(cols = -epsilon, names_to = "method", values_to = "Wdist") %>%
    mutate(method = gsub("_Wdist", "", method))
  
  p <- ggplot(plot_data, aes(x = epsilon, y = Wdist, color = method, shape = method)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(linewidth = 1, alpha = 0.7) +
    theme_minimal(base_size = 12) +
    labs(
      title = "KL Divergence vs Epsilon",
      x = "Epsilon (Skewness Parameter)",
      y = "KL Divergence",
      color = "Method",
      shape = "Method"
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    scale_x_continuous(breaks = epsilon_values) +
    scale_shape_manual(values = rep(c(16, 17, 15, 18, 19, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2), length.out = 15)) +
    guides(color = guide_legend(ncol = 3),
           shape = guide_legend(ncol = 3))
  
  return(list(compare_res = compare_res, p = p))
}











