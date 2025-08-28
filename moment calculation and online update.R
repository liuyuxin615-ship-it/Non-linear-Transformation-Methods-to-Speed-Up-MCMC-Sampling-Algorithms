library(moments)
# Define a function to online update the moments
# 辅助函数：从样本计算中心矩
compute_central_moments <- function(sample) {
  n <- length(sample)
  mean_x <- mean(sample)
  
  # 中心化数据
  centered <- sample - mean_x
  
  # 计算中心矩
  m2 <- sum(centered^2) / n  # 二阶中心矩
  m3 <- sum(centered^3) / n  # 三阶中心矩
  m4 <- sum(centered^4) / n  # 四阶中心矩
  
  return(list(
    x_bar = mean_x,
    m2 = m2,
    m3 = m3,
    m4 = m4
  ))
}

update_moments <- function(x_new,n,mean_prev,variance_prev,skew_prev,kurtosis_prev=0){
  
  m2 <- as.matrix(variance_prev * (n - 1))  # second central moment（sum((x_i - mean)^2)）
  m3 <- skew_prev * (diag(m2))^(1.5) / n^0.5  # third central moment（sum((x_i - mean)^3)）
  m4 <- (kurtosis_prev) * diag(m2)^2 / n  # fourth central moment（sum((x_i - mean)^4)）
  #print(m4)
  
  n  = n + 1
  delta = x_new - mean_prev
  delta_n = delta / n
  
  
  #update the cumulative central moment
  mu = mean_prev + delta_n
  m2 = as.matrix(m2 + delta %*% t( delta - delta_n ))
  m3 = m3 - 3.0 * delta_n * diag(m2) + delta * ( delta^2 - delta_n^2 )
  m4 = m4 - 4.0 * delta_n * m3 - 6.0 * delta_n^2 * diag(m2)  + delta * ( delta^3 - delta_n ^3)
  #print(m4)
  
  # calculate the central moment
  variance = m2 / (n-1)
  skew = m3 / diag(m2)^1.5 * n^0.5
  kurt = m4 / diag(m2)^2 * n
  
  # in case return NA or Inf
  mu = ifelse(is.finite(mu),mu,mean_prev)
  variance = ifelse(is.finite(variance),variance,variance_prev)
  skew = ifelse(is.finite(skew),skew,skew_prev)
  kurt = ifelse(is.finite(kurt),kurt,kurtosis_prev)
  
  
  return(list(
    mean = mu,
    variance = variance,
    skew = skew,
    kurt = kurt
  )
  )
}

test_update_moments <- function() {
  # 生成测试数据
  set.seed(123)
  d = 3
  data <- rmvnorm(100,mean = rep(0,d))
  
  # 初始化
  n <- 1
  mean_val <- data[1,]
  variance_val <- diag(rep(0,d),ncol = d)
  skew_val <- rep(0,d)
  kurtosis_val <- rep(0,d)
  
  # 逐个添加数据点
  for (i in 2:nrow(data)) {
    
    result <- update_moments(data[i,], n, mean_val, variance_val, skew_val, kurtosis_val)
    n <- n + 1
    mean_val <- result$mean
    variance_val <- result$variance
    skew_val <- result$skew
    kurtosis_val <- result$kurt
  }
  
  # 与直接计算对比
  library(moments)
  cat("Online update result:\n")
  cat("Mean:", mean_val, "\n")
  cat("Variance:", variance_val, "\n")
  cat("Skewness:", skew_val, "\n")
  cat("Kurtosis:", kurtosis_val, "\n\n")
  
  cat("Direct Calculation:\n")
  cat("Mean:", mean(data), "\n")
  cat("Variance:", var(data), "\n")
  cat("Skeness:", skewness(data), "\n")
  cat("Kurtosis:", kurtosis(data), "\n")  # moments包的kurtosis是总峰度，需要减3得到超额峰度
}

# run the test
# test_update_moments()