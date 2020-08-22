#OLS simulation

var_Y_hat_OLS <- function(population_Y, sample_size, confidence, numb_it, print_it = FALSE) {
    store_Y_vars <- c()
    for (index_sample in sample_size) {
      store_Y_hat <- c()
      store_beta <- c()
      for (index in 1:numb_it) {
        sample <- as.matrix(population_Y[sample(nrow(population_Y), sample_size), ])
        Y <- sample[,1]
        X <- sample[,-1]
        beta <- solve(t(X) %*% X)%*%t(X)%*%Y
        Y_hat <- X %*% beta
        store_beta <- rbind(store_beta, t(beta))
        store_Y_hat <- cbind(store_Y_hat, Y_hat)
        if (isTRUE(print_it)) {
          print(index_sample)
          print(index)
        }
      }
      store_Y_vars <- cbind(store_Y_vars, colVars(store_Y_hat)) #variances from 50, last column variances from 5000
    }
    low <- (1-confidence)/2
    up <- (1-confidence)/2 + confidence
    mean <- colMeans(store_Y_vars)
    sds <- colVars(store_Y_vars)^0.5/numb_it^0.5
    upper <- mean + qnorm(up, mean, sds)*sds
    lower <- mean - qnorm(up, mean, sds)*sds
    Y_var_CI <- cbind(upper, mean, lower)
  return(Y_var_CI)
}

