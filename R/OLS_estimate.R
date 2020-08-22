#OLS simulation

var_Y_hat_OLS <- function(population_Y, sample_size, confidence, numb_it, type_CI = 'mean', line = 'mean', print_it = FALSE) {
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
      store_Y_vars <- cbind(store_Y_vars, colVars(store_Y_hat))
    }
    low <- (1-confidence)/2
    up <- (1-confidence)/2 + confidence
    mean <- colMeans(store_Y_vars)
    sds <- colVars(store_Y_vars)^0.5/numb_it^0.5
    if (type_CI == 'mean') {
      upper <- mean + qt(up, numb_it - 1)*sds
      lower <- mean - qt(up, numb_it - 1)*sds
    } else{
      upper <-apply(store_Y_vars, MARGIN = 2, FUN = quantile, probs = up)
      lower <- apply(store_Y_vars, MARGIN = 2, FUN = quantile, probs = low)
    }

    plot_line <- mean

    if (line != 'mean') {
      plot_line <- store_Y_vars[1,]
    }

    if (line == 'none') {
      plot_line <- rep(NaN, length(mean))
    }

    Y_var_CI <- cbind(upper, plot_line, lower)
  return(Y_var_CI)
}
