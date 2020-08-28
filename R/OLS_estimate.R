#OLS simulation

var_Y_hat_OLS <- function(population_Y, sample_size, confidence, numb_it, type_CI = 'mean', line = 'mean', print_it = FALSE) {
  #input: - population_Y: matrix containing the dependent variable as first column and the regressors in the subsequent columns
  #       - sample_size: vector containing the sample sizes
  #       - confidence: confidence level in which the variance estimate should be
  #       - type_CI/line: if 'mean' the confidence intervals/plotted line of the mean are plotted. If != 'mean' quantiles are plotted
  #       - numb_it: umber of iterations used to compute the means of the standard errors
  #       - print_it: if TRUE the algorithm shows the number of iterations while executing
  #output: - data frame containing the lower/upper bounds and the means of the variances
    store_Y_vars <- c()
    store_Y_vars_PCR <- c()
    for (index_sample in sample_size) {
      store_Y_hat <- c()
      store_beta <- c()
      store_Y_PCR <- c()
      for (index in 1:numb_it) {
        sample <- as.matrix(population_Y[sample(nrow(population_Y), sample_size), ])
        Y <- sample[,1]
        X <- sample[,-1]
        beta <- solve(t(X) %*% X)%*%t(X)%*%Y
        Y_hat <- X %*% beta
        store_beta <- rbind(store_beta, t(beta))
        store_Y_hat <- cbind(store_Y_hat, Y_hat)

        Y_hat_PCR <- PCR(Y = Y, X = X, M = ncol(X))$Y_hat
        store_Y_PCR <- cbind(store_Y_PCR, Y_hat_PCR)
        if (isTRUE(print_it)) {
          print(index_sample)
          print(index)
        }
      }
      store_Y_vars <- cbind(store_Y_vars, colVars(store_Y_hat))
      store_Y_vars_PCR <- cbind(store_Y_vars_PCR, colVars(store_Y_PCR))
    }
    low <- (1-confidence)/2
    up <- (1-confidence)/2 + confidence
    mean <- colMeans(store_Y_vars)
    sds <- colVars(store_Y_vars)^0.5/numb_it^0.5
    mean_PCR <- colMeans(store_Y_vars_PCR)
    sds_PCR <- colVars(store_Y_vars_PCR)^0.5/numb_it^0.5
    if (type_CI == 'mean') {
      upper <- mean + qt(up, numb_it - 1)*sds
      lower <- mean - qt(up, numb_it - 1)*sds
      upper_PCR <- mean + qt(up, numb_it - 1)*sds_PCR
      lower_PCR <- mean - qt(up, numb_it - 1)*sds_PCR
    } else{
      upper <-apply(store_Y_vars, MARGIN = 2, FUN = quantile, probs = up)
      lower <- apply(store_Y_vars, MARGIN = 2, FUN = quantile, probs = low)
      upper_PCR <-apply(store_Y_vars_PCR, MARGIN = 2, FUN = quantile, probs = up)
      lower_PCR <- apply(store_Y_vars_PCR, MARGIN = 2, FUN = quantile, probs = low)
    }

    plot_line <- mean
    plot_line_PCR <- mean_PCR
    if (line != 'mean') {
      plot_line <- store_Y_vars[1,]
      plot_line_PCR <-  store_Y_vars_PCR[1,]
    }

    if (line == 'none') {
      plot_line <- rep(NaN, length(mean))
      plot_line_PCR <- rep(NaN, length(mean))
    }

    Y_var_CI <- cbind(upper, plot_line, lower)
    Y_var_CI_PCR <- cbind(upper_PCR, plot_line_PCR, lower_PCR)
    list_return <- list('OLS' = Y_var_CI, 'PCR' = Y_var_CI_PCR)
  return(list_return)
}
