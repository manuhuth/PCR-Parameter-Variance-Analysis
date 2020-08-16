coef_var_analysis <- function(population, true_phi, M, sample_size, transform, iterations) {
  #input: - population: matrix containing the y variable as first column and the regressors as the other columns
  #       - true_phi: maatrix of the true eigenvectors
  #       - M: number of principal components that should be used
  #       - sample_size: size of the samples used for the analysis
  #       - transform: either 'no', 'normalized' or 'standardized'
  #       - iterations: number of iterations per estimate
  #output: - variances_beta_prac/theo: variances of the simulated beta in the stochastic/non-stochastic case using the empirical distribution
  #        - variances_beta_prac/theo_formula: variances of the simulated beta in the stochastic/non-stochastic case using formula 3.21
  #        - variances_Y_prac: Variances of Y hat in the stochastic case
  #        - variances_Y_theo: Variances of Y hat in the non-stochastic case
  #        - diff_Y_hat_changed: sum of differences in Y hat between stochastic and non-stochastic case
  #        - diff_beta_hat_changed: sum of differences in beta hat between stochastic and non-stochastic case

  store_beta_theo <- c()
  store_beta_prac <- c()
  store_Y_hat_theo <- c()
  store_Y_hat_prac <- c()
  store_Y_hat_prac_changed <- c()
  store_beta_prac_changed <- c()
  store_variance_beta_prac_formula <- c()
  store_variance_beta_prac_formula2 <- c()
  store_variance_beta_theo_formula <- c()

  for (k in 1:iterations) {
    #practical case with estimated phi at each iteration empirical variance
    sample_prac <- population[sample(nrow(population), sample_size), ]
    Y_sample_prac <- sample_prac[ , 1]
    X_sample_prac <- sample_prac[ , c(-1)]
    PCR_sample_prac <- PCR(Y = Y_sample_prac, X = X_sample_prac, M = M, transform_Y = trans, transform_X = trans)
    beta_prac <- PCR_sample_prac$beta_Z
    Y_hat_prac <- PCR_sample_prac$Y_hat

    #check if Y_hat changes with other representant of [\hat{phi}]
    phi_prac_changed <- PCR_sample_prac$phi %*% diag(c(-1, rep( 1, ncol(PCR_sample_prac$phi)-1) ) ) #first column times minus one
    Z_prac_changed <- PCR_sample_prac$X_transformed %*% phi_prac_changed
    beta_prac_changed <- solve(t(Z_prac_changed) %*% Z_prac_changed ) %*% t(Z_prac_changed) %*% Y_hat_prac
    Y_hat_prac_changed <- Z_prac_changed %*% beta_prac_changed

    #fixed true phi empirical variance
    sample_theo <- population[sample(nrow(population), l), ]
    if (transform == 'normalized') {
      X_Sample_theo_demeaned <- scale(sample_theo[,c(-1)], TRUE, FALSE)
      Y_sample_theo_demeaned <- scale(sample_theo[,1], TRUE, FALSE)
    }

    if (transform == 'standardized') {
      X_Sample_theo_demeaned <- scale(sample_theo[,c(-1)], TRUE, TRUE)
      Y_sample_theo_demeaned <- scale(sample_theo[,1], TRUE, TRUE)
    }

    Z_sample_theo <- X_Sample_theo_demeaned %*% true_phi
    beta_theo <- solve(t(Z_sample_theo) %*% Z_sample_theo) %*% t(Z_sample_theo) %*% Y_sample_theo_demeaned
    Y_hat_theo <- Z_sample_theo %*% beta_theo

    #practical case with estimated phi at each iteration formula 2.21
    eigenvalues_hat_prac <-  PCR_sample_prac$eigenval
    error_estimate_prac <- PCR_sample_prac$Y_transformed - Y_hat_prac
    variance_error_estimate_prac <- sum(error_estimate_prac^2) / (l - ncol(X_sample_prac)) #adjust by degrees of freedom
    variance_beta_prac_formula <- 1/eigenvalues_hat_prac * variance_error_estimate_prac / (sample_size-1) #adjust in code because
    # the function PCR returns the eigenvalues from X'X/(n-1)
    variance_beta_prac_formula2 <- variance_error_estimate_prac * diag(solve(t(PCR_sample_prac$Z) %*% PCR_sample_prac$Z)) #check if computation is alright


    #fixed true phi formula 2.16
    error_estimate_theo <- Y_sample_theo_demeaned - Y_hat_theo
    variance_error_estimate_theo <- sum(error_estimate_theo^2) / ( (l - ncol(X_sample_prac)))
    variance_beta_theo_formula <- variance_error_estimate_theo * diag(solve(t(Z_sample_theo) %*% Z_sample_theo))

    #store results
    store_Y_hat_prac <- cbind(store_Y_hat_prac, Y_hat_prac) # ys in columns for each sample -> colVars works
    store_Y_hat_theo <- cbind(store_Y_hat_theo, Y_hat_theo)
    store_Y_hat_prac_changed <- cbind(store_Y_hat_prac_changed, Y_hat_prac_changed)

    store_beta_prac <- rbind(store_beta_prac, t(beta_prac)) # betas in columns -> colVars works
    store_beta_theo <- rbind(store_beta_theo, t(beta_theo))
    store_beta_prac_changed <- rbind(store_beta_prac_changed, t(beta_prac_changed))

    store_variance_beta_prac_formula <- rbind(store_variance_beta_prac_formula, t(variance_beta_prac_formula))
    store_variance_beta_prac_formula2 <- rbind(store_variance_beta_prac_formula2, t(variance_beta_prac_formula2))
    store_variance_beta_theo_formula <- rbind(store_variance_beta_theo_formula, t(variance_beta_theo_formula))
  }

  #evaluate results
  variances_beta_prac <- colVars(store_beta_prac)
  variances_beta_prac_formula <- colMeans(store_variance_beta_prac_formula)
  variances_beta_prac_formula2 <- colMeans(store_variance_beta_prac_formula2)

  variances_beta_theo <- colVars(store_beta_theo)
  variances_beta_theo_formula <- colMeans(store_variance_beta_theo_formula)

  variances_Y_prac <- colVars(store_Y_hat_prac)
  variances_Y_theo <- colVars(store_Y_hat_theo)

  diff_Y_hat_changed <- colSums(store_Y_hat_prac_changed - store_Y_hat_prac) #means from every iteration
  diff_beta_hat_changed <- colSums(store_beta_prac_changed - store_beta_prac)

  list_return <- list('variances_beta_prac' = variances_beta_prac, 'variances_beta_prac_formula' = variances_beta_prac_formula,
                      'variances_beta_theo' =  variances_beta_theo, ' variances_beta_theo_formula' =  variances_beta_theo_formula,
                      'variances_Y_prac' = variances_Y_prac, 'variances_Y_theo' = variances_Y_theo, 'diff_Y_hat_changed' = diff_Y_hat_changed,
                      'diff_beta_hat_changed ' = diff_beta_hat_changed )
  return(list_return)
}
