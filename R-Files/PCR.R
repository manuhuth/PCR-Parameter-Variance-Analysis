PCR <- function(Y, X, M, transform_Y = 'normalized', transform_X = 'normalized', EV_scal = 10, intercept = FALSE) {
  #input: - Y: dependent variable
  #       - X: covariates matrix
  #       - M: number of principal components
  #       - transform_Y/X: if 'normalize the matrix of covariates/vector of the independent variable is demeaned, if 'standardized the matrix of covariates /vector of the independent variable is 
  #                       demeaned and the variances are scaled to one if 'no' no transformation is performed
  #       - EV_scal: how many digits must be the same such that the length of the eigenvalues is accepted as being one (computational issue), passed to PCA
  
  #output: A list containing
  #       - Y_hat: fitted values
  #       - resid: residuals
  #       - mse: mean squarred error
  #       - beta_Z: estimated beta
  #       - Z: Pricnipal components
  #       - phi: matrix of eigenvectors
  #       - Y: untransformed dependent variable
  #       - Y_transformed: transformed dependent variable
  #       - X: untransformed independent variables
  #       - X_transformed: transformed independent variablea
  #       - transform_Y/X, intercept, M: see input
  
  #dependencies: PCA function
  
  if ( ( isFALSE( (transform_Y %in% c(0,1,2)) ) ) | ( isFALSE( (transform_X %in% c(0,1,2)) ) ) ){
    stop('transform_X and transform_Y must either be 0, 1 or 2')
  } 
  
  #save untransformed
    X_untransformed <- X
    Y_untransformed <- Y
  
  #transform Y and X
    if (transform_Y == 'normalize') { #demean Y
      Y <- scale(Y, center = TRUE, scale = FALSE)
    }
    if (transform_Y == 'standardize') { #standardize Y
      Y <- scale(Y, TRUE, TRUE)              
    }
    

  #Do PCA
    pca <- PCA(X = X, transform = transform_X, EV_scal = EV_scal) 
    Z <- pca$Z
    Z_without_intercept <- Z 
    if (intercept == TRUE) {
      Z <- cbind(rep (1,nrow(Z)), Z)
    } 
  
  #compute PC beta and store results  
    beta_Z <- solve(t(Z) %*% Z) %*% t(Z) %*% Y
  
  #compute Y hat
    Y_hat <- Z[,1:M] %*% beta_Z[1:M]
  #compute residuals and MSE
    residuals <- Y - Y_hat
    mse <- mean(residuals^2)
  
  list_return <- list('Y_hat' = Y_hat, 'resid' = residuals, 'mse' = mse, 'beta_Z' = beta_Z, 'Z' = Z, 'phi' = pca$phi, 'Y' = Y_untransformed, 'Y_transformed' = Y,
                      'X' = X_untransformed, 'X_transformed' = pca$X_transformed, 'transform_X' = transform_X, 'transform_Y' = transform_Y, 'intercept' = intercept, 'M' = M) 
  
  class(list_return) <- 'PCR' 
  return(list_return)
} 