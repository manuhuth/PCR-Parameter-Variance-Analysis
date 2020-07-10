PCR_predict <- function(pcr, X_test, Y_test = NULL) {
  #inpute: - PC: an object created by the PCR function
  #        - X_test: matrix of Covariates for which the values should be predicted
  #        - Y_test: Vector of the true observations for the independent variable of test data; (optional); prediction errors are computed
  #output: A list containing
  #        - Y_pred: predicted Y
  #        - Y_pred_transformed: transformed predicted Y
  #        - MSPE: Mean Squarred Prediction Error, if Y_test is specified
  #        - pred_error: prediction error, if Y_test is specified
  
  if ( isFALSE(is(pcr, 'PC_reg')) ){
    stop('The input `pcr` must be an object of class PCR, which is created by the PCR function')
  } 
  #get values from pcr
    M <- pcr$M
    Y <- pcr$Y
    X <- pcr$X
    Y_transformed <- pcr$Y_transformed
    X_transformed <- pcr$X_transformed
    beta_Z <- pcr$beta_Z[1:M]
    phi <- pcr$phi[ , (1:M)]
    X_pred <- X_test
    
  #transform X if necessary
    if ( (pcr$transform_X == 'normalize') | (pcr$transform_X == 'standardize') ) {
      X_pred <- t(t(X_test) - colMeans(X))
      if (pcr$transform_X == 'standardize') {
        X_pred <- X_pred %*% diag(1/colVars(X)^0.5)
      }
    }

  
  #adjust X_pred if only one PC is used (R issue)
      if ( isFALSE( is.matrix(X_test) ) ) { 
        X_pred <- t(matrix(X_pred)) #adjust vector to matrix with right dimensions, if only one PC is used
      }
  #compute Z_pred and add intercept if desired
    Z_pred <- X_pred %*% phi
  
    if (isTRUE(PC$intercept)) {
      Z_pred <- cbind(rep(1,ncol(Z_pred)), Z_pred)
    }
  
  Y_pred <- Z_pred %*% beta_Z
  
  pred_error <- NULL
  MSPE <- NULL
  
  #compute prediction error and mean squarred average predictione error if true Y is specified
    if ( isFALSE( is.null(Y_test) ) ){
      Y_true <- Y_test
      #transform Y if necessary
      if ( (pcr$transform_Y == 'normalize') | (pcr$transform_Y == 'standardize') ) {
        Y_true <- Y_test - mean(Y)
        if (pcr$transform_X == 'standardize') {
          Y_true <- Y_true / var(Y)
        }
      }
    pred_error <- Y_true - Y_pred
    MSPE <- mean(pred_error^2)
  }
  #transform Y back to get predictions
    Y_pred_untransformed = Y_pred * var(Y + mean(Y))
  
  list_return = list('Y_pred' = Y_pred_untransformed , 'Y_pred_transformed' = Y_pred, 'MSPE' = MSPE, 'pred_error' = pred_error)
  return(list_return)
}