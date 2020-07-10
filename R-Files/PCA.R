#load needed packages
library(MASS)
library(glmnet)
#install.packages('pracma')
library(pracma);
#install.packages('matrixStats')
library(matrixStats);

PCA <- function(X, transform = 'normalize',  EV_scal = 10) {
  #input: - X: matrix for which the PC should be computed
  #       - EV_scal: how many digits must be the same such that the length of the eigenvalues is accepted as being one (needed because of computational issues of apply)
  #       - transform: if 1 the matrix of covariates is demeaned, if 2 the matrix of covariates is demeaned and the variances are scaled to zero, if 0 no transformation is done
  #output: A list containing
  #       - Z: matrix of Principal Components
  #       - phi: matrix of eigenvectors
  #       - eigenval: sorted eigenvectors
  #       - X: untransformed X
  #       - X_transformed: transformed X. Can either be untransformed, normalized or standardized
  #       - adjusted_EV: Boolean - TRUE if length of Eigenvectpors was adjusted FALSE if length was not adjusted
  
  X_untransformed <- X
  
  if ( transform == 'normalize' ) { #demean
    X <- scale(X, TRUE, FALSE)   
  }
  
  if ( transform == 'standardize' ) {
    X <- scale(X, TRUE, TRUE) 
  }
  
  adjusted_EV <- FALSE
  phi <- eigen(cov(X))$vectors #get eigenvectors from X'X. (already sorted by eigenvalues)
  eigenvalues <- eigen(cov(X))$values
  
  if (any( round(apply(phi, 2, function(col_phi){sum(col_phi^2)}),EV_scal) != 1))  { #check whether length of Eigenvectors is one
    (apply(phi, 2, function(col_phi){col_phi /( sum(col_phi^2)^0.5 )})) 
    adjusted_EV <- TRUE
  }
  
  Z = X %*% phi #compute Z matrix for M regressors
  list_return <- list('Z' = Z, 'phi' = phi, 'eigenval' = eigenvalues, 'X' = X_untransformed, 'X_transformed' = X, 'adjusted_EV' = adjusted_EV,)
  return(list_return)
}