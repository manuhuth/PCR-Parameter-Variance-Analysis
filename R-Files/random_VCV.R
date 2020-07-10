#load needed packages
library(MASS)
library(glmnet)
#install.packages('pracma')
library(pracma);
#install.packages('matrixStats')
library(matrixStats);

frobenius_norm <- function(A) {
  #input: matrix for which the Frobenius norm/ Hilbert-Schmidt norm should be computed
  #output: real number yielding the value of the Frobenius Norm 
  H = t(A) %*% A
  h = Trace(H)^0.5
  return(h)
}

random_VCV <- function(lambda, low_bound_cov = 0.2, max_bound_cov = 10, low_bound_av = 0, max_it = 500, diagonal = NULL, pos = FALSE) {
  #input: - lambda: vector of eigenvalues for the VCV matrix. Length of Vector determines number of covariates
  #       - max/low_bound_cov: upper/lower bound for the absolute values of the covariances
  #       - max/low_bound_av: upper/lower bound for the average of the squared non-diagonal entries
  #       - max_it: maximum numbers of iterations
  #       - diagonal: specififes the elements on the diagonal (optional argument)
  #output: A list containing
  #       - VCV: VAriance-Covariance matrix 
  #       - lambda: chosen eigenvalues
  #       - Q: matrix of eigenvectors
  #       - it: number of iterations
  #dependencies: function frob_norm, pracma package
  p = length(lambda)
  D = diag(lambda)
  
  it = 0
  frob_av = 0
  low_cov = FALSE
  max_cov = FALSE
  pos_cov = FALSE
  VCV = matrix(0,p,p)
  
  while( (it < max_it) & ( (frob_av < low_bound_av) | (isFALSE(low_cov)) | (isFALSE(pos_cov))  | (isFALSE(max_cov)) ) ){ # (frob_av > max_bound_av) 
    Q = randortho(p) #generate random orthogonal matrix -> need to get more control. QR decomposition?
    VCV = Q %*% D %*% t(Q)
    if (is.null(diagonal) == FALSE) {
      H <- diag( (diagonal/diag(VCV))^0.5 )
      VCV <- H %*% VCV %*% H #multiply each cell by the factors to obtain the right values on the diagonal 
    }
    
    frob_av = ( frobenius_norm(VCV)^2 - sum(diag(VCV)^2) )  / (p*(p-1)) #exclude diagonal elements since they are covariances
    low_cov = sum(abs(VCV) >= low_bound_cov) == p*p #check whether all entries are bigger than the specified lower bound 'low_bound_cov'
    max_cov = sum(abs(VCV) <= max_bound_cov) - sum(abs(diag(VCV)) <= max_bound_cov) == p*(p-1) #check whether all entries are smaller than the specified upper bound 'max_bound_cov'
    pos_cov = sum(VCV >= 0) == p*p 
    if (isFALSE(pos)) {
      pos_cov = TRUE 
    }
    it = it + 1
  }  
  if  ( (frob_av < low_bound_av) | (isFALSE(low_cov)) | (isFALSE(pos_cov)) | (isFALSE(max_cov))) { #check the matrix with which we exit; if it does not fulfill conditions, we have not found a solution (maybe easier using max_it)
    warning('Could not find VCV that satisfies the requirements')
  }
  
  list_return = list('VCV' =VCV, 'lambda' = lambda, 'Q' = Q, 'it' = it)
  return(list_return)
}  