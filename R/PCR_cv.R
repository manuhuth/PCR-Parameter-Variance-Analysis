PCR_cv <- function(Y, X, K, random = FALSE, transform_Y = 'normalized', transform_X = 'normalized', intercept = FALSE) {
  #input: - Y: input vector for cross validation
  #       - X: input matrix for cross validation
  #       - K: K-fold cross validation
  #       - random: if TRUE the rows of Y and X are chosen randomly
  #       - transform_Y/X: passed to the PCR function
  #       - intercept: passed to the PCR function
  
  #output: cross validated errors for every number of Principal Components
  #dependencies: PCR function, PCR_predict function
  
  if ( ( isFALSE( (transform_Y %in% c('no','normalized','standardized')) ) ) | ( isFALSE( (transform_X %in% c('no','normalized','standardized')) ) ) ){
    stop('transform_X and transform_Y must either be "no","normalized","standardized" ')
  } 
  
  X <- as.matrix(X)
  M <- ncol(X)
  n <- length(Y)
  sample_n <- ceil(n / K)
  data <- as.data.frame(cbind(Y,X))
  
  if (isTRUE(random)) { #randomizes the order of the rows
    names <- c()
    for(k in 1:ncol(X) ) {
      names[k] <- paste('X',k, '')
    }
    colnames(data) <- c('Y', names)
    
    data <- data[sample(nrow(data)), ] #randomize position of data
  }
  
  samples <- list() #split data in K samples, if the number of obs is not divisble by the number of samples, the last sample has fewer entries (should not be a problem)
  for (j in 1:K) {
    index_low <- ((j-1)*sample_n + 1)
    index_high <- (j*sample_n)
    if (index_high > n) {
      index_high = n
    }
    samples[[j]] <- data[index_low:index_high,]             
  }
  
  index <- 1:K
  M_cve <- c()
  for (m in 1:M) {
    pe <- c()
    for (j in 1:K) { #compute average prediction error
      train_index <- index[-j] #delete index for test data
      sample_train <- c()
      for (i in train_index) {
        sample_train <- rbind(sample_train, samples[[i]]) #create matrix of train data
      }
      sample_test <- samples[[j]] #create matrx of test data
      fitted <- PCR(Y = as.matrix(Y = sample_train[,1]), X = as.matrix(sample_train[,c(-1)]), M = m, transform_Y = transform_Y, transform_X = transform_X, intercept = intercept) #get fitted object for prediction
      pe[j] <- sum((PCR_predict(pcr = fitted, X_test = as.matrix(sample_test[,c(-1)]), Y_test = sample_test[,1])$pred_error)^2) # compute sum of squared predicion errors for the sample
    }
    cve <- sum(pe)/n #average sum over all prediction errors
    M_cve[m] <- cve
  }
  cve <- cbind(1:M, M_cve)
  colnames(cve) <- c('M', 'cve')
  return(cve)
}