rdiscr_position <- function(y, cum_dens) {
  #input: - y is a vector containing the functional values of the probability density function of a discrete random variable
  #       - cum_dens is a vector with the cumulative density of the random variable
  #output:- the place of the first number that is greater than the respective value in y

  return(Position(function(x) x >= y, cum_dens))
}

rdiscrete <- function(n, prob_mass, X) {
  #input: - n: number of random values
  #       - prob_mass is a vector containing the probability mass of a discrete random variable
  #       - X: values of the random variable X in an ascending order.
  #output:- n random numbers that are distributed according to the given CDF
  #dependencies: rdiscr_position function

  random <- runif(n)
  cum_dens <- cumsum(prob_mass)
  positions <- sapply(random, rdiscr_position, cum_dens = cum_dens)
  random_value <- X[positions]

  return(random_value)
}


library(LaplacesDemon)

rgenPois <- function(n, lambda, eta){
  random_numbers <- numeric(n)

  for (j in 1:n) {
    index <- 0
    unif <- runif(1)
    nain <- dgpois(0, lambda = lambda, omega = eta)
    while(nain < unif) {
      index <- index + 1
      nain <- nain + dgpois(index, lambda = lambda, omega = eta)
    }
    random_numbers[j] <- index
  }
  output <- list('rnumb' = random_numbers, 'cumDensity')
  return(output)
}



