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
  cum_dens <- c(0,cumsum(prob_mass))
  positions <- sapply(random, rdiscr_position, cum_dens = cum_dens)
  random_value <- X[positions]

  return(random_value)
}


dgPois <- function(x, lambda, eta) {

  prob <- lambda*(lambda + x*eta)^(x-1)*exp(-lambda - x*eta)/factorial(x)
  return(prob)
}

rgenPois <- function(n, lambda, eta, length_cumDensity = 100){
  #input: - n: number of random values
  #       - lambda: first parameter that characterizes mean and variance
  #       - eta: second parameter that characterizes mean and variance.
  #output:- random_value: vector of generalized Poisson distributed random variables

  distribution_value <- runif(n)
  cum_dens <- cumsum(dgPois(0:length_cumDensity, lambda = lambda, eta = eta)) #omega logic is vice versa
  positions <- sapply(distribution_value, rdiscr_position, cum_dens = cum_dens)
  numb_help <- 0:length_cumDensity
  random_value <- numb_help[positions]
  output <- list('random_value' = random_value, 'distribution_value' = distribution_value)
  return(output)
}



rbinom_wrapper <- function(x, p){
  #input: - see for a full overview: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Binomial.html
  # wrapper to use sapply
  output <- rbinom(1, x, p)
  return(output)
}


