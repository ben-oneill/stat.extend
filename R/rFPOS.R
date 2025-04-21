#' Random generation from the Finite Population Order Statistic (FPOS) distribution
#'
#' \code{rFPOS} returns a vector of pseudo-random values from the distribution.
#'
#' This function generates pseudo-random values from the cumulative distribution function of the Finite
#' Population Order Statistic (FPOS) distribution.  This is the distribution of the order statistics when
#' sampling without replacement from a finite population.  Further details on the distribution and its properties
#' can be found in the following paper:
#'
#' O'Neill, B. (2023) The distribution of order statistics under sampling without replacement.
#'
#' If the user specifies a population size N then the population is taken to be the values 1,...,N.  The user may
#' specify an alternative population of numeric values, which may include repeat values.  If the population is
#' specified in the input then the user need not specify the population size N.  (In the event that both parameters
#' are specified and are inconsistent, the input N is ignored and the function gives a warning.)
#'
#' @usage \code{rFPOS(r, population = NULL, k, n, N = length(population))}
#' @param r The number of pseudo-random values to generate
#' @param population Optional vector specifying the population values
#' @param k The rank-order (a positive integer no greater than the sample size)
#' @param n The sample size (a positive integer no greater than the population size)
#' @param N The population size (a positive integer)
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be a vector of r pseudo-random values from the FPOS distribution

rFPOS <- function(r, population = NULL, k, n, N = length(population)) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(r))                       stop('Error: Argument nn is not numeric')
  if (!is.numeric(N))                       stop('Error: Input parameter N is not numeric')
  if (!is.numeric(n))                       stop('Error: Input parameter n is not numeric')
  if (!is.numeric(k))                       stop('Error: Input parameter k is not numeric')

  #Check that parameters are atomic
  if (length(N) != 1)                       stop('Error: Input parameter N should be a single number')
  if (length(n) != 1)                       stop('Error: Input parameter n should be a single number')
  if (length(k) != 1)                       stop('Error: Input parameter k should be a single number')

  #Check that parameters are in allowable range
  NN <- as.integer(N)
  nn <- as.integer(n)
  kk <- as.integer(k)
  if (N != NN)                              stop('Error: Input parameter N not an integer')
  if (n != nn)                              stop('Error: Input parameter n not an integer')
  if (k != kk)                              stop('Error: Input parameter k not an integer')
  if (k <= 0)                               stop('Error: Input parameter k must be positive')
  if (n < k)                                stop('Error: Input parameter k cannot be larger than n')
  if (N < n)                                stop('Error: Input parameter n cannot be larger than N')

  #Check the input r
  rr <- as.integer(r)
  if (r != rr)                              stop('Error: Input parameter r not an integer')
  if (r < 0)                                stop('Error: Input parameter r must be non-negative')
  if (r == 0) { return(numeric(0)) }

  #Generate values
  UU   <- rbeta(r, shape1 = k, shape2 = n-k+1)
  RAND <- k + rbinom(r, size = N-n, prob = UU)

  #Give output
  RAND }
