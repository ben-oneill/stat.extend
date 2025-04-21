#' Mass function of the Finite Population Order Statistic (FPOS) distribution
#'
#' \code{dFPOS} returns the probability or log-probability values for the arguments.
#'
#' This function computes probabilities or log-probabilities from the mass function of the Finite Population
#' Order Statistic (FPOS) distribution.  This is the distribution of the order statistics when sampling without
#' replacement from a finite population.  Further details on the distribution and its properties can be found 
#' in the following paper:
#'
#' O'Neill, B. (2023) The distribution of order statistics under sampling without replacement.
#'
#' If the user specifies a population size N then the population is taken to be the values 1,...,N.  The user may
#' specify an alternative population of numeric values, which may include repeat values.  If the population is 
#' specified in the input then the user need not specify the population size N.  (In the event that both parameters
#' are specified and are inconsistent, the input N is ignored and the function gives a warning.)
#'
#' @usage \code{dFPOS(x, population = NULL, k, n, N = length(population), log = FALSE)}
#' @param x A vector of numeric values to be used as arguments for the mass function
#' @param population Optional vector specifying the population values
#' @param k The rank-order (a positive integer no greater than the sample size)
#' @param n The sample size (a positive integer no greater than the population size)
#' @param N The population size (a positive integer)
#' @param log A logical value specifying whether results should be returned as log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be a vector of probabilities/log-probabilities corresponding to the vector argument x

dFPOS <- function(x, population = NULL, k, n, N = length(population), log = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(x))                       stop('Error: Argument x is not numeric')
  if (!is.numeric(N))                       stop('Error: Input parameter N is not numeric')
  if (!is.numeric(n))                       stop('Error: Input parameter n is not numeric')
  if (!is.numeric(k))                       stop('Error: Input parameter k is not numeric')
  if (!is.logical(log))                     stop('Error: log option is not a logical value')

  #Check that parameters are atomic
  if (length(N) != 1)                       stop('Error: Input parameter N should be a single number')
  if (length(n) != 1)                       stop('Error: Input parameter n should be a single number')
  if (length(k) != 1)                       stop('Error: Input parameter k should be a single number')
  if (length(log) != 1)                     stop('Error: log option should be a single logical value')

  #Check that parameters are in allowable range
  NN <- as.integer(N)
  nn <- as.integer(n)
  kk <- as.integer(k)
  if (N != NN)                              stop('Error: Input parameter N not an integer')
  if (n != nn)                              stop('Error: Input parameter n not an integer')
  if (k != kk)                              stop('Error: Input parameter k not an integer')
  if (k < 1)                                stop('Error: Input parameter k must be positive')
  if (n < k)                                stop('Error: Input parameter k cannot be larger than n')
  if (N < n)                                stop('Error: Input parameter n cannot be larger than N')
  
  #Set ordered population
  if (is.null(population)) { 
    POP <- 1:N
    POP.UNIQUE <- 1:N
  } else { 
    POP <- sort(population)
    POP.UNIQUE <- unique(POP)
    if (length(population) != N) {
      warning('Input N not used since it is inconsistent with input population') } }
  
  #Create log-probabilities for elements 1,...,N
  LOGFPOS <- rep(-Inf, N)
  for (i in 0:(N-n)) {
    LOGFPOS[k+i] <- lchoose(k+i-1, k-1) + lchoose(N-k-i, n-k) - lchoose(N,n) }
  
  #Create log-probabilities for population values
  NNpop <- length(POP.UNIQUE)
  LOGPROBS <- rep(-Inf, NNpop)
  for (i in 1:NNpop) {
    VALS <- which(POP == POP.UNIQUE[i])
    LOGPROBS[i] <- matrixStats::logSumExp(LOGFPOS[VALS]) }
  LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS)
  
  #Create output vector
  nx <- length(x)
  LOGOUT <- rep(-Inf, nx)
  for (i in 1:nx) {
    xx <- x[i]
    if (xx %in% POP.UNIQUE) { 
      VAL <- which(POP.UNIQUE == xx)
      LOGOUT[i] <- LOGPROBS[VAL] } }

  #Return output
  if (log) { LOGOUT } else { exp(LOGOUT) } }
