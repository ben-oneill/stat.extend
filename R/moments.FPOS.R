#' Moments of the Finite Population Order Statistic (FPOS) distribution
#'
#' \code{moments.FPOS} returns some representative moments from the distribution.
#'
#' This function computes some representative moments from the Finite Population Order Statistic (FPOS)
#' distribution.  This is the distribution of the order statistics when sampling without replacement from
#' a finite population.  Further details on the distribution and its properties can be found in the following
#' paper:
#'
#' O'Neill, B. (2023) The distribution of order statistics under sampling without replacement.
#'
#' If the user specifies a population size N then the population is taken to be the values 1,...,N.  The user may
#' specify an alternative population of numeric values, which may include repeat values.  If the population is
#' specified in the input then the user need not specify the population size N.  (In the event that both parameters
#' are specified and are inconsistent, the input N is ignored and the function gives a warning.)
#'
#' @usage \code{moments.FPOS(population = NULL, k, n, N = length(population), include.sd = FALS)}
#' @param population Optional vector specifying the population values
#' @param k The rank-order (a positive integer no greater than the sample size)
#' @param n The sample size (a positive integer no greater than the population size)
#' @param N The population size (a positive integer)
#' @include.sd Logical value; if \code{TRUE} the output includes the standard deviation
#' @return  If all inputs are correctly specified (i.e., parameters are in allowable range) then the output
#' will be a data frame of moments

moments.FPOS <- function(population = NULL, k, n, N = length(population), include.sd = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(N))                       stop('Error: Input parameter N is not numeric')
  if (!is.numeric(n))                       stop('Error: Input parameter n is not numeric')
  if (!is.numeric(k))                       stop('Error: Input parameter k is not numeric')
  if (!is.logical(include.sd))              stop('Error: Input include.sd is not a logical value')

  #Check that parameters are atomic
  if (length(N) != 1)                       stop('Error: Input parameter N should be a single number')
  if (length(n) != 1)                       stop('Error: Input parameter n should be a single number')
  if (length(k) != 1)                       stop('Error: Input parameter k should be a single number')
  if (length(include.sd) != 1)              stop('Error: Input include.sd should be a single logical value')

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

  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################

  #Compute moments
  if (is.null(population)) {

    #Set moments
    MEAN <- k*(N+1)/(n+1)
    VAR  <- k*(n-k+1)*(N+1)*(N-n)/(((n+1)^2)*(n+2))
    SKEW <- (n-2*k+1)*(1+2*(N-n-1)/(n+3))*sqrt((n+2)/((N+1)*(N-n)*k*(n-k+1)))
    KURT <- 3 + (n*(n+1)^3*(n+2))/((N+1)*(N-n)*k*(n-k+1)*(n+3)*(n+4)) -
                (6*(n+1)^2*(n+2))/((N+1)*(N-n)*(n+3)*(n+4)) +
                (6*(n+1)^2*(n+2))/((n+3)*(n+4)*k*(n-k+1)) -
                (6*(5*n+11))/((n+3)*(n+4))

  } else {

    #Determine unique values in population
    POP <- sort(population)
    POP.UNIQUE <- unique(POP)
    if (length(population) != N) {
      warning('Input N not used since it is inconsistent with input population') }

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

    #Set moments
    MEAN <- sum(POP.UNIQUE*exp(LOGPROBS))
    VAR  <- sum(((POP.UNIQUE - MEAN)^2)*exp(LOGPROBS))
    SKEW <- sum(((POP.UNIQUE - MEAN)^3)*exp(LOGPROBS))/(VAR^(3/2))
    KURT <- sum(((POP.UNIQUE - MEAN)^4)*exp(LOGPROBS))/(VAR^2) }

  #Generate output
  if (include.sd) {
    OUT <- data.frame(mean = MEAN, var = VAR, sd = sqrt(VAR), skew = SKEW, 
                      kurt = KURT, excess.kurt = KURT-3) } else {
    OUT <- data.frame(mean = MEAN, var = VAR, skew = SKEW, 
                      kurt = KURT, excess.kurt = KURT-3) }

  #Return the output
  OUT }
