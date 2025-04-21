#' Quantile function of the Finite Population Order Statistic (FPOS) distribution
#'
#' \code{qFPOS} returns the quantile values for the arguments.
#'
#' This function computes the quantiles from the cumulative distribution function of the Finite Population Order
#' Statistic (FPOS) distribution.  This is the distribution of the order statistics when sampling without
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
#' @usage \code{qFPOS(p, population = NULL, k, n, N = length(population), lower.tail = TRUE, log.p = FALSE)}
#' @param p A vector of probabilities or log-probabilities to be used as arguments for the mass function
#' @param population Optional vector specifying the population values
#' @param k The rank-order (a positive integer no greater than the sample size)
#' @param n The sample size (a positive integer no greater than the population size)
#' @param N The population size (a positive integer)
#' @param lower.tail A logical value specifying whether the input probabilities are for the  cumulative distribution
#' function or the corresponding survival function
#' @param log.p A logical value specifying whether the input values in p are log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be a vector of quantiles corresponding to the vector argument p

qFPOS <- function(p, population = NULL, k, n, N = length(population), lower.tail = TRUE, log.p = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(p))                       stop('Error: Argument p is not numeric')
  if (!is.numeric(N))                       stop('Error: Input parameter N is not numeric')
  if (!is.numeric(n))                       stop('Error: Input parameter n is not numeric')
  if (!is.numeric(k))                       stop('Error: Input parameter k is not numeric')
  if (!is.logical(lower.tail))              stop('Error: lower.tail option is not a logical value')
  if (!is.logical(log.p))                   stop('Error: log.p option is not a logical value')

  #Check that parameters are atomic
  if (length(N) != 1)                       stop('Error: Input parameter N should be a single number')
  if (length(n) != 1)                       stop('Error: Input parameter n should be a single number')
  if (length(k) != 1)                       stop('Error: Input parameter k should be a single number')
  if (length(lower.tail) != 1)              stop('Error: lower.tail option should be a single logical value')
  if (length(log.p) != 1)                   stop('Error: log.p option should be a single logical value')

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

  #Check the input p and set the log-probabilities
  if (log.p) {
    if (max(p) > 0)                         stop('Error: Input p contains log-probabilies that are greater than zero')
    LOGP <- p
  } else {
    if (min(p) < 0)                         stop('Error: Input p contains probabilities below zero')
    if (max(p) > 1)                         stop('Error: Input p contains probabilities above one')
    LOGP <- log(p) }
  if (!lower.tail) { LOGP <- log1mexp(-LOGP) }
  
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
  
  #Create cumulative log-probabilities for population values
  LOGCDF  <- rep(-Inf, NNpop)
  LOGCDF[1] <- LOGPROBS[1]
  if (NNpop > 1) {
    for (i in 2:NNpop) {
      LOGCDF[i]  <- matrixStats::logSumExp(c(LOGCDF[i-1], LOGPROBS[i])) } }
  LOGCDF <- LOGCDF - LOGCDF[NNpop]
  
  #Create output vector
  np <- length(p)
  QUANTILES <- rep(-Inf, np)
  for (i in 1:np) { 
    if (LOGP[i] == -Inf) { QUANTILES[i] <- POP.UNIQUE[min(which(LOGCDF > -Inf))] }
    if (LOGP[i] == 0)    { QUANTILES[i] <- POP.UNIQUE[min(which(LOGCDF == 0))] }
    if ((LOGP[i] > -Inf)&(LOGP[i] < 0)) { 
      QUANTILES[i] <- POP.UNIQUE[1 + sum(LOGP[i] >= LOGCDF)] } }

  #Return output
  QUANTILES }
