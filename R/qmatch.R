#' Quantile function of the generalised matching distribution
#'
#' \code{qmatch} returns the quantiles for the arguments.
#'
#' This function computes quantiles from the quantile function of the generalised matching distribution.  Further 
#' details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2021) A generalised matching distribution for the problem of coincidences.
#'
#' @usage \code{pmatch(p, size, prob = 0, lower.tail = TRUE, log.p = FALSE)}
#' @param p A vector of numeric probability/log-probability values to be used as arguments for the quantile function
#' @param size The size parameter for the generalised matching distribution (number of objects to match)
#' @param prob The probability parameter for the generalised matching distribution (probability of known match)
#' @param lower.tail A logical value specifying whether the input represents lower or upper tail probabilities
#' @param log.p A logical value specifying whether the input values in p are log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be a vector of quantiles corresponding to the vector argument p

qmatch <- function(p, size, prob = 0, lower.tail = TRUE, log.p = FALSE) {

  #Check inputs
  if (!is.numeric(p))                       stop('Error: Input p should be numeric')
  if (!is.numeric(size))                    stop('Error: Size parameter should be a non-negative integer')
  if (!is.numeric(prob))                    stop('Error: Probability parameter should be numeric')
  if (!is.logical(log.p))                   stop('Error: Input log.p should be a single logical value')
  if (!is.logical(lower.tail))              stop('Error: Input lower.tail should be a single logical value')
  if (length(size) != 1)                    stop('Error: Size parameter should be a single non-negative integer')
  if (length(prob) != 1)                    stop('Error: Probability parameter should be a single numeric value')
  if (length(log.p) != 1)                   stop('Error: Input log.p should be a single logical value')
  if (length(lower.tail) != 1)              stop('Error: Input lower.tail should be a single logical value')
  if (log.p) {
    if (max(p) > 0)                         stop('Error: Log-probabilities in p must not be positive') }
  if (!log.p) {
    if (min(p) < 0)                         stop('Error: Probabilities in p must be between zero and one')
    if (max(p) > 1)                         stop('Error: Probabilities in p must be between zero and one') }
  if (size < Inf) { n <- as.integer(size) } else { n <- Inf }
  if (n != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                stop('Error: Size parameter should be a non-negative integer')
  if (prob < 0)                             stop('Error: Probability parameter should be between zero and one')
  if (prob > 1)                             stop('Error: Probability parameter should be between zero and one')

  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################
  
  #Deal with trivial case where x is empty
  if (length(p) == 0) {
    return(numeric(0)) }
  
  #Deal with trivial case where n = 0
  #Distribution is a point-mass on zero
  if (n == 0) {
    OUT <- rep(0, length(p))
    return(OUT) }
  
  #Deal with trivial case where n = 1
  #Distribution is a point-mass on one
  if (n == 1) {
    OUT <- rep(1, length(p))
    return(OUT) }
  
  #Deal with special case where prob = 1
  #Distribution is a point-mass on n
  if (prob == 1) {
    OUT <- rep(n, length(p))
    return(OUT) }

  #Deal with special case where n = Inf and prob = 0
  #Distribution is the Poisson distribution with unit rate
  if ((n == Inf)&(prob == 0)) {
    OUT <- qpois(p, lambda = 1, lower.tail = lower.tail, log.p = TRUE)
    return(OUT) }

  #Deal with special case where n = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((n == Inf)&(prob > 0)) {
    OUT <- rep(Inf, length(p))
    return(OUT) }

  ########################################################################################
  ############################### NON-TRIVIAL CASES ######################################
  ########################################################################################

  #Deal with non-trivial case where 1 < n < Inf and prob = 0
  #Set up vector of log-probabilities for the distribution
  if (prob == 0) {
    LOGPROBS <- rep(-Inf, n+1)
    LOGPROBS[n+1] <- -lfactorial(n)
    for (i in 1:n) {
      k <- n-i
      T1 <- log(n-k-1) + LOGPROBS[k+2]
      T2 <- ifelse(k < n-1, log(k+2) + LOGPROBS[k+3], -Inf)
      LOGPROBS[k+1] <- log(k+1) - log(n-k) + matrixStats::logSumExp(c(T1, T2)) }
    LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS) }

  #Deal with non-trivial where 0 < prob < 1
  #Set up vector of log-probabilities for the distribution
  if ((prob > 0)&(prob < 1)) {

    #Compute a matrix of classical log-probability values
    BASE.LOGPROBS <- matrix(-Inf, nrow = n+1, ncol = n+1)
    rownames(BASE.LOGPROBS) <- sprintf('size[%s]',  0:n)
    colnames(BASE.LOGPROBS) <- sprintf('match[%s]', 0:n)
    BASE.LOGPROBS[1,1] <- 0
    for (nn in 1:n) {
      BASE.LOGPROBS[nn+1, nn+1] <- -lfactorial(nn)
      for (i in 1:nn) {
        k <- nn-i
        T1 <- log(nn-k-1) + BASE.LOGPROBS[nn+1, k+2]
        T2 <- ifelse(k < nn-1, log(k+2) + BASE.LOGPROBS[nn+1, k+3], -Inf)
        BASE.LOGPROBS[nn+1, k+1] <- log(k+1) - log(nn-k) +
          matrixStats::logSumExp(c(T1, T2)) }
      BASE.LOGPROBS[nn+1, ] <- BASE.LOGPROBS[nn+1, ] -
        matrixStats::logSumExp(BASE.LOGPROBS[nn+1, ]) }

    #Compute vector of log-probability values
    LOGPROBS <- rep(-Inf, n+1)
    for (k in 0:n) {
      T1 <- dbinom(0:k, size = n, prob = prob, log = TRUE)
      T2 <- rep(-Inf, k+1)
      for (l in 0:k) { T2[l+1] <- BASE.LOGPROBS[n-l+1, k-l+1] }
      LOGPROBS[k+1] <- matrixStats::logSumExp(T1 + T2) }
    LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS) }
  
  #Compute vector of cumulative log-probabilities
  CUMLOGPROBS <- rep(-Inf, n+1)
  CUMLOGPROBS[1] <- LOGPROBS[1]
  for (i in 1:n) {
    CUMLOGPROBS[i+1] <- matrixStats::logSumExp(c(CUMLOGPROBS[i], LOGPROBS[i+1])) }
  
  #Compute quantiles
  QUANTILES <- rep(0, length(p))
  if (log.p) { LOGP <- p } else { LOGP <- log(p) }
  if (!lower.tail) { LOGP <- VGAM::log1mexp(-LOGP) }
  for (i in 1:length(p)) { 
    LL <- LOGP[i]
    QUANTILES[i] <- sum(LL > CUMLOGPROBS) }

  #Return output
  QUANTILES }
