#' Mass function of the generalised matching distribution
#'
#' \code{dmatch} returns the probability or log-probability values for the arguments.
#'
#' This function computes probabilities or log-probabilities from the mass function of the generalised matching
#' distribution.  Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2021) A generalised matching distribution for the problem of coincidences.
#'
#' @param x A vector of numeric values to be used as arguments for the mass function
#' @param size The size parameter for the generalised matching distribution (number of objects to match)
#' @param prob The probability parameter for the generalised matching distribution (probability of known match)
#' @param log A logical value specifying whether results should be returned as log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be a vector of probabilities/log-probabilities corresponding to the vector argument x

dmatch <- function(x, size, prob = 0, log = FALSE) {

  #Check inputs
  if (!is.numeric(x))                       stop('Error: Input x should be numeric')
  if (!is.numeric(size))                    stop('Error: Size parameter should be a non-negative integer')
  if (!is.numeric(prob))                    stop('Error: Probability parameter should be numeric')
  if (!is.logical(log))                     stop('Error: Input log should be a single logical value')
  if (length(size) != 1)                    stop('Error: Size parameter should be a single non-negative integer')
  if (length(prob) != 1)                    stop('Error: Probability parameter should be a single numeric value')
  if (length(log) != 1)                     stop('Error: Input log should be a single logical value')
  if (size < Inf) { n <- as.integer(size) } else { n <- Inf }
  if (n != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                stop('Error: Size parameter should be a non-negative integer')
  if (prob < 0)                             stop('Error: Probability parameter should be between zero and one')
  if (prob > 1)                             stop('Error: Probability parameter should be between zero and one')

  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################

  #Deal with trivial case where x is empty
  if (length(x) == 0) {
    return(numeric(0)) }
  
  #Deal with trivial case where n = 0
  #Distribution is a point-mass on zero
  if (n == 0) {
    OUT <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      if (x[i] == 0) { OUT[i] <- 0 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with trivial case where n = 1
  #Distribution is a point-mass on one
  if (n == 1) {
    OUT <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      if (x[i] == 1) { OUT[i] <- 0 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with special case where prob = 1
  #Distribution is a point-mass on n
  if (prob == 1) {
    OUT <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      if (x[i] == n) { OUT[i] <- 0 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with special case where n = Inf and prob = 0
  #Distribution is the Poisson distribution with unit rate
  if ((n == Inf)&(prob == 0)) {
    OUT <- dpois(x, lambda = 1, log = TRUE)
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with special case where n = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((n == Inf)&(prob > 0)) {
    OUT <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      if (x[i] == Inf) { OUT[i] <- 0 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }

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

  #Compute output
  OUT <- rep(-Inf, length(x))
  for (i in 1:length(x)) {
    xx <- x[i]
    if (xx %in% 0:n) { OUT[i] <- LOGPROBS[xx+1] } }

  #Return output
  if (log) { OUT } else { exp(OUT) } }
