#' Highest density region (HDR) for the generalised matching distribution
#'
#' \code{HDR.match} returns the highest density region (HDR) for the arguments.
#'
#' This function computes the highest density region (HDR) for the generalised matching distribution at a 
#' specified coverage probability level.  Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2021) A generalised matching distribution for the problem of coincidences.
#'
#' @usage \code{HDR.match(cover.prob, size, prob = 0)}
#' @param cover.prob The minimum coverage probability for the HDR
#' @param size The size parameter for the generalised matching distribution (number of objects to match)
#' @param prob The probability parameter for the generalised matching distribution (probability of known match)
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be an object of class \code{'hdr'} giving the highest density region for the distribution

HDR.match <- function(cover.prob, size, prob = 0) {

  #Check inputs
  if (!is.numeric(size))                    stop('Error: Size parameter should be a non-negative integer')
  if (!is.numeric(prob))                    stop('Error: Probability parameter should be numeric')
  if (length(size) != 1)                    stop('Error: Size parameter should be a single non-negative integer')
  if (length(prob) != 1)                    stop('Error: Probability parameter should be a single numeric value')
  if (size < Inf) { n <- as.integer(size) } else { n <- Inf }
  if (n != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                stop('Error: Size parameter should be a non-negative integer')
  if (prob < 0)                             stop('Error: Probability parameter should be between zero and one')
  if (prob > 1)                             stop('Error: Probability parameter should be between zero and one')

  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################

  #Deal with trivial case where n = 0
  #Distribution is a point-mass on zero
  if (n == 0) { 
    LOGPROBS <- 0 }

  #Deal with trivial case where n = 1
  #Distribution is a point-mass on one
  if (n == 1) { 
    LOGPROBS <- c(-Inf, 0) }

  #Deal with special case where prob = 1
  #Distribution is a point-mass on n
  if (prob == 1) {
    LOGPROBS <- rep(-Inf, n+1)
    LOGPROBS[n+1] <- 0 }

  #Deal with special case where n = Inf and prob = 0
  #Distribution is the Poisson distribution with unit rate
  if ((n == Inf)&(prob == 0)) {
    LOGPROBS <- dpois(0:n, lambda = 1, log = TRUE) }

  #Deal with special case where n = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((n == Inf)&(prob > 0)) {
    stop('Error: Distribution is a point-mass on infinity') }

  ########################################################################################
  ############################### NON-TRIVIAL CASES ######################################
  ########################################################################################

  #Deal with non-trivial case where 1 < n < Inf and prob = 0
  #Set up vector of log-probabilities for the distribution
  if ((n > 1)&(n < Inf)) {
  if (prob == 0) {
    LOGPROBS <- rep(-Inf, n+1)
    LOGPROBS[n+1] <- -lfactorial(n)
    for (i in 1:n) {
      k <- n-i
      T1 <- log(n-k-1) + LOGPROBS[k+2]
      T2 <- ifelse(k < n-1, log(k+2) + LOGPROBS[k+3], -Inf)
      LOGPROBS[k+1] <- log(k+1) - log(n-k) + matrixStats::logSumExp(c(T1, T2)) }
    LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS) } }

  #Deal with non-trivial where 1 < n < Inf and 0 < prob < 1
  #Set up vector of log-probabilities for the distribution
  if ((n > 1)&(n < Inf)) {
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
    LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS) } }

  #Set density function
  DENS <- function(x) {
    if (x %in% 0:n) { exp(LOGPROBS[x+1]) } else { 0 } }
  
  #Set distribution name
  DIST <- ifelse(prob == 0, 
                 paste0('classical matching distribution with ', n, ' trials'),
                 paste0('matching distribution with ', n, 
                        ' trials and matching probability ', round(prob, 4)))
    
  #Return HDR output
  HDR <- stat.extend::HDR.discrete(cover.prob = cover.prob, f = DENS, 
                                   supp.min = 0, supp.max = n,
                                   distribution = DIST)
  HDR }
