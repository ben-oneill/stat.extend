#' The Generalized Matching Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the generalized matching distribution with parameters \code{size} and \code{prob}. 
#'
#'
#'
#' @param x,q A vector of numeric values to be used as arguments for the mass function
#' @param p vector of probabilities
#' @param n number of observations
#' @param size The size parameter for the generalised matching distribution (number of objects to match)
#' @param prob The probability parameter for the generalised matching distribution (probability of known match)
#' @param log,log.p A logical value specifying whether results should be returned as log-probabilities
#' @param lower.tail A logical value specifying whether the input represents lower or upper tail probabilities
#' 
#' @returns 
#' 
#' \code{dmatching} gives the density, \code{pmatching} gives the distribution function,
#' \code{qmatching} gives the quantile function and \code{rmatching} generates random deviates.
#' 
#' @references 
#' O'Neill, B. (2021) A generalised matching distribution for the problem of coincidences.
#' 
#' @examples 
#' x <- rmatching(1000, 5)
#' table(x)
#' # No Fours!
#' dmatching(0:5, 5)
#' 
#' @name Matching
#' @rdname Matching
dmatching <- function(x, size, prob = 0, log = FALSE) {

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

    #Set matrix M
    M <- matrix(-Inf, nrow = n+1, ncol = n+1)
    rownames(M) <- sprintf('k[%s]', 0:n)
    colnames(M) <- sprintf('l[%s]', 0:n)
    for (k in 0:n) {
    for (l in 0:k) {
      M[k+1, l+1] <- BASE.LOGPROBS[n-l+1, k-l+1] } }

    #Compute vector of log-probability values
    LOGPROBS <- rep(-Inf, n+1)
    LOGBINOM <- dbinom(0:n, size = n, prob = prob, log = TRUE)
    for (k in 0:n) {
      LOGPROBS[k+1] <- matrixStats::logSumExp(LOGBINOM + M[k+1, ]) }
    LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS) }

  #Compute output
  OUT <- rep(-Inf, length(x))
  for (i in 1:length(x)) {
    xx <- x[i]
    if (xx %in% 0:n) { OUT[i] <- LOGPROBS[xx+1] } }

  #Return output
  if (log) { OUT } else { exp(OUT) } }

#' @rdname Matching
pmatching <- function(q, size, prob = 0, lower.tail = TRUE, log.p = FALSE) {

  #Check inputs
  x <- q
  if (!is.numeric(x))                       stop('Error: Input q should be numeric')
  if (!is.numeric(size))                    stop('Error: Size parameter should be a non-negative integer')
  if (!is.numeric(prob))                    stop('Error: Probability parameter should be numeric')
  if (!is.logical(log.p))                   stop('Error: Input log.p should be a single logical value')
  if (!is.logical(lower.tail))              stop('Error: Input lower.tail should be a single logical value')
  if (length(size) != 1)                    stop('Error: Size parameter should be a single non-negative integer')
  if (length(prob) != 1)                    stop('Error: Probability parameter should be a single numeric value')
  if (length(log.p) != 1)                   stop('Error: Input log.p should be a single logical value')
  if (length(lower.tail) != 1)              stop('Error: Input lower.tail should be a single logical value')
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
      if (x[i] >= 0) { OUT[i] <- 0 } }
    if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }
  
  #Deal with trivial case where n = 1
  #Distribution is a point-mass on one
  if (n == 1) {
    OUT <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      if (x[i] >= 1) { OUT[i] <- 0 } }
    if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }
  

  #Deal with special case where prob = 1
  #Distribution is a point-mass on n
  if (prob == 1) {
    OUT <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      if (x[i] >= n) { OUT[i] <- 0 } }
    if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with special case where n = Inf and prob = 0
  #Distribution is the Poisson distribution with unit rate
  if ((n == Inf)&(prob == 0)) {
    OUT <- ppois(x, lambda = 1, lower.tail = lower.tail, log.p = TRUE)
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with special case where n = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((n == Inf)&(prob > 0)) {
    OUT <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      if (x[i] == Inf) { OUT[i] <- 0 } }
    if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }

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

    #Set matrix M
    M <- matrix(-Inf, nrow = n+1, ncol = n+1)
    rownames(M) <- sprintf('k[%s]', 0:n)
    colnames(M) <- sprintf('l[%s]', 0:n)
    for (k in 0:n) {
    for (l in 0:k) {
      M[k+1, l+1] <- BASE.LOGPROBS[n-l+1, k-l+1] } }

    #Compute vector of log-probability values
    LOGPROBS <- rep(-Inf, n+1)
    LOGBINOM <- dbinom(0:n, size = n, prob = prob, log = TRUE)
    for (k in 0:n) {
      LOGPROBS[k+1] <- matrixStats::logSumExp(LOGBINOM + M[k+1, ]) }
    LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS) }
  
  #Compute vector of cumulative log-probabilities
  CUMLOGPROBS <- rep(-Inf, n+1)
  CUMLOGPROBS[1] <- LOGPROBS[1]
  for (i in 1:n) {
    CUMLOGPROBS[i+1] <- matrixStats::logSumExp(c(CUMLOGPROBS[i], LOGPROBS[i+1])) }
  
  #Compute output
  OUT <- rep(-Inf, length(x))
  for (i in 1:length(x)) {
    xx <- floor(x[i])
    if ((xx >= 0)&(xx < n)) { OUT[i] <- CUMLOGPROBS[xx+1] }
    if (xx >= n) { OUT[i] <- 0 } }
  if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
  
  #Return output
  if (log.p) { OUT } else { exp(OUT) } }

#' @rdname Matching
qmatching <- function(p, size, prob = 0, lower.tail = TRUE, log.p = FALSE) {

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

    #Set matrix M
    M <- matrix(-Inf, nrow = n+1, ncol = n+1)
    rownames(M) <- sprintf('k[%s]', 0:n)
    colnames(M) <- sprintf('l[%s]', 0:n)
    for (k in 0:n) {
    for (l in 0:k) {
      M[k+1, l+1] <- BASE.LOGPROBS[n-l+1, k-l+1] } }

    #Compute vector of log-probability values
    LOGPROBS <- rep(-Inf, n+1)
    LOGBINOM <- dbinom(0:n, size = n, prob = prob, log = TRUE)
    for (k in 0:n) {
      LOGPROBS[k+1] <- matrixStats::logSumExp(LOGBINOM + M[k+1, ]) }
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

#' @rdname Matching
rmatching <- function(n, size, prob = 0) {
  
  #Check inputs
  if (!is.numeric(n))                        stop('Error: Input n should be numeric')
  if (!is.numeric(size))                     stop('Error: Size parameter should be a non-negative integer')
  if (length(n) != 1)                        stop('Error: Input n should be a single non-negative integer')
  if (length(size) != 1)                     stop('Error: Size parameter should be a single non-negative integer')
  if (size < Inf) { nn <- as.integer(size) } else { nn <- Inf }
  if (nn != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                 stop('Error: Input n should be a single non-negative integer')
  
  #Deal with trivial case where n = 0
  if (n == 0) { return(numeric(0)) }
  
  #Generate output vector
  OUT   <- rep(0, n)
  KNOWN <- rbinom(n, size = size, prob = prob)
  for (i in 1:n) {
    nnn    <- nn-KNOWN[i]
    PERM   <- sample.int(nnn, size = nnn, replace = FALSE)
    OUT[i] <- KNOWN[i] + sum(PERM == 1:nnn) }
  
  #Return output
  OUT }
