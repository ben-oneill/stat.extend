#' The Generalized Matching Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the generalized matching distribution with parameters \code{size}, \code{trials} and \code{prob}. 
#' The distribution is for the total number of matches over all trials.  In each trial the player 
#' initially matches objects independently with probability \code{prob} and then allocates remaining
#' objects using a random permutation.
#'
#' @param x,q A vector of numeric values to be used as arguments for the mass function
#' @param p vector of probabilities
#' @param n number of observations
#' @param size The size parameter for the generalised matching distribution (number of objects to match)
#' @param trials The trials parameter for the generalised matching distribution (number of times the matching game is repeated)
#' @param prob The probability parameter for the generalised matching distribution (probability of known match)
#' @param log,log.p A logical value specifying whether results should be returned as log-probabilities
#' @param lower.tail A logical value specifying whether the input represents lower or upper tail probabilities
#' @param approx A logical value specifying whether to use the normal approximation to the distribution
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
#' tabulate(x)
#' # No Fours!
#' # This is actually one of the key properties of the matching distribution.  
#' # With size parameter n the distribution has support 0,1,2,...,n-2,n (i.e., it 
#' # cannot give outcome n-1). The reason for this is that in a permutation it 
#' # is impossible to give n-1 matches. 
#' # If there are n-1 matches then the last object in the permutation must also be a match.
#' dmatching(0:5, 5)
#' 
#' @name Matching
#' @rdname Matching
dmatching <- function(x, size, trials = 1, prob = 0, log = FALSE, approx = (trials > 100)) {

  #Check inputs
  if (!is.numeric(x))                       stop('Error: Input x should be numeric')
  if (!is.numeric(trials))                  stop('Error: Input trials should be a positive integer')
  if (!is.numeric(size))                    stop('Error: Size parameter should be a non-negative integer')
  if (!is.numeric(prob))                    stop('Error: Probability parameter should be numeric')
  if (!is.logical(log))                     stop('Error: Input log should be a single logical value')
  if (!is.logical(approx))                  stop('Error: Input approx should be a single logical value')
  if (length(trials) != 1)                  stop('Error: Input trials should be a single positive integer')
  if (length(size) != 1)                    stop('Error: Size parameter should be a single non-negative integer')
  if (length(prob) != 1)                    stop('Error: Probability parameter should be a single numeric value')
  if (length(log) != 1)                     stop('Error: Input log should be a single logical value')
  if (length(approx) != 1)                  stop('Error: Input approx should be a single logical value')
  if (trials == Inf)                        stop('Error: Input trials must be finite')
  m <- as.integer(trials)
  if (m != trials)                          stop('Error: Input trials should be a positive integer')
  if (m <= 0)                               stop('Error: Input trials should be a positive integer')
  if (size < Inf) { n <- as.integer(size) } else { n <- Inf }
  if (n != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                stop('Error: Size parameter should be a non-negative integer')
  if (prob < 0)                             stop('Error: Probability parameter should be between zero and one')
  if (prob > 1)                             stop('Error: Probability parameter should be between zero and one')

  #Set sample total
  T <- x
  
  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################

  #Deal with trivial case where T is empty
  if (length(T) == 0) {
    return(numeric(0)) }
  
  #Deal with trivial case where n = 0
  #Distribution is a point-mass on zero
  if (n == 0) {
    OUT <- rep(-Inf, length(T))
    for (i in 1:length(T)) {
      if (T[i] == 0) { OUT <- 0 } else { OUT <- -Inf } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with trivial case where n = 1
  #Distribution is a point-mass on m
  if (n == 1) {
    OUT <- rep(-Inf, length(T))
    for (i in 1:length(T)) {
      if (T[i] == m) { OUT[i] <- 0 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with special case where prob = 1
  #Distribution is a point-mass on n*m
  if (prob == 1) {
    OUT <- rep(-Inf, length(T))
    for (i in 1:length(T)) {
      if (T[i] == n*m) { OUT[i] <- 0 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }
  
  #Deal with case where n = Inf and prob = 0
  #Distribution is Poisson with rate m
  if ((n == Inf)&(prob == 0)) {
    OUT <- dpois(T, lambda = m, log = TRUE)
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with special case where n = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((n == Inf)&(prob > 0)) {
    OUT <- rep(-Inf, length(T))
    for (i in 1:length(T)) {
      if (T[i] == Inf) { OUT[i] <- 0 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }
  
  ########################################################################################
  ############################### NON-TRIVIAL CASES ######################################
  ########################################################################################

  #Compute distribution of sample total
  if (approx) {
    
    #Compute normal approximation to generalised matching distribution
    MEAN <- 1 + n*prob - prob^n
    VAR  <- 1 - prob^(2*n) + n*(prob - prob^2 - prob^(n-1) - prob^n + 2*prob^(n+1))
    LOGPROBS.TOTAL <- dnorm(0:(n*m), mean = m*MEAN, sd = sqrt(m*VAR), log = TRUE)
    LOGPROBS.TOTAL <- LOGPROBS.TOTAL - matrixStats::logSumExp(LOGPROBS.TOTAL)
    
  } else {
    
    #Deal with non-trivial case where 1 < n < Inf and prob = 0
    #Set up vector of log-probabilities for the distribution
    if ((n < Inf)&(prob == 0)) {
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
    
    #Compute log-probabilities for sample total
    LOGPROBS.TOTAL.ALL <- matrix(-Inf, nrow = m, ncol = n*m+1)
    LOGPROBS.TOTAL.ALL[1, 0:n+1] <- LOGPROBS
    if (m > 1) {
      for (t in 2:m) {
        for (s in 0:(n*t)) {
          k  <- 0:min(s, n)
          s0 <- s - k
          T1 <- LOGPROBS.TOTAL.ALL[t-1, s0+1]
          T2 <- LOGPROBS[k+1]
          LOGPROBS.TOTAL.ALL[t, s+1] <- matrixStats::logSumExp(T1 + T2) }
        LOGPROBS.TOTAL.ALL[t, ] <- LOGPROBS.TOTAL.ALL[t, ] - 
          matrixStats::logSumExp(LOGPROBS.TOTAL.ALL[t, ]) } } 
    LOGPROBS.TOTAL <- LOGPROBS.TOTAL.ALL[m, ] }
  
  #Compute output
  OUT <- rep(-Inf, length(T))
  for (i in 1:length(T)) {
    TT <- T[i]
    if (TT %in% 0:(n*m)) { OUT[i] <- LOGPROBS.TOTAL[TT+1] } }

  #Return output
  if (log) { OUT } else { exp(OUT) } }

#' @rdname Matching
pmatching <- function(q, size, trials = 1, prob = 0, lower.tail = TRUE, log.p = FALSE, approx = (trials > 100)) {

  #Check inputs
  if (!is.numeric(q))                       stop('Error: Input q should be numeric')
  if (!is.numeric(trials))                  stop('Error: Input trials should be a positive integer')
  if (!is.numeric(size))                    stop('Error: Size parameter should be a non-negative integer')
  if (!is.numeric(prob))                    stop('Error: Probability parameter should be numeric')
  if (!is.logical(lower.tail))              stop('Error: Input lower.tail should be a single logical value')
  if (!is.logical(log.p))                   stop('Error: Input log.p should be a single logical value')
  if (!is.logical(approx))                  stop('Error: Input approx should be a single logical value')
  if (length(trials) != 1)                  stop('Error: Input trials should be a single positive integer')
  if (length(size) != 1)                    stop('Error: Size parameter should be a single non-negative integer')
  if (length(prob) != 1)                    stop('Error: Probability parameter should be a single numeric value')
  if (length(lower.tail) != 1)              stop('Error: Input lower.tail should be a single logical value')
  if (length(log.p) != 1)                   stop('Error: Input log.p should be a single logical value')
  if (length(approx) != 1)                  stop('Error: Input approx should be a single logical value')
  if (trials == Inf)                        stop('Error: Input trials must be finite')
  m <- as.integer(trials)
  if (m != trials)                          stop('Error: Input trials should be a positive integer')
  if (m <= 0)                               stop('Error: Input trials should be a positive integer')
  if (size < Inf) { n <- as.integer(size) } else { n <- Inf }
  if (n != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                stop('Error: Size parameter should be a non-negative integer')
  if (prob < 0)                             stop('Error: Probability parameter should be between zero and one')
  if (prob > 1)                             stop('Error: Probability parameter should be between zero and one')

  #Set sample total
  T <- q
  
  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################

  #Deal with trivial case where T is empty
  if (length(T) == 0) {
    return(numeric(0)) }
  
  #Deal with trivial case where n = 0
  #Distribution is a point-mass on zero
  if (n == 0) {
    OUT <- rep(-Inf, length(T))
    for (i in 1:length(T)) {
      if (T[i] >= 0) { OUT[i] <- 0 } }
    if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }
  
  #Deal with trivial case where n = 1
  #Distribution is a point-mass on m
  if (n == 1) {
    OUT <- rep(-Inf, length(T))
    for (i in 1:length(T)) {
      if (T[i] >= m) { OUT[i] <- 0 } }
    if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }
  
  #Deal with special case where prob = 1
  #Distribution is a point-mass on n*m
  if (prob == 1) {
    OUT <- rep(-Inf, length(T))
    for (i in 1:length(T)) {
      if (T[i] >= n*m) { OUT[i] <- 0 } }
    if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }
  
  #Deal with case where n = Inf and prob = 0
  #Distribution is Poisson with rate m
  if ((n == Inf)&(prob == 0)) {
    OUT <- ppois(T, lambda = m, lower.tail = lower.tail, log.p = TRUE)
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with special case where n = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((n == Inf)&(prob > 0)) {
    OUT <- rep(-Inf, length(T))
    for (i in 1:length(T)) {
      if (T[i] == Inf) { OUT[i] <- 0 } }
    if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }
    if (log.p) { return(OUT) } else { return(exp(OUT)) } }
  
  ########################################################################################
  ############################### NON-TRIVIAL CASES ######################################
  ########################################################################################

  #Compute distribution of sample total
  if (approx) {
    
    #Compute normal approximation to generalised matching distribution
    MEAN <- 1 + n*prob - prob^n
    VAR  <- 1 - prob^(2*n) + n*(prob - prob^2 - prob^(n-1) - prob^n + 2*prob^(n+1))
    LOGPROBS.TOTAL <- dnorm(0:(n*m), mean = m*MEAN, sd = sqrt(m*VAR), log = TRUE)
    LOGPROBS.TOTAL <- LOGPROBS.TOTAL - matrixStats::logSumExp(LOGPROBS.TOTAL)
    
  } else {
    
    #Deal with case where n = Inf and prob = 0
    if ((n == Inf)&(prob == 0)) {
      LOGPROBS <- dpois(q, lambda = 1, log = TRUE) }
    
    #Deal with non-trivial case where 1 < n < Inf and prob = 0
    #Set up vector of log-probabilities for the distribution
    if ((n < Inf)&(prob == 0)) {
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
    
    #Compute log-probabilities for sample total
    LOGPROBS.TOTAL.ALL <- matrix(-Inf, nrow = m, ncol = n*m+1)
    LOGPROBS.TOTAL.ALL[1, 0:n+1] <- LOGPROBS
    if (m > 1) {
      for (t in 2:m) {
        for (s in 0:(n*t)) {
          k  <- 0:min(s, n)
          s0 <- s - k
          T1 <- LOGPROBS.TOTAL.ALL[t-1, s0+1]
          T2 <- LOGPROBS[k+1]
          LOGPROBS.TOTAL.ALL[t, s+1] <- matrixStats::logSumExp(T1 + T2) }
        LOGPROBS.TOTAL.ALL[t, ] <- LOGPROBS.TOTAL.ALL[t, ] - 
          matrixStats::logSumExp(LOGPROBS.TOTAL.ALL[t, ]) } } 
    LOGPROBS.TOTAL <- LOGPROBS.TOTAL.ALL[m, ] }
  
  #Compute vector of cumulative log-probabilities
  CUMLOGPROBS.TOTAL <- rep(-Inf, n*m+1)
  CUMLOGPROBS.TOTAL[1] <- LOGPROBS.TOTAL[1]
  for (i in 1:(n*m)) {
    CUMLOGPROBS.TOTAL[i+1] <- matrixStats::logSumExp(c(CUMLOGPROBS.TOTAL[i], LOGPROBS.TOTAL[i+1])) }
  
  #Compute output
  OUT <- rep(-Inf, length(T))
  for (i in 1:length(T)) {
    TT <- floor(T[i])
    if ((TT >= 0)&(TT < n*m)) { OUT[i] <- CUMLOGPROBS.TOTAL[TT+1] }
    if (TT >= n*m) { OUT[i] <- 0 } }
  if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }

  #Return output
  if (log.p) { OUT } else { exp(OUT) } }

#' @rdname Matching
qmatching <- function(p, size, trials = 1, prob = 0, lower.tail = TRUE, log.p = FALSE, approx = (trials > 100)) {

  #Check inputs
  if (!is.numeric(p))                       stop('Error: Input p should be numeric')
  if (!is.numeric(trials))                  stop('Error: Input trials should be a positive integer')
  if (!is.numeric(size))                    stop('Error: Size parameter should be a non-negative integer')
  if (!is.numeric(prob))                    stop('Error: Probability parameter should be numeric')
  if (!is.logical(lower.tail))              stop('Error: Input lower.tail should be a single logical value')
  if (!is.logical(log.p))                   stop('Error: Input log.p should be a single logical value')
  if (log.p) {
    if (max(p) > 0)                         stop('Error: Log-probabilities in p must not be positive') }
  if (!log.p) {
    if (min(p) < 0)                         stop('Error: Probabilities in p must be between zero and one')
    if (max(p) > 1)                         stop('Error: Probabilities in p must be between zero and one') }
  if (!is.logical(approx))                  stop('Error: Input approx should be a single logical value')
  if (length(trials) != 1)                  stop('Error: Input trials should be a single positive integer')
  if (length(size) != 1)                    stop('Error: Size parameter should be a single non-negative integer')
  if (length(prob) != 1)                    stop('Error: Probability parameter should be a single numeric value')
  if (length(lower.tail) != 1)              stop('Error: Input lower.tail should be a single logical value')
  if (length(log.p) != 1)                   stop('Error: Input log.p should be a single logical value')
  if (length(approx) != 1)                  stop('Error: Input approx should be a single logical value')
  if (trials == Inf)                        stop('Error: Input trials must be finite')
  m <- as.integer(trials)
  if (m != trials)                          stop('Error: Input trials should be a positive integer')
  if (m <= 0)                               stop('Error: Input trials should be a positive integer')
  if (size < Inf) { n <- as.integer(size) } else { n <- Inf }
  if (n != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                stop('Error: Size parameter should be a non-negative integer')
  if (prob < 0)                             stop('Error: Probability parameter should be between zero and one')
  if (prob > 1)                             stop('Error: Probability parameter should be between zero and one')
  
  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################

  #Deal with trivial case where p is empty
  if (length(p) == 0) {
    return(numeric(0)) }
  
  #Deal with trivial case where n = 0
  #Distribution is a point-mass on zero
  if (n == 0) {
    OUT <- rep(0, length(p))
    return(OUT) }
  
  #Deal with trivial case where n = 1
  #Distribution is a point-mass on m
  if (n == 1) {
    OUT <- rep(m, length(p))
    return(OUT) }
  
  #Deal with special case where prob = 1
  #Distribution is a point-mass on n*m
  if (prob == 1) {
    OUT <- rep(n*m, length(p))
    return(OUT) }
  
  #Deal with case where n = Inf and prob = 0
  #Distribution is Poisson with rate m
  if ((n == Inf)&(prob == 0)) {
    OUT <- qpois(p, lambda = m, lower.tail = lower.tail, log.p = log.p)
    return(OUT) }
  
  #Deal with special case where n = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((n == Inf)&(prob > 0)) {
    OUT <- rep(Inf, length(p))
    return(OUT) }
  
  ########################################################################################
  ############################### NON-TRIVIAL CASES ######################################
  ########################################################################################

  #Compute distribution of sample total
  if (approx) {
    
    #Compute normal approximation to generalised matching distribution
    MEAN <- 1 + n*prob - prob^n
    VAR  <- 1 - prob^(2*n) + n*(prob - prob^2 - prob^(n-1) - prob^n + 2*prob^(n+1))
    LOGPROBS.TOTAL <- dnorm(0:(n*m), mean = m*MEAN, sd = sqrt(m*VAR), log = TRUE)
    LOGPROBS.TOTAL <- LOGPROBS.TOTAL - matrixStats::logSumExp(LOGPROBS.TOTAL)
    
  } else {
    
    #Deal with non-trivial case where 1 < n < Inf and prob = 0
    #Set up vector of log-probabilities for the distribution
    if ((n < Inf)&(prob == 0)) {
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
    
    #Compute log-probabilities for sample total
    LOGPROBS.TOTAL.ALL <- matrix(-Inf, nrow = m, ncol = n*m+1)
    LOGPROBS.TOTAL.ALL[1, 0:n+1] <- LOGPROBS
    if (m > 1) {
      for (t in 2:m) {
        for (s in 0:(n*t)) {
          k  <- 0:min(s, n)
          s0 <- s - k
          T1 <- LOGPROBS.TOTAL.ALL[t-1, s0+1]
          T2 <- LOGPROBS[k+1]
          LOGPROBS.TOTAL.ALL[t, s+1] <- matrixStats::logSumExp(T1 + T2) }
        LOGPROBS.TOTAL.ALL[t, ] <- LOGPROBS.TOTAL.ALL[t, ] - 
          matrixStats::logSumExp(LOGPROBS.TOTAL.ALL[t, ]) } } 
    LOGPROBS.TOTAL <- LOGPROBS.TOTAL.ALL[m, ] }
  
  #Compute vector of cumulative log-probabilities
  CUMLOGPROBS.TOTAL <- rep(-Inf, n*m+1)
  CUMLOGPROBS.TOTAL[1] <- LOGPROBS.TOTAL[1]
  for (i in 1:(n*m)) {
    CUMLOGPROBS.TOTAL[i+1] <- matrixStats::logSumExp(c(CUMLOGPROBS.TOTAL[i], LOGPROBS.TOTAL[i+1])) }
  
  #Compute quantiles
  QUANTILES <- rep(0, length(p))
  if (log.p) { LOGP <- p } else { LOGP <- log(p) }
  if (!lower.tail) { LOGP <- VGAM::log1mexp(-LOGP) }
  for (i in 1:length(p)) { 
    LL <- LOGP[i]
    QUANTILES[i] <- sum(LL > CUMLOGPROBS.TOTAL) }

  #Return output
  QUANTILES }

#' @rdname Matching
rmatching <- function(n, size, trials = 1, prob = 0) {
  
  #Check inputs
  if (!is.numeric(n))                        stop('Error: Input n should be numeric')
  if (!is.numeric(trials))                   stop('Error: Input trials should be a positive integer')
  if (!is.numeric(size))                     stop('Error: Size parameter should be a non-negative integer')
  if (!is.numeric(prob))                     stop('Error: Probability parameter should be numeric')
  if (length(n) != 1)                        stop('Error: Input n should be a single non-negative integer')
  if (length(trials) != 1)                   stop('Error: Input trials should be a single positive integer')
  if (length(size) != 1)                     stop('Error: Size parameter should be a single non-negative integer')
  if (length(prob) != 1)                     stop('Error: Probability parameter should be a single numeric value')
  if (trials == Inf)                         stop('Error: Input trials must be finite')
  m <- as.integer(trials)
  if (m != trials)                           stop('Error: Input trials should be a positive integer')
  if (m <= 0)                                stop('Error: Input trials should be a positive integer')
  if (size < Inf) { nn <- as.integer(size) } else { nn <- Inf }
  if (nn != size)                            stop('Error: Size parameter should be a non-negative integer')
  if (n < 0)                                 stop('Error: Input n should be a single non-negative integer')
  if (prob < 0)                              stop('Error: Probability parameter should be between zero and one')
  if (prob > 1)                              stop('Error: Probability parameter should be between zero and one')
  
  ########################################################################################
  ################################# TRIVIAL CASES ########################################
  ########################################################################################
  
  #Deal with trivial case where n = 0
  if (n == 0) { return(numeric(0)) }
  
  #Deal with trivial case where size = 0
  #Distribution is a point-mass on zero
  if (nn == 0) {
    OUT <- rep(0, n)
    return(OUT) }
  
  #Deal with trivial case where n = 1
  #Distribution is a point-mass on m
  if (nn == 1) {
    OUT <- rep(m, n)
    return(OUT) }
  
  #Deal with special case where prob = 1
  #Distribution is a point-mass on size*m
  if (prob == 1) {
    OUT <- rep(nn*m, n)
    return(OUT) }
  
  #Deal with case where size = Inf and prob = 0
  #Distribution is Poisson with rate m
  if ((nn == Inf)&(prob == 0)) { 
    OUT <- rpois(n, lambda = m)
    return(OUT) }
  
  #Deal with case where size = Inf and prob > 0
  #Distribution is a point-mass on infinity
  if ((nn == Inf)&(prob > 0)) { 
    OUT <- rep(Inf, n)
    return(OUT) }
  
  ########################################################################################
  ############################### NON-TRIVIAL CASES ######################################
  ########################################################################################
  
  #Generate output vector
  RAND  <- matrix(0, nrow = m, ncol = n)
  for (t in 1:m) {
    KNOWN <- rbinom(n, size = size, prob = prob)
    for (i in 1:n) {
      nnn  <- nn-KNOWN[i]
      PERM <- sample.int(nnn, size = nnn, replace = FALSE)
      RAND[t, i] <- KNOWN[i] + sum(PERM == 1:nnn) } }
  OUT <- colSums(RAND)
  
  #Return output
  OUT }
