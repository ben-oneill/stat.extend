#' Matching test
#'
#' \code{matching.test} performs the matching test.
#'
#' The matching test considers a situation where a person attempts to match a set of objects to a corresponding set of
#' positions.  The null hypothesis is that the matching is performed at random and the alternative hypothesis is that the
#' matching occurs with some ability on the part of the player.
#' 
#' For data vectors \code{x} that are not too large we perform an exact test by computing the exact distribution of the 
#' mean number of matches.  If the number of data points is too large then we perform an approximate test using the normal
#' approximation to the distribution of the mean number of matches.  The parameter \code{max.m} sets the maximum number of 
#' data points where we perform an exact test.
#'
#' @param x Sample vector containing values from the generalised matching distribution
#' @param size The size parameter (number of objects to match)
#' @param null.prob The null value of the probability parameter
#' @param alternative The alternative hypothesis
#' @param approx A logical value specifying whether to use the normal approximation to the distribution
#' @return An htest object giving the output of the matching test

matching.test <- function(x, size, null.prob = 0, alternative = 'greater', approx = (length(x) > 100)) {
  
  #Set description of test and data
  method    <- 'Matching test'
  data.name <- deparse(substitute(x))
  
  #Check inputs x and size
  if (!is.numeric(size))                      stop('Error: Size parameter should be a non-negative integer') 
  if (length(size) != 1)                      stop('Error: Size parameter should be a single non-negative integer') 
  n <- as.integer(size)
  if (!(n == size))                           stop('Error: Size parameter should be a non-negative integer') 
  if (size < 0)                               stop('Error: Size parameter should be a non-negative integer') 
  if (!is.numeric(x))                         stop('Error: Input x should be numeric') 
  xx <- as.integer(x)
  m <- length(x)
  for (i in 1:m) {
    if (!(x[i] == xx[i]))                     stop(paste0('Error: Element ', i, ' in input x is not an integer')) 
    if (xx[i] < 0)                            stop(paste0('Error: Element ', i, ' in input x is negative')) 
    if (xx[i] > n)                            stop(paste0('Error: Element ', i, ' in input x is greater than size parameter')) 
    if (xx[i] == n-1)                         stop(paste0('Error: Element ', i, ' in input x is not in support of distribution')) }
  
  #Check input approx
  if (!is.logical(approx))                  stop('Error: Input approx should be a single logical value')
  if (length(approx) != 1)                  stop('Error: Input approx should be a single logical value')
  if ((approx)&(m < 30)) {
    warning(paste0('Using approximate computation when number of data points in ', data.name, 
                   ' is low can lead to poor approximation of p-value --- \n',
                   '    we recommend you use an exact test')) }
  
  #Check inputs null and alternative
  if (!is.numeric(null.prob))                 stop('Error: Input null.prob should be a probability value') 
  if (length(null.prob) != 1)                 stop('Error: Input null.prob should be a single probability value') 
  if (min(null.prob) < 0)                     stop('Error: Input null.prob should be a single probability value') 
  if (max(null.prob) > 1)                     stop('Error: Input null.prob should be a single probability value') 
  ALT.VALS <- c('two.sided', 'less', 'greater')
  if (!(alternative %in% ALT.VALS))                 stop('Error: Alternative hypothesis not recognised')
  if ((alternative == 'less')&(null.prob == 0))     stop('Error: Alternative hypothesis is empty')
  if ((alternative == 'greater')&(null.prob == 1))  stop('Error: Alternative hypothesis is empty')
  
  #Calculate test statistic and set test details
  S <- sum(x)
  null.value  <- null.prob
  names(null.value) <- 'prob'
  statistic  <- mean(x)
  attr(statistic, 'names') <- 'mean matches'
  
  #Compute p-value for trivial case where n <= 1
  if (n <= 1) { 
    p.value <- 1
    warning('For size <= 1 the matching test is trivial --- p-value is always one and power is zero') }
  
  #Compute p-value for non-trivial case where n > 1
  if (n > 1) {
  if (approx) {
    
    #Compute normal approximation to generalised matching distribution
    MEAN <- 1 + n*null.prob - null.prob^n
    VAR  <- 1 - null.prob^(2*n) + 
            n*(null.prob - null.prob^2 - null.prob^(n-1) - null.prob^n + 2*null.prob^(n+1))
    LOGPROBS.TOTAL <- dnorm(0:(n*m), mean = m*MEAN, sd = sqrt(m*VAR), log = TRUE)
    LOGPROBS.TOTAL <- LOGPROBS.TOTAL - matrixStats::logSumExp(LOGPROBS.TOTAL)
    
    #Compute p-value
    if (alternative == 'greater') {
      p.value  <- exp(matrixStats::logSumExp(LOGPROBS.TOTAL[S:(n*m)+1])) }
    if (alternative == 'less') {
      p.value  <- exp(matrixStats::logSumExp(LOGPROBS.TOTAL[0:S+1])) }
    if (alternative == 'two.sided') {
      TTT <- (LOGPROBS.TOTAL >= LOGPROBS.TOTAL[S+1])
      p.value  <- exp(matrixStats::logSumExp(LOGPROBS.TOTAL[TTT])) }
  
  } else {
    
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
    LOGBINOM <- dbinom(0:n, size = n, prob = null.prob, log = TRUE)
    for (k in 0:n) {
      LOGPROBS[k+1] <- matrixStats::logSumExp(LOGBINOM + M[k+1, ]) }
    LOGPROBS <- LOGPROBS - matrixStats::logSumExp(LOGPROBS)

    #Compute exact distribution of sample total
    LOGPROBS.TOTAL <- matrix(-Inf, nrow = m, ncol = n*m+1)
    LOGPROBS.TOTAL[1, 0:n+1] <- LOGPROBS
    if (m > 1) {
      for (t in 2:m) {
        for (s in 0:(n*t)) {
          k  <- 0:min(s, n)
          s0 <- s - k
          T1 <- LOGPROBS.TOTAL[t-1, s0+1]
          T2 <- LOGPROBS[k+1]
          LOGPROBS.TOTAL[t, s+1] <- matrixStats::logSumExp(T1 + T2) }
        LOGPROBS.TOTAL[t, ] <- LOGPROBS.TOTAL[t, ] - matrixStats::logSumExp(LOGPROBS.TOTAL[t, ]) } }
    
    #Compute p-value
    if (alternative == 'greater') {
      p.value  <- exp(matrixStats::logSumExp(LOGPROBS.TOTAL[m, S:(n*m)+1])) }
    if (alternative == 'less') {
      p.value  <- exp(matrixStats::logSumExp(LOGPROBS.TOTAL[m, 0:S+1])) }
    if (alternative == 'two.sided') {
      TTT <- (LOGPROBS.TOTAL[m, ] >= LOGPROBS.TOTAL[m, S+1])
      p.value  <- exp(matrixStats::logSumExp(LOGPROBS.TOTAL[TTT])) } } }
  
  #Create htest object
  TEST <- list(method = method, data.name = data.name, data = x, size = size, approx = approx,
               null.value = null.value, alternative = alternative, 
               statistic = statistic, p.value = p.value)
  class(TEST) <- 'htest'
  
  #Return htest
  TEST }
